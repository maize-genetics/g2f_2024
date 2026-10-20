#!/usr/bin/env python3
"""Structural summary of AnchorWave proali anchors (and optionally its MAF) for
one contig assembly against a reference: collinear blocks, contig assignment,
junctions between blocks along each contig, reference loci matched by more than
one assembly copy, cross-referenced with read-depth and pairwise evidence from
g2f_asmqc when available, plus dot plots.

  g2f_anchorwave_synteny.py --anchors A [--variant NAME=ANCHORS ...]
        --query-fai Q.fai --ref-fai R.fai --out DIR
        [--maf ALIGNMENT.maf] [--asmqc-assembly-dir DIR] [--asmqc-pair-dir DIR]
        [--query-fasta Q.fa --self-align-minimap2 /path/minimap2]
        [--zoom chr2:80000000-130000000 ...] [--label TEXT]

Every output is an observation from gene anchors. Junction types describe how
adjacent blocks relate; they are not misassembly calls. Copy classes combine
heuristics (documented in the TSV header) and are candidates only.
"""
import argparse
import bisect
import gzip
import json
import os
import statistics
import subprocess
import sys
import tempfile

# ------------------------------------------------------------------ thresholds
MIN_BLOCK_ANCHORS = 3          # blocks with fewer gene anchors are reported but marked weak
JUMP_BP = 1_000_000            # same chromosome/strand: |reference step - query step| above this -> length_discordant_gap
REVERSAL_BP = 50_000           # same chromosome/strand: reference moving backwards by more than this -> order reversal
HALF = (0.30, 0.70)            # depth_ratio (depth / c) of one separately represented allele
NORMAL = (0.70, 1.40)          # depth_ratio of sequence represented once with both haplotypes' reads
ALLELE_MIN_IDENTITY = 0.98     # copy-vs-copy identity for "alternative alleles"


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_fai(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def read_anchors(path):
    """Gene anchors only (interanchor/local alignment rows dropped)."""
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            r = dict(zip(header, f))
            if r["gene"] == "interanchor" or r["gene"].startswith("localAlignment"):
                continue
            rows.append({"ref": r["refChr"], "rs": int(r["referenceStart"]), "re": int(r["referenceEnd"]),
                         "q": r["queryChr"], "qs": int(r["queryStart"]), "qe": int(r["queryEnd"]),
                         "strand": r["strand"], "gene": r["gene"], "block": r["blockIndex"]})
    return rows


def write_tsv(path, header, rows, comment=None):
    with open(path, "w") as fh:
        if comment:
            for c in comment.strip().splitlines():
                fh.write(f"# {c}\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(h) is None else (f"{r[h]:.4g}" if isinstance(r[h], float) else str(r[h]))
                               for h in header) + "\n")


# ------------------------------------------------------------------ blocks

def build_blocks(anchors):
    by = {}
    for a in anchors:
        by.setdefault(a["block"], []).append(a)
    blocks = []
    for bid, rs in by.items():
        rs.sort(key=lambda a: a["qs"])
        b = {"block": bid, "ref": rs[0]["ref"], "q": rs[0]["q"], "strand": rs[0]["strand"], "n": len(rs),
             "rs": min(a["rs"] for a in rs), "re": max(a["re"] for a in rs),
             "qs": min(a["qs"] for a in rs), "qe": max(a["qe"] for a in rs), "genes": [a["gene"] for a in rs],
             "anchors": rs}
        if len({a["ref"] for a in rs}) != 1 or len({a["q"] for a in rs}) != 1:
            raise SystemExit(f"block {bid} spans several sequences; unexpected anchors format")
        b["weak"] = b["n"] < MIN_BLOCK_ANCHORS
        blocks.append(b)
    return blocks


def contig_assignment(blocks, qlen):
    per = {}
    for b in blocks:
        d = per.setdefault(b["q"], {"anchors": 0, "chroms": {}, "plus": 0, "minus": 0, "blocks": 0, "span_bp": 0})
        d["anchors"] += b["n"]
        d["blocks"] += 1
        d["chroms"][b["ref"]] = d["chroms"].get(b["ref"], 0) + b["n"]
        d["plus" if b["strand"] == "+" else "minus"] += b["n"]
        d["span_bp"] += b["qe"] - b["qs"] + 1
    rows = []
    for q, L in sorted(qlen.items(), key=lambda kv: -kv[1]):
        d = per.get(q)
        if not d:
            rows.append({"contig": q, "length": L, "anchors": 0, "blocks": 0, "major_chrom": None})
            continue
        major = max(d["chroms"], key=d["chroms"].get)
        strong = sorted(c for c, n in d["chroms"].items() if n >= MIN_BLOCK_ANCHORS)
        rows.append({"contig": q, "length": L, "anchors": d["anchors"], "blocks": d["blocks"],
                     "major_chrom": major, "major_frac": d["chroms"][major] / d["anchors"],
                     "chroms_ge3_anchors": ",".join(strong), "n_chroms_ge3": len(strong),
                     "minus_frac": d["minus"] / d["anchors"], "block_span_frac": min(1.0, d["span_bp"] / L)})
    return rows


def junctions(blocks):
    """Adjacent blocks along each contig (ordered by query start)."""
    by = {}
    for b in blocks:
        by.setdefault(b["q"], []).append(b)
    out = []
    for q, bs in by.items():
        bs.sort(key=lambda b: b["qs"])
        for x, y in zip(bs, bs[1:]):
            # reference coordinate at the facing ends of the two blocks
            xa = x["anchors"][-1]
            ya = y["anchors"][0]
            if x["ref"] != y["ref"]:
                kind = "inter_chromosomal"
            elif x["strand"] != y["strand"]:
                kind = "strand_switch"
            else:
                # collinear continuation expects the reference to keep moving in the block direction
                step = (ya["rs"] - xa["re"]) if x["strand"] == "+" else (xa["rs"] - ya["re"])
                if step < -REVERSAL_BP:
                    kind = "intra_order_reversal"
                elif abs(step - (ya["qs"] - xa["qe"])) > JUMP_BP:
                    kind = "length_discordant_gap"
                else:
                    kind = "collinear_gap"
            overlap = y["qs"] <= x["qe"]
            out.append({"contig": q, "q_break_start": xa["qe"], "q_break_end": ya["qs"],
                        "query_interval_bp": ya["qs"] - xa["qe"], "nested_or_overlapping": overlap,
                        "left_block": x["block"], "left_ref": f"{x['ref']}:{x['rs']}-{x['re']}", "left_strand": x["strand"],
                        "left_anchors": x["n"], "left_ref_at_break": f"{xa['ref']}:{xa['rs'] if x['strand'] == '-' else xa['re']}",
                        "right_block": y["block"], "right_ref": f"{y['ref']}:{y['rs']}-{y['re']}", "right_strand": y["strand"],
                        "right_anchors": y["n"], "right_ref_at_break": f"{ya['ref']}:{ya['re'] if y['strand'] == '-' else ya['rs']}",
                        "type": kind, "min_flank_anchors": min(x["n"], y["n"]),
                        "ref_step_bp": ((ya["rs"] - xa["re"]) if x["strand"] == "+" else (xa["rs"] - ya["re"]))
                        if x["ref"] == y["ref"] and x["strand"] == y["strand"] else None})
    return out


# ------------------------------------------------------------------ reference copy depth

def gene_copies(anchors):
    d = {}
    for a in anchors:
        d.setdefault(a["gene"], []).append(a)
    return d


def ref_windows(anchors, rlen, win=1_000_000):
    rows = {}
    for g, cs in gene_copies(anchors).items():
        a = cs[0]
        key = (a["ref"], a["rs"] // win)
        r = rows.setdefault(key, {"chrom": a["ref"], "start": key[1] * win, "end": min(rlen.get(a["ref"], 0), (key[1] + 1) * win),
                                  "genes": 0, "c1": 0, "c2": 0, "c3plus": 0})
        r["genes"] += 1
        r["c1" if len(cs) == 1 else ("c2" if len(cs) == 2 else "c3plus")] += 1
    return [rows[k] for k in sorted(rows, key=lambda k: (chrom_key(k[0]), k[1]))]


def chrom_key(c):
    import re
    m = re.match(r"chr(\d+)$", c)
    return (0, int(m.group(1)), c) if m else (1, 0, c)


def duplicated_regions(blocks):
    """Pairs of blocks anchored to the same reference genes = reference loci with two assembly copies."""
    gene_block = {}
    for b in blocks:
        for a in b["anchors"]:
            gene_block.setdefault(a["gene"], []).append((b, a))
    pairs = {}
    for g, lst in gene_block.items():
        if len(lst) < 2:
            continue
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                (b1, a1), (b2, a2) = sorted([lst[i], lst[j]], key=lambda t: (t[0]["q"], t[0]["qs"]))
                k = (b1["block"], b2["block"])
                p = pairs.setdefault(k, {"b1": b1, "b2": b2, "a1": [], "a2": [], "genes": 0})
                p["a1"].append(a1)
                p["a2"].append(a2)
                p["genes"] += 1
    out = []
    for (k1, k2), p in pairs.items():
        b1, b2 = p["b1"], p["b2"]
        rs = min(a["rs"] for a in p["a1"])
        re_ = max(a["re"] for a in p["a1"])
        q1s, q1e = min(a["qs"] for a in p["a1"]), max(a["qe"] for a in p["a1"])
        q2s, q2e = min(a["qs"] for a in p["a2"]), max(a["qe"] for a in p["a2"])
        out.append({"ref_chrom": b1["ref"], "ref_start": rs, "ref_end": re_, "shared_genes": p["genes"],
                    "copy1": f"{b1['q']}:{q1s}-{q1e}", "copy1_contig": b1["q"], "copy1_start": q1s, "copy1_end": q1e,
                    "copy1_strand": b1["strand"], "copy1_block": k1, "copy1_block_anchors": b1["n"],
                    "copy2": f"{b2['q']}:{q2s}-{q2e}", "copy2_contig": b2["q"], "copy2_start": q2s, "copy2_end": q2e,
                    "copy2_strand": b2["strand"], "copy2_block": k2, "copy2_block_anchors": b2["n"],
                    "same_contig": b1["q"] == b2["q"]})
    out.sort(key=lambda r: (chrom_key(r["ref_chrom"]), r["ref_start"]))
    return out


# ------------------------------------------------------------------ evidence from g2f_asmqc

class Windows:
    """Interval lookup over 10-kb window tables (contig, start, end, ...)."""

    def __init__(self, path, fields):
        self.d = {}
        with open_text(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rec = {k: f[idx[k]] for k in fields if k in idx}
                self.d.setdefault(f[0], ([], []))
                self.d[f[0]][0].append(int(f[1]))
                self.d[f[0]][1].append((int(f[2]), rec))

    def over(self, contig, s, e):
        if contig not in self.d:
            return []
        starts, recs = self.d[contig]
        i = max(0, bisect.bisect_right(starts, s) - 1)
        out = []
        while i < len(starts) and starts[i] < e:
            if recs[i][0] > s:
                out.append(recs[i][1])
            i += 1
        return out


def depth_summary(win, contig, s, e):
    ws = win.over(contig, s, e) if win else []
    if not ws:
        return {}
    r = [float(w["depth_ratio"]) for w in ws if w.get("depth_ratio") not in (None, "", "nan")]
    cls = {}
    for w in ws:
        cls[w["class"]] = cls.get(w["class"], 0) + 1
    top = max(cls, key=cls.get)
    return {"depth_ratio_median": statistics.median(r) if r else None, "depth_class": top,
            "depth_class_frac": cls[top] / len(ws),
            "multimap_median": statistics.median(float(w["multimap_frac"]) for w in ws)}


class PairClasses:
    """-l0 vs -l3 classes per -l0 base (redundant_copies + a_unaligned)."""

    def __init__(self, pair_dir):
        self.iv, self.ident = {}, {}
        p = os.path.join(pair_dir, "redundant_copies.tsv")
        if os.path.exists(p):
            with open(p) as fh:
                h = fh.readline().rstrip("\n").split("\t")
                ix = {k: i for i, k in enumerate(h)}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    self.ident.setdefault(f[ix["a_contig"]], []).append(
                        (int(f[ix["a_start"]]), int(f[ix["a_end"]]), float(f[ix["identity"]])))
            for v in self.ident.values():
                v.sort()
        for fn, cls_col in (("redundant_copies.tsv", "class"), ("a_unaligned.tsv", "class")):
            p = os.path.join(pair_dir, fn)
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                h = fh.readline().rstrip("\n").split("\t")
                ix = {k: i for i, k in enumerate(h)}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    c = f[ix.get("a_contig", ix.get("contig", 0))]
                    s = int(f[ix.get("a_start", ix.get("start"))])
                    e = int(f[ix.get("a_end", ix.get("end"))])
                    self.iv.setdefault(c, []).append((s, e, f[ix[cls_col]]))
        for v in self.iv.values():
            v.sort()

    def fractions(self, contig, s, e):
        tot = e - s
        out = {}
        for a, b, c in self.iv.get(contig, []):
            if b <= s:
                continue
            if a >= e:
                break
            out[c] = out.get(c, 0) + min(b, e) - max(a, s)
        cov = sum(out.values())
        out["shared"] = max(0, tot - cov)  # bases not listed are shared with the -l3 primary
        return {k: v / tot for k, v in out.items()} if tot > 0 else {}

    def counterpart_identity(self, contig, s, e):
        """bp-weighted identity of this interval's redundant parts to their counterpart in the other assembly."""
        w = num = 0
        for a, b, i in self.ident.get(contig, []):
            if b <= s:
                continue
            if a >= e:
                break
            n = min(b, e) - max(a, s)
            w += n
            num += n * i
        return num / w if w else None


def self_identity(minimap2, samtools, fasta, r):
    """Align copy2 against copy1 (asm20); returns aligned fraction of the shorter copy and identity."""
    def region(c, s, e):
        return f"{c}:{s}-{e}"
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as td:
        f1, f2 = os.path.join(td, "1.fa"), os.path.join(td, "2.fa")
        for fn, reg in ((f1, region(r["copy1_contig"], r["copy1_start"], r["copy1_end"])),
                        (f2, region(r["copy2_contig"], r["copy2_start"], r["copy2_end"]))):
            with open(fn, "w") as fh:
                subprocess.run([samtools, "faidx", fasta, reg], stdout=fh, check=True)
        p = subprocess.run([minimap2, "-x", "asm20", "-c", "--secondary=no", "-t", "2", f1, f2],
                           capture_output=True, text=True, check=True)
    m = blen = q_al = 0
    for line in p.stdout.splitlines():
        f = line.split("\t")
        m += int(f[9])
        blen += int(f[10])
        q_al += int(f[3]) - int(f[2])
    shorter = min(r["copy1_end"] - r["copy1_start"], r["copy2_end"] - r["copy2_start"]) + 1
    return {"copy_aligned_frac": round(min(1.0, q_al / shorter), 3) if shorter else None,
            "copy_identity": round(m / blen, 4) if blen else None}


def classify_copies(r):
    d1, d2 = r.get("copy1_depth_ratio"), r.get("copy2_depth_ratio")
    ident = r.get("copy_identity")
    if ident is None:
        known = [v for v in (r.get("copy1_l3_counterpart_identity"), r.get("copy2_l3_counterpart_identity")) if v]
        ident = max(known) if known else None
    if d1 is None or d2 is None:
        return "unresolved", "no read-depth evidence"
    half = lambda d: HALF[0] <= d < HALF[1]
    normal = lambda d: NORMAL[0] <= d < NORMAL[1]
    if half(d1) and half(d2):
        if ident is not None and ident < ALLELE_MIN_IDENTITY:
            return "unresolved", f"both copies at one-allele depth but copy identity {ident:.3f}"
        return "candidate_alternative_alleles", "both copies at one-allele depth" + (
            f", copy identity {ident:.3f}" if ident is not None else ", copy identity not measured")
    if normal(d1) and normal(d2):
        return "candidate_genuine_duplication", "both copies at full (two-haplotype) depth"
    if (half(d1) and normal(d2)) or (normal(d1) and half(d2)):
        return "unresolved", "one copy at one-allele depth, the other at full depth"
    return "unresolved", f"depth ratios {d1:.2f} / {d2:.2f} outside the half/normal ranges"


# ------------------------------------------------------------------ MAF

def maf_coverage(path, rlen):
    """Reference bases covered by >=1 and >=2 alignments, query bases aligned per contig, identity."""
    import numpy as np
    ref_iv, q_bp, matches, cols = {}, {}, 0, 0
    with open(path) as fh:
        block = []
        for line in fh:
            if line.startswith("s"):
                block.append(line.split())
            elif not line.strip() and block:
                r, q = block[0], block[1]
                ref_iv.setdefault(r[1], []).append((int(r[2]), int(r[2]) + int(r[3])))
                q_bp[q[1]] = q_bp.get(q[1], 0) + int(q[3])
                a = np.frombuffer(r[6].upper().encode(), dtype=np.uint8)
                b = np.frombuffer(q[6].upper().encode(), dtype=np.uint8)
                both = (a != 45) & (b != 45)
                matches += int(np.count_nonzero((a == b) & both))
                cols += int(np.count_nonzero(both))
                block = []
    per_chrom = []
    for c, ivs in ref_iv.items():
        ev = []
        for s, e in ivs:
            ev += [(s, 1), (e, -1)]
        ev.sort()
        depth, last, c1, c2 = 0, None, 0, 0
        for pos, d in ev:
            if last is not None and depth >= 1:
                c1 += pos - last
                if depth >= 2:
                    c2 += pos - last
            depth += d
            last = pos
        per_chrom.append({"chrom": c, "length": rlen.get(c), "covered_ge1": c1, "covered_ge2": c2,
                          "frac_ge1": c1 / rlen[c] if rlen.get(c) else None,
                          "frac_ge2": c2 / rlen[c] if rlen.get(c) else None})
    per_chrom.sort(key=lambda r: chrom_key(r["chrom"]))
    return per_chrom, q_bp, (matches / cols if cols else None)


# ------------------------------------------------------------------ plots

def plots(anchors, rlen, qlen, assign, out, label, zooms):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chroms = sorted([c for c in rlen if c.startswith("chr")], key=chrom_key)
    roff, o = {}, 0
    for c in chroms:
        roff[c] = o
        o += rlen[c]
    # contigs ordered by major chromosome, then median reference position
    med = {}
    for a in anchors:
        med.setdefault(a["q"], []).append((a["ref"], a["rs"]))
    major = {r["contig"]: r["major_chrom"] for r in assign if r.get("major_chrom")}
    def ckey(q):
        m = major.get(q)
        pos = statistics.median(p for c, p in med[q] if c == m) if m and q in med else 0
        return (chrom_key(m) if m else (2, 0, ""), pos)
    order = sorted([q for q in qlen if q in med], key=ckey)
    qoff, o = {}, 0
    for q in order:
        qoff[q] = o
        o += qlen[q]
    fig, ax = plt.subplots(figsize=(11, 10))
    for strand, col in (("+", "#2a6fbb"), ("-", "#d1495b")):
        xs = [qoff[a["q"]] + a["qs"] for a in anchors if a["strand"] == strand and a["ref"] in roff and a["q"] in qoff]
        ys = [roff[a["ref"]] + a["rs"] for a in anchors if a["strand"] == strand and a["ref"] in roff and a["q"] in qoff]
        ax.scatter([x / 1e6 for x in xs], [y / 1e6 for y in ys], s=0.3, c=col, label=f"strand {strand}", rasterized=True)
    for c in chroms:
        ax.axhline(roff[c] / 1e6, lw=0.3, c="0.7")
    ax.set_yticks([(roff[c] + rlen[c] / 2) / 1e6 for c in chroms])
    ax.set_yticklabels(chroms)
    ax.set_xlabel(f"{label} contigs (ordered by assigned chromosome), Mb")
    ax.set_ylabel("reference, Mb")
    ax.set_title(f"{label}: AnchorWave gene anchors")
    ax.legend(markerscale=20, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"dotplot_genome.{ext}"), dpi=150)
    plt.close(fig)
    # one panel per chromosome: reference position vs contig position, contigs with >= 3 anchors on it
    fig, axes = plt.subplots(5, 2, figsize=(12, 22))
    for ax, c in zip(axes.flat, chroms[:10]):
        sub = [a for a in anchors if a["ref"] == c]
        cnt = {}
        for a in sub:
            cnt[a["q"]] = cnt.get(a["q"], 0) + 1
        keep = [q for q in order if cnt.get(q, 0) >= MIN_BLOCK_ANCHORS]
        lo, o = {}, 0
        for q in keep:
            lo[q] = o
            o += qlen[q]
        for strand, col in (("+", "#2a6fbb"), ("-", "#d1495b")):
            pts = [(lo[a["q"]] + a["qs"], a["rs"]) for a in sub if a["q"] in lo and a["strand"] == strand]
            ax.scatter([p[0] / 1e6 for p in pts], [p[1] / 1e6 for p in pts], s=0.5, c=col, rasterized=True)
        for q in keep:
            ax.axvline(lo[q] / 1e6, lw=0.2, c="0.8")
        ax.set_title(f"{c}: {len(keep)} contigs with >= {MIN_BLOCK_ANCHORS} anchors", fontsize=9)
        ax.set_ylabel(f"{c}, Mb", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"dotplot_chromosomes.{ext}"), dpi=130)
    plt.close(fig)
    for z in zooms:
        c, _, rng = z.partition(":")
        s, e = (int(x) for x in rng.split("-"))
        sub = [a for a in anchors if a["ref"] == c and s <= a["rs"] <= e]
        cnt = {}
        for a in sub:
            cnt[a["q"]] = cnt.get(a["q"], 0) + 1
        keep = sorted([q for q in cnt if cnt[q] >= MIN_BLOCK_ANCHORS], key=lambda q: statistics.median(a["rs"] for a in sub if a["q"] == q))
        if not keep:
            continue
        fig, axes = plt.subplots(1, len(keep), figsize=(3.2 * len(keep), 5), squeeze=False)
        for ax, q in zip(axes[0], keep):
            for strand, col in (("+", "#2a6fbb"), ("-", "#d1495b")):
                pts = [(a["qs"], a["rs"]) for a in sub if a["q"] == q and a["strand"] == strand]
                ax.scatter([p[0] / 1e6 for p in pts], [p[1] / 1e6 for p in pts], s=2, c=col)
            ax.set_title(f"{q} ({cnt[q]} anchors)", fontsize=8)
            ax.set_xlabel(f"{q}, Mb", fontsize=8)
            ax.set_ylim(s / 1e6, e / 1e6)
            ax.tick_params(labelsize=7)
        axes[0][0].set_ylabel(f"{c}, Mb")
        fig.suptitle(f"{label}: {z}", fontsize=10)
        fig.tight_layout()
        tag = z.replace(":", "_").replace("-", "_")
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(out, f"zoom_{tag}.{ext}"), dpi=150)
        plt.close(fig)


# ------------------------------------------------------------------ main

def variant_stats(anchors):
    gc = gene_copies(anchors)
    dist = {}
    for cs in gc.values():
        k = str(len(cs)) if len(cs) < 3 else "3+"
        dist[k] = dist.get(k, 0) + 1
    return {"gene_anchors": len(anchors), "distinct_genes": len(gc), "copies_per_gene": dict(sorted(dist.items())),
            "minus_strand_anchors": sum(1 for a in anchors if a["strand"] == "-"),
            "blocks": len({a["block"] for a in anchors})}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True, help="baseline anchors file")
    ap.add_argument("--variant", action="append", default=[], help="NAME=ANCHORS for comparison")
    ap.add_argument("--query-fai", required=True)
    ap.add_argument("--ref-fai", required=True)
    ap.add_argument("--maf")
    ap.add_argument("--asmqc-assembly-dir", help="g2f_asmqc samples/<S>/<A> directory (regions, read_support)")
    ap.add_argument("--asmqc-pair-dir", help="g2f_asmqc pairs/<A>_vs_<B>/analysis (this assembly = A)")
    ap.add_argument("--query-fasta")
    ap.add_argument("--self-align-minimap2", help="minimap2 binary: align the two copies of each duplicated region")
    ap.add_argument("--self-align-min-genes", type=int, default=5)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--zoom", action="append", default=[])
    ap.add_argument("--label", default="assembly")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; choose a new output directory (results are never overwritten)")
    os.makedirs(a.out)
    qlen, rlen = read_fai(a.query_fai), read_fai(a.ref_fai)
    anchors = read_anchors(a.anchors)
    blocks = build_blocks(anchors)
    summary = {"label": a.label, "anchors_file": os.path.abspath(a.anchors), "baseline": variant_stats(anchors),
               "thresholds": {"MIN_BLOCK_ANCHORS": MIN_BLOCK_ANCHORS, "JUMP_BP": JUMP_BP, "REVERSAL_BP": REVERSAL_BP, "HALF": HALF,
                              "NORMAL": NORMAL, "ALLELE_MIN_IDENTITY": ALLELE_MIN_IDENTITY}}
    summary["variants"] = {}
    for v in a.variant:
        name, _, path = v.partition("=")
        va = read_anchors(path)
        summary["variants"][name] = {"file": os.path.abspath(path), **variant_stats(va)}
        vb = build_blocks(va)
        write_tsv(os.path.join(a.out, f"junctions_{name}.tsv"),
                  ["contig", "q_break_start", "q_break_end", "query_interval_bp", "type", "left_ref", "left_strand",
                   "left_anchors", "left_ref_at_break", "right_ref", "right_strand", "right_anchors",
                   "right_ref_at_break", "min_flank_anchors", "ref_step_bp", "nested_or_overlapping"],
                  [j for j in junctions(vb) if j["type"] != "collinear_gap"])
    # blocks
    write_tsv(os.path.join(a.out, "blocks.tsv"),
              ["block", "q", "qs", "qe", "strand", "ref", "rs", "re", "n", "weak"],
              sorted(blocks, key=lambda b: (b["q"], b["qs"])))
    assign = contig_assignment(blocks, qlen)
    write_tsv(os.path.join(a.out, "contig_assignment.tsv"),
              ["contig", "length", "anchors", "blocks", "major_chrom", "major_frac", "n_chroms_ge3",
               "chroms_ge3_anchors", "minus_frac", "block_span_frac"], assign)
    total = sum(qlen.values())
    with_anchor = sum(r["length"] for r in assign if r["anchors"])
    multi = [r for r in assign if (r.get("n_chroms_ge3") or 0) > 1]
    summary["contigs"] = {"n": len(qlen), "bp": total, "with_anchors": sum(1 for r in assign if r["anchors"]),
                          "bp_in_contigs_with_anchors": with_anchor,
                          "block_span_bp": sum(b["qe"] - b["qs"] + 1 for b in blocks if not b["weak"]),
                          "contigs_multi_chrom_ge3": len(multi), "bp_multi_chrom_ge3": sum(r["length"] for r in multi)}
    js = junctions(blocks)
    comment = ("Adjacent collinear blocks along each contig. type: inter_chromosomal | strand_switch | "
               f"intra_order_reversal (reference backwards > {REVERSAL_BP} bp) | length_discordant_gap "
               f"(|ref_step_bp - query_interval_bp| > {JUMP_BP}); collinear_gap rows are omitted.\n"
               "q_break_start/end: last gene anchor of the left block .. first of the right (query interval holding "
               "the junction). Observations relative to the reference, not misassembly calls.")
    jrows = [j for j in js if j["type"] != "collinear_gap"]
    # evidence at junction intervals
    win = rsw = None
    if a.asmqc_assembly_dir:
        p = os.path.join(a.asmqc_assembly_dir, "regions", "window_classes.tsv.gz")
        if os.path.exists(p):
            win = Windows(p, ["depth_ratio", "class", "multimap_frac", "mixed_per_kb"])
        p = os.path.join(a.asmqc_assembly_dir, "read_support", "windows.tsv.gz")
        if os.path.exists(p):
            rsw = Windows(p, ["depth_all_median", "span_min"])
    for j in jrows:
        s, e = sorted((j["q_break_start"], j["q_break_end"]))
        s, e = s - 5000, e + 5000
        if rsw:
            ws = rsw.over(j["contig"], s, e)
            if ws:
                j["hifi_min_spanning"] = min(int(float(w["span_min"])) for w in ws)
                j["hifi_median_depth"] = statistics.median(float(w["depth_all_median"]) for w in ws)
        if win:
            d = depth_summary(win, j["contig"], s, e)
            j["depth_class"] = d.get("depth_class")
            j["multimap_median"] = d.get("multimap_median")
    write_tsv(os.path.join(a.out, "junctions.tsv"),
              ["contig", "q_break_start", "q_break_end", "query_interval_bp", "type", "left_ref", "left_strand",
               "left_anchors", "left_ref_at_break", "right_ref", "right_strand", "right_anchors", "right_ref_at_break",
               "min_flank_anchors", "ref_step_bp", "nested_or_overlapping", "hifi_min_spanning", "hifi_median_depth", "depth_class",
               "multimap_median"], jrows, comment)
    summary["junctions"] = {}
    for j in jrows:
        k = j["type"] + ("_weak" if j["min_flank_anchors"] < MIN_BLOCK_ANCHORS else "")
        summary["junctions"][k] = summary["junctions"].get(k, 0) + 1
    write_tsv(os.path.join(a.out, "reference_windows.tsv"), ["chrom", "start", "end", "genes", "c1", "c2", "c3plus"],
              ref_windows(anchors, rlen))
    # duplicated syntenic regions
    dups = duplicated_regions(blocks)
    pc = PairClasses(a.asmqc_pair_dir) if a.asmqc_pair_dir else None
    for r in dups:
        for k in ("copy1", "copy2"):
            d = depth_summary(win, r[f"{k}_contig"], r[f"{k}_start"], r[f"{k}_end"])
            r[f"{k}_depth_ratio"] = d.get("depth_ratio_median")
            r[f"{k}_depth_class"] = d.get("depth_class")
            r[f"{k}_multimap"] = d.get("multimap_median")
            if pc:
                f = pc.fractions(r[f"{k}_contig"], r[f"{k}_start"], r[f"{k}_end"])
                r[f"{k}_l3_counterpart_identity"] = pc.counterpart_identity(r[f"{k}_contig"], r[f"{k}_start"], r[f"{k}_end"])
                r[f"{k}_l3_fate"] = ";".join(f"{c}={v:.2f}" for c, v in sorted(f.items(), key=lambda kv: -kv[1]) if v >= 0.05)
        if a.self_align_minimap2 and a.query_fasta and r["shared_genes"] >= a.self_align_min_genes:
            r.update(self_identity(a.self_align_minimap2, a.samtools, a.query_fasta, r))
        r["class"], r["class_reason"] = classify_copies(r)
    comment = ("Reference loci anchored by two assembly copies (pairs of collinear blocks sharing reference genes).\n"
               f"depth_ratio = median HiFi depth / c over the copy (g2f_asmqc regions); half {HALF}, normal {NORMAL}.\n"
               "class (heuristic, candidates only): candidate_alternative_alleles = both copies at one-allele depth "
               f"(and copy identity >= {ALLELE_MIN_IDENTITY} when measured: self-alignment, else the higher l3 counterpart identity); candidate_genuine_duplication = both at "
               "full depth; otherwise unresolved. Syntenic overlap alone does not establish allelic redundancy.\n"
               "l3_fate: fraction of the copy that is shared with the -l3 primary / redundant (counterpart in the "
               "primary) / in_B_other_outputs (hap1/hap2 only) / A_unique, from the -l0 vs -l3 pair analysis;\n"
               "l3_counterpart_identity: bp-weighted identity of the redundant part to its -l3 counterpart (for an allelic "
               "pair merged by purging this approximates copy-vs-copy identity).")
    write_tsv(os.path.join(a.out, "duplicated_regions.tsv"),
              ["ref_chrom", "ref_start", "ref_end", "shared_genes", "copy1", "copy1_strand", "copy1_block_anchors",
               "copy1_depth_ratio", "copy1_depth_class", "copy1_multimap", "copy1_l3_fate", "copy1_l3_counterpart_identity",
               "copy2", "copy2_strand", "copy2_block_anchors", "copy2_depth_ratio", "copy2_depth_class",
               "copy2_multimap", "copy2_l3_fate", "copy2_l3_counterpart_identity", "same_contig", "copy_aligned_frac", "copy_identity",
               "class", "class_reason"], dups, comment)
    cl = {}
    for r in dups:
        c = cl.setdefault(r["class"], {"pairs": 0, "shared_genes": 0, "ref_bp": 0})
        c["pairs"] += 1
        c["shared_genes"] += r["shared_genes"]
        c["ref_bp"] += r["ref_end"] - r["ref_start"] + 1
    summary["duplicated_regions"] = cl
    if a.maf:
        per_chrom, q_bp, ident = maf_coverage(a.maf, rlen)
        write_tsv(os.path.join(a.out, "maf_reference_coverage.tsv"),
                  ["chrom", "length", "covered_ge1", "covered_ge2", "frac_ge1", "frac_ge2"], per_chrom)
        write_tsv(os.path.join(a.out, "maf_query_aligned.tsv"), ["contig", "length", "aligned_bp", "aligned_frac"],
                  [{"contig": q, "length": L, "aligned_bp": q_bp.get(q, 0), "aligned_frac": q_bp.get(q, 0) / L}
                   for q, L in sorted(qlen.items(), key=lambda kv: -kv[1])])
        chrs = [r for r in per_chrom if r["chrom"].startswith("chr")]
        summary["maf"] = {"file": os.path.abspath(a.maf), "query_aligned_bp": sum(q_bp.values()),
                          "query_aligned_frac": sum(q_bp.values()) / total, "aligned_column_identity": ident,
                          "ref_chrom_covered_ge1": sum(r["covered_ge1"] for r in chrs),
                          "ref_chrom_covered_ge2": sum(r["covered_ge2"] for r in chrs),
                          "ref_chrom_bp": sum(r["length"] or 0 for r in chrs)}
    plots(anchors, rlen, qlen, assign, a.out, a.label, a.zoom)
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
