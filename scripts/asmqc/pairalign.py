"""Classify the sequence of assembly A (e.g. unpurged) against assembly B
(e.g. purged) of the same reads, from `minimap2 -cx asm5 B A` (PAF with
base-level alignment), optionally with A also aligned to B's other outputs
(e.g. hifiasm hap1 + hap2 contigs).

Every base of A falls in one class:
  shared                 the best A copy of a B region (B keeps this sequence)
  redundant_identical    a further A copy of a B region, identity >= near_identical
  redundant_divergent    a further A copy, identity < near_identical (allelic copy)
  in_B_other_outputs     not in B's primary but in B's other outputs (purged into haplotigs)
  A_unique               not found in any B output
B bases not covered by any A alignment are B_only.

If read-support window classes exist for both assemblies, redundant and
A_unique sequence is cross-tabulated by read-depth class (half / normal / dup):
an A copy at half depth whose B counterpart is at normal depth is a retained
allele that B correctly merged; a B counterpart at double depth means B
collapsed two real copies (over-purging).
"""
import gzip
import json
from bisect import bisect_right

from .paf import merge_intervals, read_paf, subtract


def best_copy_pieces(alns):
    """Split target coordinates into elementary pieces; per piece, the alignment
    with the highest identity (then longest) is 'best', the others redundant.
    Returns list of (aln_index, t_start, t_end, is_best)."""
    by_t = {}
    for i, a in enumerate(alns):
        by_t.setdefault(a["t"], []).append(i)
    pieces = []
    for t, idx in by_t.items():
        pts = sorted({alns[i]["ts"] for i in idx} | {alns[i]["te"] for i in idx})
        events = {}
        for i in idx:
            events.setdefault(alns[i]["ts"], []).append(("+", i))
            events.setdefault(alns[i]["te"], []).append(("-", i))
        active = set()
        for x, y in zip(pts, pts[1:]):
            for kind, i in events.get(x, []):
                (active.add if kind == "+" else active.discard)(i)
            if not active:
                continue
            best = max(active, key=lambda i: (alns[i]["identity"], alns[i]["alen"], -i))
            for i in active:
                pieces.append((i, x, y, i == best))
    return pieces


def target_to_query(a, ts, te):
    """Linear map of a target sub-interval to query coordinates."""
    tl = max(1, a["te"] - a["ts"])
    ql = a["qe"] - a["qs"]
    f0, f1 = (ts - a["ts"]) / tl, (te - a["ts"]) / tl
    if a["strand"] == "+":
        return int(a["qs"] + f0 * ql), int(a["qs"] + f1 * ql)
    return int(a["qe"] - f1 * ql), int(a["qe"] - f0 * ql)


class WindowClasses:
    """Read-depth class lookup from a regions-stage window_classes.tsv.gz."""

    def __init__(self, path):
        self.w = {}
        with gzip.open(path, "rt") as fh:
            fh.readline()
            for line in fh:
                c, s, e, ratio, cls = line.split("\t")[:5]
                d = self.w.setdefault(c, ([], [], []))
                d[0].append(int(s))
                d[1].append(int(e))
                d[2].append(cls)

    def classes(self, contig, s, e):
        """bp per class over [s, e)."""
        d = self.w.get(contig)
        out = {}
        if not d:
            return {"na": e - s}
        i = max(0, bisect_right(d[0], s) - 1)
        while i < len(d[0]) and d[0][i] < e:
            ov = min(e, d[1][i]) - max(s, d[0][i])
            if ov > 0:
                out[d[2][i]] = out.get(d[2][i], 0) + ov
            i += 1
        return out or {"na": e - s}

    def major(self, contig, s, e):
        c = self.classes(contig, s, e)
        return max(c, key=c.get)


def run(paf_path, a_lengths, b_lengths, sci, out_dir, extra_paf=None, a_classes=None, b_classes=None):
    p = sci["pair_align"]
    alns = list(read_paf(paf_path, min_len=p["min_aln_len"]))
    pieces = best_copy_pieces(alns)
    q_iv = {"shared": {}, "redundant_identical": {}, "redundant_divergent": {}}
    cross = {}
    redundant_rows = []
    ident_hist = {}
    for i, ts, te, best in pieces:
        a = alns[i]
        qs, qe = target_to_query(a, ts, te)
        if best:
            cls = "shared"
        else:
            cls = "redundant_identical" if a["identity"] >= p["near_identical"] else "redundant_divergent"
            redundant_rows.append((a["q"], qs, qe, a["t"], ts, te, a["strand"], round(a["identity"], 5), cls))
            b = round(min(a["identity"], 1.0) * 1000) / 1000
            ident_hist[b] = ident_hist.get(b, 0) + (te - ts)
            if a_classes and b_classes:
                key = f"{cls}|A:{a_classes.major(a['q'], qs, qe)}|B:{b_classes.major(a['t'], ts, te)}"
                cross[key] = cross.get(key, 0) + (te - ts)
        q_iv[cls].setdefault(a["q"], []).append((qs, qe))
    # A sequence covered by any alignment
    aligned = {}
    for a in alns:
        aligned.setdefault(a["q"], []).append((a["qs"], a["qe"]))
    extra = {}
    if extra_paf:
        for a in read_paf(extra_paf, min_len=p["min_aln_len"]):
            extra.setdefault(a["q"], []).append((a["qs"], a["qe"]))
    bp = {k: 0 for k in ("shared", "redundant_identical", "redundant_divergent", "in_B_other_outputs",
                         "A_unique")}
    unaligned_rows = []
    per_contig = []
    for q, qlen in a_lengths.items():
        al = merge_intervals(aligned.get(q, []))
        un = subtract([[0, qlen]], al)
        ex = merge_intervals(extra.get(q, []))
        in_other = [iv for iv in _intersect(un, ex)]
        uniq = subtract(un, ex) if ex else un
        row = {"contig": q, "length": qlen}
        # shared takes precedence over redundant where query pieces overlap
        sh = merge_intervals(q_iv["shared"].get(q, []))
        ri = subtract(merge_intervals(q_iv["redundant_identical"].get(q, [])), sh)
        rd = subtract(subtract(merge_intervals(q_iv["redundant_divergent"].get(q, [])), sh), ri)
        for k, iv in (("shared", sh), ("redundant_identical", ri), ("redundant_divergent", rd),
                      ("in_B_other_outputs", in_other), ("A_unique", uniq)):
            n = sum(e - s for s, e in iv)
            row[k] = n
            bp[k] += n
        for s, e in in_other:
            unaligned_rows.append((q, s, e, "in_B_other_outputs"))
        for s, e in uniq:
            unaligned_rows.append((q, s, e, "A_unique"))
            if a_classes:
                for c, n in a_classes.classes(q, s, e).items():
                    key = f"A_unique|A:{c}"
                    cross[key] = cross.get(key, 0) + n
        per_contig.append(row)
    b_cov = {}
    for a in alns:
        b_cov.setdefault(a["t"], []).append((a["ts"], a["te"]))
    b_only = sum(l - sum(e - s for s, e in merge_intervals(b_cov.get(t, []))) for t, l in b_lengths.items())

    with open(f"{out_dir}/redundant_copies.tsv", "w") as fh:
        fh.write("a_contig\ta_start\ta_end\tb_contig\tb_start\tb_end\tstrand\tidentity\tclass\n")
        for r in sorted(redundant_rows):
            fh.write("\t".join(map(str, r)) + "\n")
    with open(f"{out_dir}/a_unaligned.tsv", "w") as fh:
        fh.write("a_contig\tstart\tend\tclass\n")
        for r in unaligned_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    cols = ["contig", "length", "shared", "redundant_identical", "redundant_divergent", "in_B_other_outputs",
            "A_unique"]
    with open(f"{out_dir}/a_contig_classes.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(per_contig, key=lambda r: -r["length"]):
            fh.write("\t".join(str(r[k]) for k in cols) + "\n")
    with open(f"{out_dir}/redundant_identity_hist.tsv", "w") as fh:
        fh.write("identity\tbp\n")
        for k in sorted(ident_hist):
            fh.write(f"{k:.3f}\t{ident_hist[k]}\n")
    summary = {"a_bp": sum(a_lengths.values()), "b_bp": sum(b_lengths.values()), "class_bp": bp,
               "b_only_bp": b_only, "alignments": len(alns), "depth_crosstab_bp": cross,
               "extra_targets_used": bool(extra_paf)}
    with open(f"{out_dir}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def _intersect(a, b):
    out, j = [], 0
    for s, e in a:
        while j < len(b) and b[j][1] <= s:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            lo, hi = max(s, b[k][0]), min(e, b[k][1])
            if lo < hi:
                out.append([lo, hi])
            k += 1
    return out
