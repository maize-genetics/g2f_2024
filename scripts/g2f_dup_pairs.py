#!/usr/bin/env python3
"""Investigate same-FASTA duplicate candidates: two records of one assembly (haplotype) that carry the same
reference gene anchors (e.g. g2f_riparian `redundant_candidate`).

For each pair (A = smaller candidate, B = record covering it) the tool reports:
  anchoring   reference genes anchored on A, on B, on both (the same gene ID twice = repeated anchoring), and the
              coordinates of the shared-gene span on A, on B and on the reference
  alignment   minimap2 (asm20, -c) of A against the shared region of B, and, when another haplotype of the same line is
              given, of A and of the B region against that haplotype's counterpart region (C). Per alignment: identity
              (matches / alignment columns, alignments >= --min-aln-bp), gap-compressed identity (1 - de), and the
              fraction of the query and of the target covered
  depth       median read-depth ratio (observed / expected) and depth classes over A, over the B region and over the
              B region's flanks, from g2f_asmqc regions/window_classes.tsv[.gz] of the same FASTA

and a rule-based reading with a confidence:
  residual allelic duplication   A matches the other haplotype's allele (A~C identity >= --allelic-id and higher
                                 than A~B) and A~B are not near-identical: this haplotype carries both alleles
  near-identical duplicate       A~B >= --identical-id over most of A (both copies the same sequence: a collapsed /
                                 duplicated assembly copy or a very recent duplication; depth decides)
  paralogous anchoring           A and B share few genes (< --min-shared-frac of A's genes) or barely align: two
                                 diverged loci that AnchorWave anchored to the same reference genes
  true duplication candidate     A and B align over most of A at moderate identity, A does not match C better than B:
                                 a second, divergent copy in the same haplotype
  ambiguous                      anything else
The rules are explicit thresholds, not a decision; the table keeps every number so a reader can disagree.

Pairs table (TSV): pair_id, line, unit, fasta, anchors, record_a, record_b, and optional windows, expected_depth_note,
other_unit, other_fasta, other_anchors, other_records (comma-separated counterpart records in the other haplotype),
other_windows. FASTAs need a .fai; samtools and minimap2 are required for `run`.

Usage:
  g2f_dup_pairs.py plan --pairs pairs.tsv                 # what would be aligned (sizes), nothing run
  g2f_dup_pairs.py run  --pairs pairs.tsv --out OUTDIR [--threads 8] [--minimap2 PATH] [--samtools PATH]
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict

CODE_VERSION = 1
REQ = ("pair_id", "line", "unit", "fasta", "anchors", "record_a", "record_b")


def die(msg):
    sys.exit(f"g2f_dup_pairs: error: {msg}")


def read_tsv(path):
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


def read_fai(path):
    with open(path) as fh:
        return {f[0]: int(f[1]) for f in (ln.split("\t") for ln in fh)}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def anchors_by_record(path, records):
    """record -> {gene: (query_mid, ref_chr, ref_mid)} for gene anchors on the given records."""
    out = defaultdict(dict)
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("#", "refChr")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[7] == "interanchor" or f[3] not in records:
                continue
            out[f[3]][f[7]] = ((int(f[4]) + int(f[5])) // 2, f[0], (int(f[1]) + int(f[2])) // 2)
    return out


def span(genes, d):
    q = sorted(d[g][0] for g in genes)
    return (q[0], q[-1]) if q else (0, 0)


def region(center, L, flank):
    return max(1, center[0] - flank), min(L, center[1] + flank)


def read_windows(path):
    per = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("contig", "start", "end", "depth_ratio", "class")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            try:
                d = float(f[ix["depth_ratio"]])
            except ValueError:
                continue
            per[f[ix["contig"]]].append((int(f[ix["start"]]) + 1, int(f[ix["end"]]), d, f[ix["class"]]))
    return per


def depth(per, rec, s, e):
    sel = [w for w in per.get(rec, []) if w[1] >= s and w[0] <= e]
    if not sel:
        return dict(median="", classes="")
    ds = sorted(w[2] for w in sel)
    cls = defaultdict(int)
    for w in sel:
        cls[w[3]] += 1
    return dict(median=round(ds[len(ds) // 2], 3),
                classes=";".join(f"{k}:{v / len(sel):.2f}" for k, v in sorted(cls.items(), key=lambda kv: -kv[1])))


def union(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if ce is not None else 0)


def paf_stats(path, qlen, tlen, min_aln):
    m = blk = 0
    de_w = 0.0
    qiv, tiv = [], []
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tags = dict(t.split(":", 2)[0::2] for t in f[12:])
            if tags.get("tp") != "P" or int(f[10]) < min_aln:
                continue
            m += int(f[9])
            blk += int(f[10])
            if "de" in tags:
                de_w += float(tags["de"]) * int(f[10])
            qiv.append((int(f[2]), int(f[3])))
            tiv.append((int(f[7]), int(f[8])))
    if not blk:
        return dict(identity="", gc_identity="", q_cov=0.0, t_cov=0.0, aln_bp=0)
    return dict(identity=round(m / blk, 5), gc_identity=round(1 - de_w / blk, 5),
                q_cov=round(union(qiv) / qlen, 4), t_cov=round(union(tiv) / tlen, 4), aln_bp=blk)


COLS = ["pair_id", "line", "unit", "record_a", "len_a", "record_b", "len_b", "ref_span", "shared_span_a",
        "shared_span_b", "frac_a_in_shared_span", "frac_b_in_shared_span", "genes_a", "genes_b", "genes_shared",
        "shared_frac_of_a", "genes_a_only", "b_region", "c_region", "genes_shared_with_c"]
for _n in ("AB", "AC", "BC"):
    COLS += [f"{_n}_{k}" for k in ("identity", "gc_identity", "q_cov", "t_cov", "aln_bp")]
COLS += ["depth_a_median", "depth_a_classes", "depth_b_median", "depth_b_classes", "depth_b_flank_median",
         "depth_c_median", "AO_identity", "AO_gc_identity", "AO_q_cov", "AO_aln_bp", "AO_top_targets",
         "other_haplotype_state", "classification", "confidence", "basis"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def other_state(r, P):
    """What the other haplotype holds at this locus."""
    ac, acq, bc = num(r.get("AC_gc_identity")), num(r.get("AC_q_cov")), num(r.get("BC_gc_identity"))
    if ac is not None or bc is not None:
        if ac and ac >= P["allelic_id"] and (acq or 0) >= 0.5:
            return "carries A's allele"
        if bc and bc >= P["identical_id"]:
            return "carries a copy of B's allele (B~C near-identical)"
        return "counterpart diverges from both"
    aoq = num(r.get("AO_q_cov"))
    if aoq is not None:
        return (f"no anchored counterpart; A matches it over {aoq:.0%}" +
                (" (scattered, locus absent)" if aoq < 0.5 else " (locus present on unanchored records)"))
    return "not tested"


def classify(r, P):
    """Rule-based reading (thresholds in P; every number stays in the table)."""
    ab, abq = num(r.get("AB_gc_identity")), num(r.get("AB_q_cov")) or 0
    shared = num(r["shared_frac_of_a"])
    da, db, dfl = num(r.get("depth_a_median")), num(r.get("depth_b_median")), num(r.get("depth_b_flank_median"))
    state = other_state(r, P)
    r["other_haplotype_state"] = state
    if shared < P["min_shared_frac"] or ab is None or abq < 0.2:
        return ("paralogous anchoring (different loci, same reference genes)", "medium",
                f"A shares {shared:.0%} of its anchored genes with B; A~B aligns over {abq:.0%} of A")
    if state == "carries A's allele" and ab < P["identical_id"]:
        return ("residual allelic duplication: A is the other haplotype's allele, also kept here", "high",
                f"A~C {r['AC_gc_identity']} >= {P['allelic_id']}; A~B {ab:.4f}")
    split = (da is not None and db is not None and da <= P["half_max"] and db <= P["half_max"]
             and dfl is not None and dfl - db >= P["flank_gain"])
    if split and ab < P["identical_id"]:
        conf = "high" if (state.startswith("carries a copy of B") or "locus absent" in state) else "medium"
        return ("residual allelic duplication: both alleles of a heterozygous locus in this haplotype", conf,
                f"read depth splits between the copies (A {da}, B region {db}, B flanks {dfl}); A~B {ab:.4f} over "
                f"{abq:.0%} of A (divergent alleles, not identical copies); other haplotype: {state}")
    if ab >= P["identical_id"] and abq >= 0.8:
        return ("near-identical duplicate copy", "medium",
                f"A~B {ab:.4f} over {abq:.0%} of A; depth A {da}, B region {db}")
    if abq >= 0.5 and da is not None and da > P["half_max"] and (db or 0) > P["half_max"]:
        return ("true duplication candidate (both copies at full depth)", "low",
                f"A~B {ab:.4f} over {abq:.0%} of A; depth A {da}, B region {db}")
    return ("ambiguous", "low", f"A~B {ab:.4f} over {abq:.0%} of A; depth A {da}, B region {db}, B flanks {dfl}; "
                                f"other haplotype: {state}")


def reclassify(a, P):
    """Re-apply the rules to a finished run (no alignment); writes a new directory with the new table."""
    if not a.src or not a.out:
        die("reclassify needs --from and --out")
    if os.path.exists(a.out) or os.path.exists(a.out + ".partial"):
        die(f"{a.out} (or .partial) exists")
    rows = read_tsv(os.path.join(a.src, "dup_pairs.tsv"))
    for r in rows:
        r["classification"], r["confidence"], r["basis"] = classify(r, P)
    part = a.out + ".partial"
    os.makedirs(part)
    write_tsv(os.path.join(part, "dup_pairs.tsv"), rows, COLS)
    with open(os.path.join(part, "RUN.json"), "w") as fh:
        json.dump(dict(schema="g2f-dup-pairs/1", code_version=CODE_VERSION, command=[sys.executable] + sys.argv,
                       reclassified_from=os.path.abspath(a.src),
                       source_sha256=sha256(os.path.join(a.src, "dup_pairs.tsv")), params=P,
                       script_sha256=sha256(__file__),
                       start=datetime.datetime.now().astimezone().isoformat(timespec="seconds")), fh, indent=1)
    os.rename(part, a.out)
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["plan", "run", "reclassify"])
    ap.add_argument("--pairs", help="pairs table (plan / run)")
    ap.add_argument("--from", dest="src", help="reclassify: a finished run directory (re-reads dup_pairs.tsv)")
    ap.add_argument("--out")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--flank-bp", type=int, default=1_000_000, help="added around the shared-gene span on B and C")
    ap.add_argument("--min-aln-bp", type=int, default=5_000)
    ap.add_argument("--allelic-id", type=float, default=0.999)
    ap.add_argument("--identical-id", type=float, default=0.999)
    ap.add_argument("--min-shared-frac", type=float, default=0.3)
    ap.add_argument("--half-max", type=float, default=0.65, help="median depth ratio counted as half depth")
    ap.add_argument("--flank-gain", type=float, default=0.15,
                    help="B flanks must exceed the B shared region by this depth ratio for a depth split")
    ap.add_argument("--search-other", action="store_true",
                    help="when the other haplotype has no anchored counterpart, align A against that whole FASTA "
                         "(one minimap2 run per other FASTA, all such queries batched; ~RAM 15-25 GB for 2 Gb)")
    a = ap.parse_args(argv)
    P = dict(flank_bp=a.flank_bp, min_aln_bp=a.min_aln_bp, allelic_id=a.allelic_id, identical_id=a.identical_id,
             min_shared_frac=a.min_shared_frac, half_max=a.half_max, flank_gain=a.flank_gain)
    if a.command == "reclassify":
        return reclassify(a, P)
    if not a.pairs:
        die("plan / run need --pairs")
    pairs = read_tsv(a.pairs)
    for p in pairs:
        miss = [k for k in REQ if not p.get(k)]
        if miss:
            die(f"pair {p.get('pair_id', '?')} lacks {', '.join(miss)}")
        for k in ("fasta", "anchors", "windows", "other_fasta", "other_anchors", "other_windows"):
            if p.get(k) and not os.path.isfile(p[k]):
                die(f"{p['pair_id']}: {k} not found: {p[k]}")
        for k in ("fasta", "other_fasta"):
            if p.get(k) and not os.path.isfile(p[k] + ".fai"):
                die(f"{p['pair_id']}: {p[k]}.fai missing")
    if len({p["pair_id"] for p in pairs}) != len(pairs):
        die("pair_id must be unique")
    # anchoring (cheap) for every pair first
    rows = []
    for p in pairs:
        fai = read_fai(p["fasta"] + ".fai")
        A, B = p["record_a"], p["record_b"]
        for x in (A, B):
            if x not in fai:
                die(f"{p['pair_id']}: {x} not in {p['fasta']}")
        an = anchors_by_record(p["anchors"], {A, B})
        ga, gb = set(an[A]), set(an[B])
        shared = ga & gb
        if not shared:
            die(f"{p['pair_id']}: {A} and {B} share no anchored gene")
        sa, sb = span(shared, an[A]), span(shared, an[B])
        refc = sorted({an[A][g][1] for g in shared})
        refp = sorted(an[A][g][2] for g in shared)
        r = dict(pair_id=p["pair_id"], line=p["line"], unit=p["unit"], record_a=A, len_a=fai[A], record_b=B,
                 len_b=fai[B], genes_a=len(ga), genes_b=len(gb), genes_shared=len(shared),
                 shared_frac_of_a=round(len(shared) / len(ga), 4), genes_a_only=len(ga - gb),
                 shared_span_a=f"{sa[0]}-{sa[1]}", shared_span_b=f"{sb[0]}-{sb[1]}",
                 frac_a_in_shared_span=round((sa[1] - sa[0] + 1) / fai[A], 4),
                 frac_b_in_shared_span=round((sb[1] - sb[0] + 1) / fai[B], 4),
                 ref_span=f"{','.join(refc)}:{refp[0]}-{refp[-1]}")
        rb = region(sb, fai[B], a.flank_bp)
        r["b_region"] = f"{B}:{rb[0]}-{rb[1]}"
        r["_jobs"] = [("AB", (p["fasta"], A, 1, fai[A]), (p["fasta"], B, rb[0], rb[1]))]
        if p.get("other_records") and p.get("other_fasta") and p.get("other_anchors"):
            ofai = read_fai(p["other_fasta"] + ".fai")
            orecs = [x for x in p["other_records"].split(",") if x]
            oan = anchors_by_record(p["other_anchors"], set(orecs))
            best = max(orecs, key=lambda x: len(set(oan[x]) & shared))
            sc = span(set(oan[best]) & shared, oan[best])
            if sc == (0, 0):
                r["c_region"] = "no counterpart shares genes"
                r["_search"] = p["other_fasta"]
            else:
                rc = region(sc, ofai[best], a.flank_bp)
                r["c_region"] = f"{p.get('other_unit', '')} {best}:{rc[0]}-{rc[1]}"
                r["genes_shared_with_c"] = len(set(oan[best]) & ga)
                r["_jobs"] += [("AC", (p["fasta"], A, 1, fai[A]), (p["other_fasta"], best, rc[0], rc[1])),
                               ("BC", (p["fasta"], B, rb[0], rb[1]), (p["other_fasta"], best, rc[0], rc[1]))]
        elif p.get("other_fasta") and a.search_other:
            r["c_region"] = "no anchored counterpart in the other haplotype"
            r["_search"] = p["other_fasta"]
        if p.get("windows"):
            per = read_windows(p["windows"])
            da = depth(per, A, 1, fai[A])
            db = depth(per, B, sb[0], sb[1])
            fl = [w for w in per.get(B, []) if (sb[0] - 5_000_000 <= w[0] < sb[0]) or (sb[1] < w[1] <= sb[1] + 5_000_000)]
            fl_d = sorted(w[2] for w in fl)
            r.update(depth_a_median=da["median"], depth_a_classes=da["classes"], depth_b_median=db["median"],
                     depth_b_classes=db["classes"], depth_b_flank_median=round(fl_d[len(fl_d) // 2], 3) if fl_d else "")
        if p.get("other_windows") and r.get("c_region", "").count(":"):
            per = read_windows(p["other_windows"])
            cr = r["c_region"].split(" ")[-1]
            c, iv = cr.split(":")
            s, e = map(int, iv.split("-"))
            r["depth_c_median"] = depth(per, c, s, e)["median"]
        rows.append(r)
    if a.command == "plan":
        for r in rows:
            for name, q, t in r["_jobs"]:
                print(f"{r['pair_id']}\t{name}\tquery {q[1]}:{q[2]}-{q[3]} ({(q[3] - q[2] + 1) / 1e6:.1f} Mb)\t"
                      f"target {t[1]}:{t[2]}-{t[3]} ({(t[3] - t[2] + 1) / 1e6:.1f} Mb)")
        srch = defaultdict(list)
        for r in rows:
            if r.get("_search"):
                srch[r["_search"]].append(r["pair_id"])
        for fa, ids in srch.items():
            print(f"search\t{','.join(ids)}\tagainst whole {fa}")
        print(f"plan: {sum(len(r['_jobs']) for r in rows)} pairwise alignments, {len(srch)} whole-haplotype "
              f"searches; nothing run", file=sys.stderr)
        return 0
    if not a.out:
        die("run needs --out")
    if os.path.exists(a.out) or os.path.exists(a.out + ".partial"):
        die(f"{a.out} (or .partial) exists")
    for tool in (a.minimap2, a.samtools):
        if not shutil.which(tool):
            die(f"{tool} not found")
    start = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    part = a.out + ".partial"
    work = os.path.join(part, "work")
    os.makedirs(os.path.join(part, "paf"))
    os.makedirs(work)
    log = open(os.path.join(part, "commands.log"), "w")

    def extract(fa, rec, s, e, path):
        with open(path, "w") as fh:
            cmd = [a.samtools, "faidx", fa, f"{rec}:{s}-{e}"]
            log.write(" ".join(cmd) + f" > {path}\n")
            subprocess.run(cmd, stdout=fh, check=True)

    for r in rows:
        for name, q, t in r.pop("_jobs"):
            qf = os.path.join(work, f"{r['pair_id']}.{name}.q.fa")
            tf = os.path.join(work, f"{r['pair_id']}.{name}.t.fa")
            extract(*q, qf)
            extract(*t, tf)
            paf = os.path.join(part, "paf", f"{r['pair_id']}.{name}.paf")
            cmd = [a.minimap2, "-x", "asm20", "-c", "--secondary=no", "-t", str(a.threads), tf, qf]
            log.write(" ".join(cmd) + f" > {paf}\n")
            log.flush()
            with open(paf, "w") as fh:
                subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)
            st = paf_stats(paf, q[3] - q[2] + 1, t[3] - t[2] + 1, a.min_aln_bp)
            for k, v in st.items():
                r[f"{name}_{k}"] = v
            for f in (qf, tf):  # extracted sequence is a disposable intermediate of this run
                os.remove(f)
    # whole-other-haplotype search for candidates without an anchored counterpart (batched per FASTA)
    srch = defaultdict(list)
    for r in rows:
        if r.get("_search"):
            srch[r.pop("_search")].append(r)
    for k, (fa, rs) in enumerate(srch.items()):
        qf = os.path.join(work, f"search{k}.q.fa")
        with open(qf, "w") as fh:
            for r in rs:
                cmd = [a.samtools, "faidx", [p for p in pairs if p["pair_id"] == r["pair_id"]][0]["fasta"],
                       r["record_a"]]
                log.write(" ".join(cmd) + f" >> {qf}\n")
                subprocess.run(cmd, stdout=fh, check=True)
        paf = os.path.join(part, "paf", f"search{k}.{os.path.basename(fa)}.paf")
        cmd = [a.minimap2, "-x", "asm20", "-c", "--secondary=no", "-I", "8G", "-t", str(a.threads), fa, qf]
        log.write(" ".join(cmd) + f" > {paf}\n")
        log.flush()
        with open(paf, "w") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.remove(qf)
        for r in rs:
            sub = os.path.join(work, f"{r['pair_id']}.AO.paf")
            best = defaultdict(int)
            with open(paf) as fh, open(sub, "w") as out:
                for ln in fh:
                    f = ln.split("\t")
                    if f[0] == r["record_a"]:
                        out.write(ln)
                        if int(f[10]) >= a.min_aln_bp and "tp:A:P" in ln:
                            best[f[5]] += int(f[10])
            st = paf_stats(sub, r["len_a"], 1, a.min_aln_bp)
            os.remove(sub)
            r.update(AO_identity=st["identity"], AO_gc_identity=st["gc_identity"], AO_q_cov=st["q_cov"],
                     AO_aln_bp=st["aln_bp"],
                     AO_top_targets=";".join(f"{t}:{v / 1e6:.2f}Mb" for t, v in
                                             sorted(best.items(), key=lambda kv: -kv[1])[:3]))
    for r in rows:
        r["classification"], r["confidence"], r["basis"] = classify(r, P)
    os.rmdir(work)
    log.close()
    cols = COLS
    write_tsv(os.path.join(part, "dup_pairs.tsv"), rows, cols)
    ver = subprocess.run([a.minimap2, "--version"], capture_output=True, text=True).stdout.strip()
    with open(os.path.join(part, "RUN.json"), "w") as fh:
        json.dump(dict(schema="g2f-dup-pairs/1", code_version=CODE_VERSION, start=start,
                       end=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
                       command=[sys.executable] + sys.argv, cwd=os.getcwd(), params=P, minimap2=ver,
                       pairs=os.path.abspath(a.pairs), pairs_sha256=sha256(a.pairs), script_sha256=sha256(__file__)),
                  fh, indent=1)
    os.rename(part, a.out)
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
