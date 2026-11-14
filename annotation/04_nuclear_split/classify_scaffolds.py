#!/usr/bin/env python3
"""
Step 04a (called by organelle_screen.sbatch): label every scaffold of one assembly as nuclear,
plastid or mitochondrion from its minimap2 alignments to the organelle genomes.

  covered fraction = bases of the scaffold in alignments with identity >= 0.95 and length >= 500 bp,
                     divided by the scaffold length (overlapping alignments counted once)
  plastid fraction >= 0.80 -> plastid ; mitochondrial fraction >= 0.80 -> mitochondrion
  plastid + mitochondrial >= 0.30 -> review ; otherwise nuclear

'review' scaffolds need a decision in review_overrides.tsv (unit, scaffold, final_class, confidence,
reason). Nuclear scaffolds with >= 100 kb matching the organelle genomes are reported in the log.

Usage: classify_scaffolds.py UNIT scaffolds.fa alignments.paf(.gz) out.tsv [overrides.tsv]
       [--plastid NC_001666] [--mito NC_007982]   (sequence IDs of the two genomes in ORGANELLE_FASTA)
Output columns: scaffold length plastid_frac mito_frac organelle_frac screen_class final_class decided_by
Prints one summary line: unit, scaffolds and Mb per final class, and any unresolved scaffolds.
"""
import argparse
import gzip
import sys
from collections import defaultdict

MIN_ID, MIN_ALN, CALL, REVIEW = 0.95, 500, 0.80, 0.30
COLS = ["scaffold", "length", "plastid_frac", "mito_frac", "organelle_frac",
        "screen_class", "final_class", "decided_by"]


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def union_len(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if ce is not None else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("unit")
    ap.add_argument("fasta")
    ap.add_argument("paf")
    ap.add_argument("out")
    ap.add_argument("overrides", nargs="?")
    ap.add_argument("--plastid", default="NC_001666", help="plastid sequence ID (version suffix ignored)")
    ap.add_argument("--mito", default="NC_007982", help="mitochondrion sequence ID (version suffix ignored)")
    a = ap.parse_args()

    lengths, name = {}, None
    with opener(a.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())

    hits = {"pt": defaultdict(list), "mt": defaultdict(list)}
    with opener(a.paf) as fh:
        for line in fh:
            f = line.split("\t")
            matches, alen = int(f[9]), int(f[10])
            if alen < MIN_ALN or matches / alen < MIN_ID:
                continue
            t = f[5].split(".")[0]
            key = "pt" if t == a.plastid else "mt" if t == a.mito else None
            if key:
                hits[key][f[0]].append((int(f[2]), int(f[3])))

    ov = {}
    if a.overrides:
        with open(a.overrides) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4 and f[0] == a.unit:
                    ov[f[1]] = (f[2], f[3])
    unknown = set(ov) - set(lengths)
    if unknown:
        sys.exit(f"{a.unit}: overrides name scaffolds not in the FASTA: {sorted(unknown)}")

    rows = []
    for s, L in lengths.items():
        p = union_len(hits["pt"].get(s, [])) / L
        m = union_len(hits["mt"].get(s, [])) / L
        b = union_len(hits["pt"].get(s, []) + hits["mt"].get(s, [])) / L
        k = "plastid" if p >= CALL else "mitochondrion" if m >= CALL else "review" if b >= REVIEW else "nuclear"
        final, by = (ov[s][0], f"override_{ov[s][1]}") if s in ov else (k, "screen")
        rows.append({"scaffold": s, "length": L, "plastid_frac": f"{p:.3f}", "mito_frac": f"{m:.3f}",
                     "organelle_frac": f"{b:.3f}", "screen_class": k, "final_class": final, "decided_by": by})
        kb = b * L / 1e3
        if final == "nuclear" and kb >= 100:
            print(f"{a.unit}: CHECK nuclear scaffold {s} has {kb:.0f} kb matching organelle genomes", file=sys.stderr)

    with open(a.out, "w") as o:
        o.write("\t".join(COLS) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in COLS) + "\n")

    agg = defaultdict(lambda: [0, 0])
    for r in rows:
        agg[r["final_class"]][0] += 1
        agg[r["final_class"]][1] += r["length"]
    line = [a.unit]
    for c in ("nuclear", "plastid", "mitochondrion"):
        line += [agg[c][0], f"{agg[c][1] / 1e6:.2f}"]
    open_ = [r["scaffold"] for r in rows if r["final_class"] not in ("nuclear", "plastid", "mitochondrion")]
    line += [len(open_), ",".join(open_) or "-"]
    print("\t".join(map(str, line)))


if __name__ == "__main__":
    main()
