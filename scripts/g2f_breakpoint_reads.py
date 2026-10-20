#!/usr/bin/env python3
"""Read evidence at candidate breakpoints in an assembly.

For each point (contig, position) and each alignment file (BAM or CRAM), count
primary alignments with MAPQ >= --min-mapq whose aligned span covers
[position - flank, position + flank] for every --flank given, reads with a
clip of >= --min-clip within --clip-window of the position, and the median
depth of MAPQ-filtered reads at the point.

  g2f_breakpoint_reads.py --points POINTS.tsv --reads NAME=FILE.bam [--reads NAME=FILE.cram ...]
        [--reference ASM.fa] [--flank 5000 --flank 20000] --out OUT.tsv

POINTS.tsv: label<TAB>contig<TAB>position (1-based), '#' comments allowed.
Spanning reads show that reads are continuous across the point in this
assembly; inside repeats longer than the reads they cannot exclude a misjoin,
and reads that were assembly input test consistency, not independent truth.
"""
import argparse
import os
import statistics
import sys

import pysam


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", required=True)
    ap.add_argument("--reads", action="append", required=True, help="NAME=alignment file")
    ap.add_argument("--reference", help="assembly FASTA (needed for CRAM)")
    ap.add_argument("--flank", type=int, action="append")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-clip", type=int, default=500)
    ap.add_argument("--clip-window", type=int, default=2000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    flanks = a.flank or [5000, 20000]
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; results are never overwritten")
    points = []
    with open(a.points) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            lab, c, p = line.rstrip("\n").split("\t")[:3]
            points.append((lab, c, int(p)))
    cols = ["label", "contig", "position", "reads", "depth_mapq", "clipped_near"] + [f"spanning_{f}" for f in flanks]
    rows = []
    for spec in a.reads:
        name, _, path = spec.partition("=")
        aln = pysam.AlignmentFile(path, reference_filename=a.reference) if path.endswith(".cram") \
            else pysam.AlignmentFile(path)
        for lab, c, p in points:
            p0 = p - 1
            span = {f: 0 for f in flanks}
            clipped = 0
            reach = max(flanks) + a.clip_window
            for r in aln.fetch(c, max(0, p0 - reach), p0 + reach):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
                    continue
                for f in flanks:
                    if r.reference_start <= p0 - f and r.reference_end >= p0 + f:
                        span[f] += 1
                cig = r.cigartuples or []
                if cig and cig[0][0] in (4, 5) and cig[0][1] >= a.min_clip and abs(r.reference_start - p0) <= a.clip_window:
                    clipped += 1
                elif cig and cig[-1][0] in (4, 5) and cig[-1][1] >= a.min_clip and abs(r.reference_end - p0) <= a.clip_window:
                    clipped += 1
            depth = [col.get_num_aligned() for col in aln.pileup(c, p0, p0 + 1, truncate=True, min_mapping_quality=a.min_mapq,
                                                                 stepper="nofilter", ignore_orphans=False)]
            rows.append([lab, c, p, name, depth[0] if depth else 0, clipped] + [span[f] for f in flanks])
        aln.close()
    tmp = a.out + ".partial"
    with open(tmp, "w") as fh:
        fh.write(f"# min_mapq={a.min_mapq} min_clip={a.min_clip} clip_window={a.clip_window}; spanning_F = primary reads "
                 "aligned over [pos-F, pos+F]\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    os.replace(tmp, a.out)
    for r in rows:
        print("\t".join(map(str, r)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
