#!/usr/bin/env python3
"""
Step 05c (called by b73_mcscan.sbatch): turn one assembly's MCScan result into B73 ortholog tables.

Reads runs/<unit>/B73.<unit>.lifted.1x1.anchors (one-to-one syntenic orthologs; blocks separated
by '###'; anchors recovered by liftover carry an 'L' after the score), plus the two .bed files.

The quota 1:1 screen selects blocks, so inside a block a gene can still have several candidate
partners. These are resolved to one partner per gene (synteny anchors before lifted pairs, then
higher score); the pairs not kept are written to <unit>.B73_orthologs.dropped.tsv.

Gene IDs: B73 transcript Zm00001eb000010_T001 -> gene Zm00001eb000010;
          Helixer transcript <unit>_chr01_000001.1 -> gene <unit>_chr01_000001.

Writes  PROJECT_DIR/b73_synteny/runs/<unit>/<unit>.B73_orthologs.tsv          B73_gene, <unit>_gene
        PROJECT_DIR/b73_synteny/runs/<unit>/<unit>.B73_orthologs.full.tsv     + transcripts, locations, block, anchor type
        PROJECT_DIR/b73_synteny/runs/<unit>/<unit>.B73_orthologs.summary.tsv
        PROJECT_DIR/b73_synteny/runs/<unit>/<unit>.B73_orthologs.dropped.tsv
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.B73_orthologs.tsv            copy of the two-column table

Can be rerun on its own (no MCScan rerun needed):
  python3 05_b73_synteny/b73_table.py --project "$PROJECT_DIR" <unit>
"""
import argparse
import re
import shutil
from collections import Counter


def b73_gene(t):
    return re.sub(r"_T\d+$", "", t)


def unit_gene(t):
    return re.sub(r"\.\d+$", "", t)


def read_bed(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        d[f[3]] = (f[0], int(f[1]) + 1, int(f[2]), f[5])
    return d


def num(s):
    try:
        return float(s)
    except ValueError:
        return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="PROJECT_DIR from config.sh")
    ap.add_argument("unit")
    a = ap.parse_args()
    unit, B = a.unit, a.project
    run = f"{B}/b73_synteny/runs/{unit}"
    bb, ub = read_bed(f"{run}/B73.bed"), read_bed(f"{run}/{unit}.bed")

    pairs, block = [], 0
    for line in open(f"{run}/B73.{unit}.lifted.1x1.anchors"):
        if line.startswith("#"):
            block += 1
            continue
        f = line.split()
        if len(f) < 2:
            continue
        x, y = f[0], f[1]
        if x not in bb and y in bb:          # in case the columns are swapped
            x, y = y, x
        score = f[2] if len(f) > 2 else ""
        pairs.append((x, y, block, "lifted" if score.endswith("L") else "synteny", score.rstrip("L")))

    missing = [p for p in pairs if p[0] not in bb or p[1] not in ub]
    if missing:
        raise SystemExit(f"{unit}: {len(missing)} anchor names not found in the .bed files")

    candidates = len(pairs)
    ca, cb = Counter(p[0] for p in pairs), Counter(p[1] for p in pairs)
    multi_a = sum(1 for v in ca.values() if v > 1)
    multi_b = sum(1 for v in cb.values() if v > 1)

    # one partner per gene: synteny anchors first, then higher score
    ranked = sorted(pairs, key=lambda p: (p[3] != "synteny", -num(p[4]), p[0], p[1]))
    used_a, used_b, kept, dropped = set(), set(), [], []
    for p in ranked:
        if p[0] in used_a or p[1] in used_b:
            dropped.append((*p, "B73 gene already paired" if p[0] in used_a else f"{unit} gene already paired"))
            continue
        used_a.add(p[0]); used_b.add(p[1]); kept.append(p)
    pairs = sorted(kept, key=lambda p: p[0])
    assert len({p[0] for p in pairs}) == len(pairs) == len({p[1] for p in pairs})

    with open(f"{run}/{unit}.B73_orthologs.dropped.tsv", "w") as o:
        o.write(f"B73_transcript\t{unit}_transcript\tblock\tanchor_type\tscore\treason\n")
        for d in dropped:
            o.write("\t".join(map(str, d)) + "\n")

    two = f"{run}/{unit}.B73_orthologs.tsv"
    with open(two, "w") as o:
        o.write(f"B73_gene\t{unit}_gene\n")
        for x, y, *_ in pairs:
            o.write(f"{b73_gene(x)}\t{unit_gene(y)}\n")

    with open(f"{run}/{unit}.B73_orthologs.full.tsv", "w") as o:
        o.write(f"B73_gene\tB73_transcript\tB73_seq\tB73_start\tB73_end\tB73_strand\t"
                f"{unit}_gene\t{unit}_transcript\t{unit}_seq\t{unit}_start\t{unit}_end\t{unit}_strand\t"
                f"block\tanchor_type\tscore\n")
        for x, y, blk, kind, score in pairs:
            o.write("\t".join(map(str, [b73_gene(x), x, *bb[x], unit_gene(y), y, *ub[y], blk, kind, score])) + "\n")

    row = [unit, len(bb), len(ub), len(pairs),
           f"{100 * len(pairs) / len(bb):.1f}", f"{100 * len(pairs) / len(ub):.1f}",
           len({p[2] for p in pairs}), sum(1 for p in pairs if p[3] == "lifted"),
           candidates, len(dropped), multi_a, multi_b]
    with open(f"{run}/{unit}.B73_orthologs.summary.tsv", "w") as o:
        o.write("unit\tB73_genes\tunit_genes\tortholog_pairs\tpct_B73_with_ortholog\tpct_unit_with_ortholog\t"
                "blocks\tlifted_pairs\tcandidate_pairs\tdropped_pairs\t"
                "B73_genes_with_several_candidates\tunit_genes_with_several_candidates\n")
        o.write("\t".join(map(str, row)) + "\n")
    shutil.copyfile(two, f"{B}/nuclear_genome/{unit}/{unit}.B73_orthologs.tsv")
    print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
