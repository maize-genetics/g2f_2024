#!/usr/bin/env python3
"""
Step 05d: combine the one-to-one B73 ortholog tables of all assemblies.

1. Combined table: one row per B73 canonical gene (in B73 order) with the matching gene in each
   assembly (empty where there is none), and in how many assemblies and lines it is found
   (haplotypes <line>_hap1 / _hap2 count as one line).
     b73_synteny/B73_orthologs_combined.tsv   (also copied to nuclear_genome/)
     b73_synteny/B73_orthologs_sharing.tsv    number of B73 genes found in 0..N assemblies / lines

2. Chromosome check: is each assembly gene on the same chromosome as its B73 partner?
     b73_synteny/chromosome_check.tsv              one row per assembly
     b73_synteny/chromosome_check_discordant.tsv   every pair on different chromosomes
   If the assemblies were scaffolded against B73, this checks agreement between scaffolding and
   synteny, not an independent truth.

Usage (from the repository folder, after all b73_mcscan.sbatch tasks have finished):
  source config.sh && python3 05_b73_synteny/b73_combined.py --project "$PROJECT_DIR"
"""
import argparse
import re
import shutil
from collections import Counter, defaultdict


def chrom_number(name):
    m = re.search(r"(?:^|_)chr0*(\d+)(?:_|$)", name)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="PROJECT_DIR from config.sh")
    B = ap.parse_args().project
    BO = f"{B}/b73_synteny"

    units = [l.split("\t")[0].strip() for l in open(f"{B}/units.tsv") if l.strip()]
    line_of = {u: re.sub(r"_hap\d+$", "", u) for u in units}
    lines = list(dict.fromkeys(line_of.values()))

    b73 = []
    for l in open(f"{BO}/shared/B73.bed"):
        f = l.rstrip("\n").split("\t")
        b73.append((re.sub(r"_T\d+$", "", f[3]), f[0], int(f[1]) + 1))
    b73.sort(key=lambda r: (chrom_number(r[1]) or 99, r[1], r[2]))

    match = defaultdict(dict)                     # B73 gene -> unit -> unit gene
    check_rows, discordant = [], []
    for u in units:
        n = placed = same = 0
        with open(f"{BO}/runs/{u}/{u}.B73_orthologs.full.tsv") as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                bg, bseq, ug, useq = f[0], f[2], f[6], f[8]
                match[bg][u] = ug
                n += 1
                uc, bc = chrom_number(useq), chrom_number(bseq)
                if uc is None or bc is None:
                    continue
                placed += 1
                if uc == bc:
                    same += 1
                else:
                    discordant.append([u, bg, bseq, f[3], ug, useq, f[9]])
        check_rows.append([u, n, placed, same, placed - same, f"{100 * same / placed:.2f}" if placed else "NA"])

    out = f"{BO}/B73_orthologs_combined.tsv"
    share_u, share_l = Counter(), Counter()
    with open(out, "w") as o:
        o.write("B73_gene\tB73_chr\tB73_start\tn_units\tn_lines\t" + "\t".join(units) + "\n")
        for g, c, s in b73:
            row = [match[g].get(u, "") for u in units]
            nu = sum(1 for x in row if x)
            nl = len({line_of[u] for u in units if match[g].get(u)})
            share_u[nu] += 1
            share_l[nl] += 1
            o.write("\t".join([g, c, str(s), str(nu), str(nl)] + row) + "\n")
    shutil.copyfile(out, f"{B}/nuclear_genome/B73_orthologs_combined.tsv")

    with open(f"{BO}/B73_orthologs_sharing.tsv", "w") as o:
        o.write("found_in\tB73_genes_by_units\tB73_genes_by_lines\n")
        for k in range(len(units), -1, -1):
            o.write(f"{k}\t{share_u.get(k, 0)}\t{share_l.get(k, 0) if k <= len(lines) else ''}\n")

    with open(f"{BO}/chromosome_check.tsv", "w") as o:
        o.write("unit\tpairs\tpairs_on_chromosomes\tsame_chromosome\tdifferent_chromosome\tpct_same\n")
        for r in check_rows:
            o.write("\t".join(map(str, r)) + "\n")
    with open(f"{BO}/chromosome_check_discordant.tsv", "w") as o:
        o.write("unit\tB73_gene\tB73_chr\tB73_start\tunit_gene\tunit_chr\tunit_start\n")
        for r in discordant:
            o.write("\t".join(map(str, r)) + "\n")

    print(f"B73 genes: {len(b73):,}   assemblies: {len(units)}   lines: {len(lines)}")
    print(f"found in all {len(units)} assemblies: {share_u[len(units)]:,}   "
          f"in all {len(lines)} lines: {share_l[len(lines)]:,}   in none: {share_u[0]:,}")
    print("\nunit\tpairs\tplaced\tsame_chr\tdiff_chr\tpct_same")
    for r in check_rows:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()
