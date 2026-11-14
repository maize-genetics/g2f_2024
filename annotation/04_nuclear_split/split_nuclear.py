#!/usr/bin/env python3
"""
Step 04b: write the nuclear genome of each assembly and its gene models, and label
every gene by the scaffold it sits on. The whole-genome files are not changed.

Chromosomes are found by name (chr1..chr10, also chr01, <unit>_chr01, and split chromosomes such
as chr4_part1). Each assembly must have all 10, and each must be labelled nuclear by the organelle
screen. Unplaced scaffolds are not part of the nuclear genome.

Gene groups (gene_classes.tsv):
  nuclear            on the chromosomes
  unplaced_nuclear   on an unplaced scaffold labelled nuclear
  plastid / mitochondrion / review   from the scaffold label

Reads   PROJECT_DIR/units.tsv
        PROJECT_DIR/inputs/<unit>.scaffolds.fa
        PROJECT_DIR/scaffold_classes/<unit>.scaffold_classification.tsv           (organelle_screen.sbatch)
        PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.psauron.gff3              (step 03)
        PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.cds.fa, .proteins.fa      (step 02)
Writes  PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.fa                       nuclear genome sequence
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.helixer.gff3             nuclear genes, with PSAURON scores
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.helixer.cds.fa
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.helixer.proteins.fa
        PROJECT_DIR/raw/<unit>/<unit>.gene_classes.tsv                            one row per gene, whole genome
        PROJECT_DIR/nuclear_genome/gene_classes_summary.tsv                       genes per group, per unit

Usage (from the repository folder):
  source config.sh && python3 04_nuclear_split/split_nuclear.py --project "$PROJECT_DIR"
"""
import argparse
import os
import re
import sys
from collections import Counter

CHROM = re.compile(r"(^|_)chr0?([1-9]|10)(_part\d+)?$")
GROUPS = ["nuclear", "unplaced_nuclear", "plastid", "mitochondrion", "review"]


def attrs(c9):
    return dict(p.split("=", 1) for p in c9.strip().rstrip(";").split(";") if "=" in p)


def fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def write_fa(path, recs, width=60):
    n = 0
    with open(path + ".part", "w") as o:
        for name, s in recs:
            o.write(f">{name}\n")
            for i in range(0, len(s), width):
                o.write(s[i:i + width] + "\n")
            n += 1
    os.replace(path + ".part", path)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=os.environ.get("PROJECT_DIR"), help="PROJECT_DIR from config.sh")
    ap.add_argument("units", nargs="*", help="units to process (default: all in units.tsv)")
    a = ap.parse_args()
    if not a.project:
        sys.exit("give --project (PROJECT_DIR in config.sh)")
    B = a.project
    units = a.units or [l.split("\t")[0].strip() for l in open(f"{B}/units.tsv") if l.strip()]
    NG = f"{B}/nuclear_genome"
    os.makedirs(NG, exist_ok=True)

    summary = []
    for u in units:
        raw = f"{B}/raw/{u}"
        gff = f"{raw}/{u}.scaffolds.helixer.psauron.gff3"
        if not os.path.exists(gff):
            sys.exit(f"{u}: missing {gff}; run step 03 (psauron_add_to_gff3.py) first")

        # scaffold labels
        lab, length = {}, {}
        with open(f"{B}/scaffold_classes/{u}.scaffold_classification.tsv") as fh:
            h = fh.readline().rstrip("\n").split("\t")
            i_s, i_l, i_k = h.index("scaffold"), h.index("length"), h.index("final_class")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                lab[f[i_s]], length[f[i_s]] = f[i_k], int(f[i_l])
        chroms = {s for s in lab if CHROM.search(s)}
        nums = {int(CHROM.search(s).group(2)) for s in chroms}
        if nums != set(range(1, 11)):
            sys.exit(f"{u}: chromosomes found for {sorted(nums)}, expected 1-10")
        bad = sorted(s for s in chroms if lab[s] != "nuclear")
        if bad:
            sys.exit(f"{u}: chromosome scaffolds not labelled nuclear: {bad}")

        def group(s):
            if s in chroms:
                return "nuclear"
            k = lab[s]
            return "unplaced_nuclear" if k == "nuclear" else k if k in GROUPS else "review"

        out = f"{NG}/{u}"
        os.makedirs(out, exist_ok=True)

        # nuclear genome FASTA, chromosomes in the order of the assembly
        n_chr = write_fa(f"{out}/{u}.nuclear.fa", ((n, s) for n, s in fasta(f"{B}/inputs/{u}.scaffolds.fa") if n in chroms))
        if n_chr != len(chroms):
            sys.exit(f"{u}: {n_chr} chromosome sequences written, {len(chroms)} expected")

        # GFF3: keep the header and every feature on the chromosomes
        genes, mrna_of, nuc_lines = {}, {}, []
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    if line.startswith("##gff-version") or (line.startswith("##sequence-region")
                                                            and line.split()[1] in chroms):
                        nuc_lines.append(line)
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[0] not in lab:
                    sys.exit(f"{u}: GFF3 sequence {f[0]} is not in the scaffold labels")
                if f[0] in chroms:
                    nuc_lines.append(line)
                x = attrs(f[8])
                if f[2] == "gene":
                    genes[x["ID"]] = (f[0], f[3], f[4], f[6])
                elif f[2] == "mRNA":
                    mrna_of[x["Parent"]] = x["ID"]
        with open(f"{out}/{u}.nuclear.helixer.gff3.part", "w") as o:
            o.writelines(nuc_lines)
        os.replace(f"{out}/{u}.nuclear.helixer.gff3.part", f"{out}/{u}.nuclear.helixer.gff3")

        # gene table, whole genome
        counts, keep = Counter(), set()
        with open(f"{raw}/{u}.gene_classes.tsv", "w") as o:
            o.write("gene\tmrna\tscaffold\tstart\tend\tstrand\tlocation\tclass\tgroup\n")
            for g, (s, b, e, st) in genes.items():
                grp = group(s)
                counts[grp] += 1
                m = mrna_of.get(g, "")
                if grp == "nuclear":
                    keep.add(m)
                o.write("\t".join([g, m, s, b, e, st, "chromosome" if s in chroms else "unplaced", lab[s], grp]) + "\n")

        # CDS and proteins of the nuclear genes
        for kind in ("cds", "proteins"):
            n = write_fa(f"{out}/{u}.nuclear.helixer.{kind}.fa",
                         ((n, s) for n, s in fasta(f"{raw}/{u}.scaffolds.helixer.{kind}.fa") if n in keep))
            if n != len(keep):
                sys.exit(f"{u}: {kind}: {n} sequences for {len(keep)} nuclear mRNAs")

        summary.append([u, len(genes)] + [counts[g] for g in GROUPS])
        print(f"{u}: {n_chr} chromosomes, {sum(length[s] for s in chroms) / 1e6:,.1f} Mb; {len(genes):,} genes, "
              f"nuclear {counts['nuclear']:,}, unplaced nuclear {counts['unplaced_nuclear']:,}, "
              f"plastid {counts['plastid']:,}, mito {counts['mitochondrion']:,}, review {counts['review']:,}")

    with open(f"{NG}/gene_classes_summary.tsv", "w") as o:
        o.write("unit\tall_genes\t" + "\t".join(GROUPS) + "\n")
        for r in summary:
            o.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {NG}/gene_classes_summary.tsv")


if __name__ == "__main__":
    main()
