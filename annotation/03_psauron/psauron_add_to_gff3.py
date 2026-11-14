#!/usr/bin/env python3
"""
Step 03b: add each transcript's PSAURON result to the whole-genome Helixer GFF3 (column 9, the
attributes), and summarise the scores per genome.

  mRNA line:  ...;psauron_score=0.9731;psauron_is_protein=True
  gene line:  the same values (Helixer predicts one mRNA per gene)
  transcripts PSAURON did not score (shorter than its minimum length) get psauron_score=NA

Reads   PROJECT_DIR/units.tsv
        PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.raw.gff3        (step 01)
        PROJECT_DIR/psauron/<unit>.helixer.cds.psauron.csv              (psauron_run.sbatch)
Writes  PROJECT_DIR/raw/<unit>/<unit>.scaffolds.helixer.psauron.gff3
        PROJECT_DIR/psauron/psauron_summary.tsv   per genome: genes, scored, % protein-coding, mean / median score

Usage (from the repository folder):
  source config.sh && python3 03_psauron/psauron_add_to_gff3.py --project "$PROJECT_DIR"
"""
import argparse
import csv
import os
import re
import statistics
import sys


def attrs(c9):
    return dict(x.split("=", 1) for x in c9.strip().rstrip(";").split(";") if "=" in x)


def read_csv(path):
    """Return {transcript: (in_frame_score, is_protein)}.
    PSAURON writes a few lines (command, overall score, a note) before the header row."""
    lines = open(path).read().splitlines()
    try:
        h = next(i for i, l in enumerate(lines) if l.lower().startswith("description,"))
    except StopIteration:
        sys.exit(f"no header row starting with 'description,' in {path}")
    rd = csv.reader(lines[h:])
    hdr = [x.strip().lower() for x in next(rd)]
    for col in ("in_frame_score", "psauron_is_protein"):
        if col not in hdr:
            sys.exit(f"column {col} not found in {path}; header = {hdr}")
    si, pi = hdr.index("in_frame_score"), hdr.index("psauron_is_protein")
    res = {}
    for r in rd:
        if r:
            res[r[0].split()[0].lstrip(">")] = (float(r[si]), r[pi].strip())
    return res


def add(c9, score, isp):
    c9 = re.sub(r";?psauron_(score|is_protein)=[^;]*", "", c9.rstrip("\n").rstrip(";"))
    return f"{c9};psauron_score={score};psauron_is_protein={isp}\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default=os.environ.get("PROJECT_DIR"), help="PROJECT_DIR from config.sh")
    ap.add_argument("units", nargs="*", help="units to process (default: all in units.tsv)")
    a = ap.parse_args()
    if not a.project:
        sys.exit("give --project (PROJECT_DIR in config.sh)")
    B = a.project

    units = a.units or [l.split("\t")[0].strip() for l in open(f"{B}/units.tsv") if l.strip()]
    summ_path = f"{B}/psauron/psauron_summary.tsv"
    with open(summ_path, "w") as summ:
        summ.write("unit\tgenes\tscored\tnot_scored\tpct_protein_coding\tmean_score\tmedian_score\n")
        for u in units:
            gff_in = f"{B}/raw/{u}/{u}.scaffolds.helixer.raw.gff3"
            gff_out = f"{B}/raw/{u}/{u}.scaffolds.helixer.psauron.gff3"
            res = read_csv(f"{B}/psauron/{u}.helixer.cds.psauron.csv")
            gff = open(gff_in).read().splitlines(keepends=True)

            # mRNA -> gene, so each gene line gets its transcript's values
            tx_gene = {}
            for l in gff:
                f = l.split("\t")
                if len(f) > 8 and f[2] == "mRNA":
                    x = attrs(f[8])
                    tx_gene[x["ID"]] = x["Parent"]
            extra = set(res) - set(tx_gene)
            if extra:
                sys.exit(f"{u}: {len(extra)} PSAURON IDs are not mRNA IDs in the GFF3, e.g. {sorted(extra)[:3]}")
            gene_val = {g: res[t] for t, g in tx_gene.items() if t in res}

            vals = []
            with open(gff_out + ".part", "w") as o:
                for l in gff:
                    f = l.split("\t")
                    if not l.startswith("#") and len(f) > 8 and f[2] in ("gene", "mRNA"):
                        key = attrs(f[8])["ID"]
                        v = res.get(key) if f[2] == "mRNA" else gene_val.get(key)
                        f[8] = add(f[8], *((f"{v[0]:.4f}", v[1]) if v else ("NA", "NA")))
                        l = "\t".join(f)
                        if f[2] == "gene":
                            vals.append(v)
                    o.write(l)
            os.replace(gff_out + ".part", gff_out)

            sc = [v[0] for v in vals if v]
            pc = sum(1 for v in vals if v and v[1].lower() == "true")
            nan = float("nan")
            summ.write(f"{u}\t{len(vals)}\t{len(sc)}\t{len(vals) - len(sc)}\t"
                       f"{100 * pc / len(sc) if sc else nan:.1f}\t"
                       f"{statistics.mean(sc) if sc else nan:.3f}\t{statistics.median(sc) if sc else nan:.3f}\n")
            print(f"{u}: {len(vals):,} genes, {len(sc):,} scored, {len(vals) - len(sc):,} not scored (NA), "
                  f"{100 * pc / max(len(sc), 1):.1f}% protein-coding")
    print(f"wrote {summ_path}")


if __name__ == "__main__":
    main()
