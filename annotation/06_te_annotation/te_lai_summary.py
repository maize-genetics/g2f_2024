#!/usr/bin/env python3
"""
Step 06f: collect the nuclear-genome LAI of each unit into one table (after te_lai.sbatch).

Reads   PROJECT_DIR/te_annotation/runs/<unit>/lai/*.LAI   (the whole_genome row)
Writes  PROJECT_DIR/te_annotation/te_LAI.tsv              unit, intact LTR %, total LTR %, raw LAI, LAI

Usage (from the repository folder):
  source config.sh && python3 06_te_annotation/te_lai_summary.py --project "$PROJECT_DIR"
"""
import argparse
import glob


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="PROJECT_DIR from config.sh")
    B = ap.parse_args().project
    TE = f"{B}/te_annotation"
    units = [l.split("\t")[0].strip() for l in open(f"{B}/units.tsv") if l.strip()]
    with open(f"{TE}/te_LAI.tsv", "w") as o:
        o.write("unit\tintact_LTR_pct\ttotal_LTR_pct\traw_LAI\tLAI\n")
        for u in units:
            row = None
            files = sorted(glob.glob(f"{TE}/runs/{u}/lai/*.LAI"))
            if files:
                for line in open(files[-1]):
                    p = line.split()
                    if p and p[0] == "whole_genome":
                        row = p          # Chr From To Intact Total raw_LAI LAI
            if row:
                o.write(f"{u}\t{row[3]}\t{row[4]}\t{row[5]}\t{row[6]}\n")
                print(f"{u:16s} raw LAI {row[5]:>6}   LAI {row[6]:>6}")
            else:
                o.write(f"{u}\tNA\tNA\tNA\tNA\n")
                print(f"{u:16s} no LAI result yet")
    print(f"wrote {TE}/te_LAI.tsv")


if __name__ == "__main__":
    main()
