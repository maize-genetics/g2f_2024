#!/usr/bin/env python3
"""
Step 06a: copy each nuclear genome (from step 04) into the TE folder with short sequence
names s00001, s00002, ..., and write the name map.

Why short names: EDTA 2.3.0 re-codes sequence names longer than 8 characters, while EDTA 2.2.2
(used for TIRs) does not, so with the original names the two versions would label the same
chromosome differently. te_restore_ids.py puts the real names back at the end.

Reads   PROJECT_DIR/units.tsv
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.fa
Writes  PROJECT_DIR/te_annotation/inputs/<unit>.nuclear.fa          short names
        PROJECT_DIR/te_annotation/inputs/<unit>.nuclear.idmap.tsv   short_id, original name, length
        PROJECT_DIR/te_annotation/inputs/prep_summary.tsv

Usage (from the repository folder):
  source config.sh && python3 06_te_annotation/te_prep_inputs.py --project "$PROJECT_DIR"
"""
import argparse
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="PROJECT_DIR from config.sh")
    ap.add_argument("units", nargs="*", help="units to prepare (default: all in units.tsv)")
    a = ap.parse_args()
    B = a.project
    out = f"{B}/te_annotation/inputs"
    os.makedirs(out, exist_ok=True)
    units = a.units or [l.split("\t")[0].strip() for l in open(f"{B}/units.tsv") if l.strip()]

    with open(f"{out}/prep_summary.tsv", "w") as summ:
        summ.write("unit\tsequences\tbp\n")
        for u in units:
            src = f"{B}/nuclear_genome/{u}/{u}.nuclear.fa"
            if not os.path.exists(src):
                sys.exit(f"{u}: missing {src}; run step 04 first")
            fa, mp = f"{out}/{u}.nuclear.fa", f"{out}/{u}.nuclear.idmap.tsv"
            n, lens, order, name = 0, {}, [], None
            with open(src) as fi, open(fa + ".part", "w") as fo:
                for line in fi:
                    if line.startswith(">"):
                        n += 1
                        name = line[1:].split()[0]
                        if name in lens:
                            sys.exit(f"{u}: duplicate sequence name {name}")
                        lens[name] = 0
                        order.append((f"s{n:05d}", name))
                        fo.write(f">s{n:05d}\n")
                    else:
                        lens[name] += len(line.strip())
                        fo.write(line)
            with open(mp + ".part", "w") as fm:
                fm.write("short_id\toriginal_name\tlength\n")
                for s, o in order:
                    fm.write(f"{s}\t{o}\t{lens[o]}\n")
            os.replace(fa + ".part", fa)
            os.replace(mp + ".part", mp)
            summ.write(f"{u}\t{n}\t{sum(lens.values())}\n")
            print(f"{u}: {n} sequences, {sum(lens.values()):,} bp")


if __name__ == "__main__":
    main()
