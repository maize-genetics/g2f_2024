#!/usr/bin/env python3
"""
Step 06e (called by te_restore_ids.sbatch): put the real sequence names back into the EDTA and
RepeatMasker results (they were run on short names s00001, ...; see te_prep_inputs.py) and collect
the final files in results/.

Reads (in TE/runs/<unit>/):
  <unit>.nuclear.fa.mod.EDTA.TElib.fa      TE library (te_edta_combine)
  <unit>.nuclear.fa.mod.EDTA.intact.gff3   intact TEs (te_edta_combine)
  repeatmasker/<unit>.tar                  RepeatMasker .out, .out.gff, .tbl (te_repeatmasker)
  <unit>.softmasked.fa                     soft-masked genome (te_repeatmasker)
  edta_raw/<unit>.ltr.tar                  LTR_retriever pass lists (te_edta_ltr)
  TE/inputs/<unit>.nuclear.idmap.tsv
Writes TE/runs/<unit>/results/:
  <unit>.nuclear.EDTA.TElib.fa
  <unit>.nuclear.EDTA.intact.gff3
  <unit>.nuclear.RepeatMasker.out          all TE hits (LAI reads this)
  <unit>.nuclear.RepeatMasker.gff3         all TE hits
  <unit>.nuclear.RepeatMasker.tbl          % of the genome in each TE class
  <unit>.nuclear.softmasked.fa             TEs in lower case
  <unit>.LTR.pass.list                     all intact LTR retrotransposons (<genome>.mod.pass.list; input for LAI)
  <unit>.LTR.nmtf.pass.list                the subset without the canonical TGCA motif

Usage: te_restore_ids.py --te PROJECT_DIR/te_annotation <unit>
"""
import argparse
import os
import re
import sys
import tarfile

ID = re.compile(r"(?<![A-Za-z0-9])s\d{5}(?![0-9])")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--te", required=True, help="PROJECT_DIR/te_annotation")
    ap.add_argument("unit")
    a = ap.parse_args()
    u, TE = a.unit, a.te
    run = f"{TE}/runs/{u}"
    res = f"{run}/results"
    os.makedirs(res, exist_ok=True)
    g = f"{u}.nuclear.fa"

    idmap = {}
    with open(f"{TE}/inputs/{u}.nuclear.idmap.tsv") as fh:
        next(fh)
        for line in fh:
            s, o, _ = line.rstrip("\n").split("\t")
            idmap[s] = o

    def sub(text):
        return ID.sub(lambda m: idmap.get(m.group(0), m.group(0)), text)

    def convert(src, dst, whole=True):
        """copy src to dst with real names; whole=False changes FASTA headers only"""
        with open(src) as fi, open(dst + ".part", "w") as fo:
            for line in fi:
                fo.write(sub(line) if whole or line.startswith(">") else line)
        os.replace(dst + ".part", dst)
        print("wrote", dst)

    def from_tar(t, member, dst, rename=True):
        tmp = f"{res}/.{os.path.basename(member)}"
        with t.extractfile(member) as fi, open(tmp, "wb") as fo:
            fo.write(fi.read())
        if rename:
            convert(tmp, dst)
            os.remove(tmp)
        else:
            os.replace(tmp, dst)
            print("wrote", dst)

    convert(f"{run}/{g}.mod.EDTA.TElib.fa", f"{res}/{u}.nuclear.EDTA.TElib.fa", whole=False)
    convert(f"{run}/{g}.mod.EDTA.intact.gff3", f"{res}/{u}.nuclear.EDTA.intact.gff3")

    with tarfile.open(f"{run}/repeatmasker/{u}.tar") as t:
        from_tar(t, f"{g}.out", f"{res}/{u}.nuclear.RepeatMasker.out")
        from_tar(t, f"{g}.out.gff", f"{res}/{u}.nuclear.RepeatMasker.gff3")
        from_tar(t, f"{g}.tbl", f"{res}/{u}.nuclear.RepeatMasker.tbl", rename=False)

    convert(f"{run}/{u}.softmasked.fa", f"{res}/{u}.nuclear.softmasked.fa", whole=False)

    # pass lists from the LTR tar: the full list (<genome>.mod.pass.list) and, separately,
    # the non-TGCA subset (<genome>.mod.nmtf.pass.list); the names differ only by "nmtf"
    with tarfile.open(f"{run}/edta_raw/{u}.ltr.tar") as t:
        names = [m.name for m in t.getmembers()]
        full = [n for n in names if n.endswith(".mod.pass.list") and not n.endswith("nmtf.pass.list")]
        nmtf = [n for n in names if n.endswith(".mod.nmtf.pass.list")]
        if len(full) != 1:
            sys.exit(f"{u}: expected one .mod.pass.list in the LTR tar, found {full}")
        from_tar(t, full[0], f"{res}/{u}.LTR.pass.list")
        if nmtf:
            from_tar(t, nmtf[0], f"{res}/{u}.LTR.nmtf.pass.list")

    n = sum(1 for l in open(f"{res}/{u}.LTR.pass.list") if l.strip() and not l.startswith("#"))
    print(f"intact LTR retrotransposons in the pass list: {n:,}")

    # check: no short names left
    for f in sorted(os.listdir(res)):
        if f.endswith((".gff3", ".out", ".pass.list")):
            with open(f"{res}/{f}") as fh:
                left = sum(1 for line in fh if ID.search(line))
            print(f"{f}: lines still with a short name = {left}")


if __name__ == "__main__":
    main()
