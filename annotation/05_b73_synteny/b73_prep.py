#!/usr/bin/env python3
"""
Step 05a: MCScan (jcvi) inputs for B73 against each assembly, using the nuclear genes.

MCScan needs, per genome, two files with the same gene names:
  <name>.bed   one line per gene: sequence, start (0-based), end, name, score, strand
  <name>.cds   coding sequences, FASTA headers = those same names

B73:   canonical transcripts of Zm00001eb.1 (one per gene); genes on organelle sequences
       (Pt/Mt) are left out, to match the nuclear gene sets.
Units: Helixer mRNAs (one per gene) of the nuclear genome, from step 04; Helixer's own IDs.

Reads   B73_GFF3, B73_CDS, B73_CANONICAL (config.sh)
        PROJECT_DIR/nuclear_genome/<unit>/<unit>.nuclear.helixer.gff3, .cds.fa
Writes  PROJECT_DIR/b73_synteny/shared/B73.bed, B73.cds
        PROJECT_DIR/b73_synteny/runs/<unit>/<unit>.bed, <unit>.cds, and links to the B73 files
        PROJECT_DIR/b73_synteny/prep_summary.tsv

Usage (from the repository folder):
  source config.sh && python3 05_b73_synteny/b73_prep.py --project "$PROJECT_DIR" \
      --b73-gff3 "$B73_GFF3" --b73-cds "$B73_CDS" --b73-canonical "$B73_CANONICAL"
"""
import argparse
import gzip
import os
import re
import sys

ORGANELLE_SEQ = re.compile(r"(^|_)(chr)?(mt|pt|mito|chloro|plastid|mitochondrion|chloroplast)", re.I)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def attrs(c9):
    return dict(p.split("=", 1) for p in c9.strip().rstrip(";").split(";") if "=" in p)


def read_fasta(path):
    name, seq = None, []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name:
        yield name, "".join(seq)


def mrna_bed(gff, keep=None, drop_seq=None):
    rows, dropped = [], 0
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            tid = attrs(f[8])["ID"]
            if keep is not None and tid not in keep:
                continue
            if drop_seq and drop_seq.search(f[0]):
                dropped += 1
                continue
            rows.append([f[0], int(f[3]) - 1, int(f[4]), tid, 0, f[6]])
    return rows, dropped


def write_pair(prefix, bed_rows, cds):
    """Write .bed and .cds with only the names present in both; return the count."""
    rows = sorted((r for r in bed_rows if r[3] in cds), key=lambda r: (r[0], r[1]))
    with open(prefix + ".bed", "w") as b:
        for r in rows:
            b.write("\t".join(map(str, r)) + "\n")
    with open(prefix + ".cds", "w") as f:
        for r in rows:
            s = cds[r[3]]
            f.write(f">{r[3]}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="PROJECT_DIR from config.sh")
    ap.add_argument("--b73-gff3", required=True)
    ap.add_argument("--b73-cds", required=True)
    ap.add_argument("--b73-canonical", required=True)
    a = ap.parse_args()
    B = a.project
    BO = f"{B}/b73_synteny"
    shared = f"{BO}/shared"
    os.makedirs(shared, exist_ok=True)

    # B73
    with opener(a.b73_canonical) as fh:
        canon = {l.split()[0] for l in fh if l.strip()}
    b73_cds = {n: s for n, s in read_fasta(a.b73_cds) if n in canon}
    b73_bed, b73_org = mrna_bed(a.b73_gff3, keep=canon, drop_seq=ORGANELLE_SEQ)
    n_b73 = write_pair(f"{shared}/B73", b73_bed, b73_cds)
    print(f"B73: {len(canon):,} canonical transcripts, {len(b73_cds):,} with CDS, "
          f"{b73_org} on organelle sequences left out, {n_b73:,} genes written")
    if n_b73 < 35000:
        sys.exit("Too few B73 genes; check the B73 files (do the IDs in the three files match?)")

    # units
    rows = [["unit", "nuclear_genes_in_gff3", "genes_written", "B73_genes_written", "check"]]
    for line in open(f"{B}/units.tsv"):
        if not line.strip():
            continue
        u = line.split("\t")[0].strip()
        run = f"{BO}/runs/{u}"
        os.makedirs(run, exist_ok=True)
        ng = f"{B}/nuclear_genome/{u}"
        gff = f"{ng}/{u}.nuclear.helixer.gff3"
        cds = dict(read_fasta(f"{ng}/{u}.nuclear.helixer.cds.fa"))
        bed, _ = mrna_bed(gff)
        n_genes = sum(1 for l in open(gff) if not l.startswith("#") and len(l.split("\t")) > 2 and l.split("\t")[2] == "gene")
        n = write_pair(f"{run}/{u}", bed, cds)
        for ext in ("bed", "cds"):
            link = f"{run}/B73.{ext}"
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(f"{shared}/B73.{ext}", link)
        ok = "OK" if n == n_genes == len(bed) else "FAIL: gene/mRNA/CDS counts differ"
        rows.append([u, n_genes, n, n_b73, ok])
        print(f"{u}: {n:,} genes written ({ok})")

    with open(f"{BO}/prep_summary.tsv", "w") as o:
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    if any(r[-1] != "OK" for r in rows[1:]):
        sys.exit("Some units failed; see prep_summary.tsv")


if __name__ == "__main__":
    main()
