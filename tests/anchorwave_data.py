"""Synthetic reference with annotated genes and a contig 'assembly' with planted
structural features for the AnchorWave workflow tests. All sequence is random;
nothing comes from real data.

Reference (genes every ~6 kb, 3 exons each, canonical GT..AG introns):
  chrA 300 kb, chrB 200 kb

Assembly (0.5-1% substitutions so CDS hits are not identical):
  ctg1  chrA[0:150k]                                   collinear
  ctg2  chrA[150k:300k] with [200k:260k] inverted      inversion inside one contig
  ctg3  chrB (whole)                                   collinear
  ctg4  chrB[50k:120k], second allele                  retained allele: needs -R 2
  ctg5  chrA[0:60k] + chrB[140k:200k]                  chimeric contig (chromosome switch)
"""
import gzip
import os
import random

STOP = ("TAA", "TAG", "TGA")
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOP]


def rseq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate, rng):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def make_chrom(name, length, rng, gene_prefix):
    """Random chromosome with genes on alternating strands; returns seq and GFF rows."""
    seq = list(rseq(length, rng))
    rows = [(name, "synthetic", "chromosome", 1, length, ".", ".", ".", f"ID={name}")]
    pos, k = 2000, 0
    while pos + 4000 < length:
        k += 1
        gid = f"{gene_prefix}{k:04d}"
        exons_cds = ["ATG" + "".join(rng.choice(CODONS) for _ in range(59))]
        exons_cds.append("".join(rng.choice(CODONS) for _ in range(60)))
        exons_cds.append("".join(rng.choice(CODONS) for _ in range(59)) + rng.choice(STOP))
        introns = ["GT" + rseq(300, rng) + "AG" for _ in range(2)]
        gene = exons_cds[0] + introns[0] + exons_cds[1] + introns[1] + exons_cds[2]
        strand = "+" if k % 2 else "-"
        start = pos
        seq[start:start + len(gene)] = list(gene if strand == "+" else revcomp(gene))
        end = start + len(gene)
        rows.append((name, "synthetic", "gene", start + 1, end, ".", strand, ".", f"ID={gid}"))
        rows.append((name, "synthetic", "mRNA", start + 1, end, ".", strand, ".", f"ID={gid}_T001;Parent={gid}"))
        # exon coordinates on the + strand of the gene, then mapped to the chromosome
        offs, o = [], 0
        for i, e in enumerate(exons_cds):
            offs.append((o, o + len(e)))
            o += len(e) + (len(introns[i]) if i < 2 else 0)
        phase = 0
        for i, (a, b) in enumerate(offs):
            if strand == "+":
                s1, e1 = start + a + 1, start + b
            else:
                s1, e1 = end - b + 1, end - a
            rows.append((name, "synthetic", "exon", s1, e1, ".", strand, ".", f"Parent={gid}_T001"))
            rows.append((name, "synthetic", "CDS", s1, e1, ".", strand, str(phase),
                         f"ID=CDS:{gid}_T001;Parent={gid}_T001"))
            phase = (3 - ((b - a - phase) % 3)) % 3
        pos = end + rng.randint(4000, 6000)
    return "".join(seq), rows


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def write_fai(path, records, width=60):
    off = 0
    with open(path + ".fai", "w") as fh:
        for name, s in records:
            off += len(name) + 2
            fh.write(f"{name}\t{len(s)}\t{off}\t{width}\t{width + 1}\n")
            off += len(s) + (len(s) + width - 1) // width


def build(root, seed=7):
    rng = random.Random(seed)
    chrA, rowsA = make_chrom("chrA", 300_000, rng, "GA")
    chrB, rowsB = make_chrom("chrB", 200_000, rng, "GB")
    ref = [("chrA", chrA), ("chrB", chrB)]
    ref_fa = os.path.join(root, "ref.fa")
    write_fasta(ref_fa, ref)
    write_fai(ref_fa, ref)
    gff = os.path.join(root, "ref.gff3.gz")
    with gzip.open(gff, "wt") as fh:
        fh.write("##gff-version 3\n")
        for r in rowsA + rowsB:
            fh.write("\t".join(map(str, r)) + "\n")
    m = lambda s, r=0.005: mutate(s, r, rng)
    ctg2 = chrA[150_000:200_000] + revcomp(chrA[200_000:260_000]) + chrA[260_000:]
    asm = [("ctg1", m(chrA[:150_000])), ("ctg2", m(ctg2)), ("ctg3", m(chrB)),
           ("ctg4", m(chrB[50_000:120_000], 0.01)), ("ctg5", m(chrA[:60_000] + chrB[140_000:]))]
    asm_fa = os.path.join(root, "asm.fa")
    write_fasta(asm_fa, asm)
    write_fai(asm_fa, asm)
    return {"ref_fa": ref_fa, "gff": gff, "asm_fa": asm_fa,
            "genes": len([r for r in rowsA + rowsB if r[2] == "gene"])}
