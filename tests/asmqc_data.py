"""Synthetic genome, reads and two 'assemblies' with planted structural
features for the assembly-validation tests. All sequence is random; nothing
comes from real data.

True diploid genome (haplotypes 1 and 2 differ only in chrB:100-140 kb, 1.5%
substitutions):
  chrA 200 kb; a 6 kb repeat R at 50 kb and again at 120 kb
  chrB 200 kb
  chrC 100 kb

Assembly A ("unpurged"):
  ctg1  chrA with the second R copy removed       -> collapsed repeat at ~50 kb (2x depth)
  ctg2  chrB[0:150k] + chrC[50k:100k]             -> misjoin at 150 kb (clipped reads, none spanning)
  ctg3  haplotype-2 chrB[100k:140k]               -> retained allele (half depth, unique)
  ctg4  chrC[60k:90k] identical copy              -> near-identical duplicate (half depth, MAPQ 0)
  ctg5  chrB[150k:200k]
  ctg6  chrC[0:50k]
Assembly B ("purged"): A without ctg3 and ctg4.
"""
import gzip
import os
import random

R_LEN = 6000


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


def write_fa(path, recs, width=80):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def build(d, seed=7, depth_per_hap=15):
    rng = random.Random(seed)
    rep = rseq(R_LEN, rng)
    chrA = rseq(50000, rng) + rep + rseq(120000 - 50000 - R_LEN, rng) + rep + rseq(200000 - 120000 - R_LEN, rng)
    chrB = rseq(200000, rng)
    chrC = rseq(100000, rng)
    chrB2 = chrB[:100000] + mutate(chrB[100000:140000], 0.015, rng) + chrB[140000:]
    hap1 = {"chrA": chrA, "chrB": chrB, "chrC": chrC}
    hap2 = {"chrA": chrA, "chrB": chrB2, "chrC": chrC}
    ctgs = {
        "ctg1": chrA[:120000] + chrA[120000 + R_LEN:],
        "ctg2": chrB[:150000] + chrC[50000:],
        "ctg3": chrB2[100000:140000],
        "ctg4": chrC[60000:90000],
        "ctg5": chrB[150000:],
        "ctg6": chrC[:50000],
    }
    os.makedirs(d, exist_ok=True)
    paths = {"ref": os.path.join(d, "ref.fa"), "asm_a": os.path.join(d, "asm_a.fa"),
             "asm_b": os.path.join(d, "asm_b.fa"), "reads": os.path.join(d, "reads.fastq.gz")}
    write_fa(paths["ref"], list(hap1.items()))
    write_fa(paths["asm_a"], list(ctgs.items()))
    write_fa(paths["asm_b"], [(k, v) for k, v in ctgs.items() if k not in ("ctg3", "ctg4")])
    n = 0
    with gzip.open(paths["reads"], "wt") as fh:
        for hi, hap in enumerate((hap1, hap2)):
            for name, seq in hap.items():
                L = len(seq)
                for _ in range(int(depth_per_hap * L / 10000)):
                    ln = rng.randint(8000, 12000)
                    s = rng.randint(-ln // 2, L - ln // 2)
                    s, e = max(0, s), min(L, s + ln)
                    if e - s < 3000:
                        continue
                    r = mutate(seq[s:e], 0.001, rng)
                    if rng.random() < 0.5:
                        r = revcomp(r)
                    n += 1
                    fh.write(f"@m00000_000000_000000/{n}/ccs\n{r}\n+\n{'I' * len(r)}\n")
    paths["n_reads"] = n
    return paths


SCIENCE_OVERRIDES = {
    "read_support": {"window": 1000, "chunk": 50000, "end_margin": 3000},
    "ref_align": {"min_block_len": 2000, "plot_min_len": 1000, "synteny_merge_gap": 20000,
                  "min_segment_len": 20000, "multi_chr_min_bp": 20000, "ref_jump": 50000,
                  "coverage_bin": 10000, "min_absent_bins": 2},
    "pair_align": {"min_aln_len": 1000},
    "evidence": {"junction_flank": 1000},
    "graph": {"top_candidates": 5},
}
