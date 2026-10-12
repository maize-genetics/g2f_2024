"""Synthetic genotypes and reads for fingerprinting tests (small, deterministic).

A random ancestral genome is mutated into several "inbred lines" (SNPs,
deleted blocks, unique inserted blocks). One pair is closely related. Reads:
HiFi-like (0.1 % substitutions) and ONT-like (~5 % substitutions/indels).
"""
import gzip
import io
import os
import random
import tarfile

BASES = "ACGT"


def rand_seq(n, rng):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(seq, rng, snp_rate, n_del=0, n_ins=0, block=2000):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < snp_rate:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    s = "".join(s)
    for _ in range(n_del):
        p = rng.randrange(0, len(s) - block)
        s = s[:p] + s[p + block:]
    for _ in range(n_ins):
        p = rng.randrange(0, len(s))
        s = s[:p] + rand_seq(block, rng) + s[p:]
    return s


def genotypes(seed=7, length=300_000):
    """Return {name: genome}. B2 is a close relative of B1; X is outside the panel."""
    rng = random.Random(seed)
    anc = rand_seq(length, rng)
    g = {}
    g["A1"] = mutate(anc, rng, 0.004, 4, 4)
    g["B1"] = mutate(anc, rng, 0.004, 4, 4)
    g["B2"] = mutate(g["B1"], rng, 0.0015, 2, 2)
    g["C1"] = mutate(anc, rng, 0.004, 4, 4)
    g["X"] = mutate(anc, rng, 0.004, 4, 4)
    return g


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def reads(genome, depth, mean_len, rng, sub=0.001, indel=0.0):
    out, total, n = [], 0, len(genome)
    while total < depth * n:
        L = max(500, int(rng.gauss(mean_len, mean_len * 0.2)))
        p = rng.randrange(0, max(1, n - L))
        s = genome[p:p + L]
        if rng.random() < 0.5:
            s = _revcomp(s)
        if sub or indel:
            b = []
            for c in s:
                r = rng.random()
                if r < sub:
                    b.append(rng.choice([x for x in BASES if x != c]))
                elif r < sub + indel / 2:
                    continue
                elif r < sub + indel:
                    b.append(c + rng.choice(BASES))
                else:
                    b.append(c)
            s = "".join(b)
        out.append(s)
        total += len(s)
    return out


def mix(parts, rng):
    out = []
    for rs in parts:
        out += rs
    rng.shuffle(out)
    return out


def write_fastq_gz(path, seqs, prefix):
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@{prefix}_{i}\n{s}\n+\n{'5' * len(s)}\n")
    return path


def write_fasta_gz(path, seqs, prefix):
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">{prefix}/{i}/ccs\n" + "\n".join(s[j:j + 80] for j in range(0, len(s), 80)) + "\n")
    return path


def write_ont_tar(path, seqs, prefix, batch=200, fail=None):
    """MinKNOW-like tar: fastq_pass/<..>_<n>.fastq.gz members (+ optional fastq_fail)."""
    with tarfile.open(path, "w") as tf:
        def add(name, data):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
        for b, i in enumerate(range(0, len(seqs), batch)):
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
                for j, s in enumerate(seqs[i:i + batch]):
                    gz.write(f"@{prefix}_{i + j} runid=r1 sample_id={prefix}\n{s}\n+\n"
                             f"{'5' * len(s)}\n".encode())
            add(f"./fastq_pass/{prefix}_pass_{b}.fastq.gz", buf.getvalue())
        if fail:
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
                for j, s in enumerate(fail):
                    gz.write(f"@{prefix}_fail_{j}\n{s}\n+\n{'!' * len(s)}\n".encode())
            add(f"./fastq_fail/{prefix}_fail_0.fastq.gz", buf.getvalue())
    return path


def build(outdir, samtools=None, seed=11):
    """Write the synthetic dataset. Returns a list of (label, tech, run_id, path, truth)."""
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    g = genotypes()
    items = []

    def hifi(label, run, depth, fmt="fastq"):
        rs = reads(g[label], depth, 8000, rng, sub=0.001)
        p = os.path.join(outdir, f"{label}.{run}.hifi.{'fasta' if fmt == 'fasta' else 'fastq'}.gz")
        (write_fasta_gz if fmt == "fasta" else write_fastq_gz)(p, rs, f"m0_{label}_{run}")
        if fmt == "bam" and samtools:
            import subprocess
            bam = os.path.join(outdir, f"{label}.{run}.hifi.bam")
            sam = ["@HD\tVN:1.6", f"@RG\tID:{run}\tPL:PACBIO\tDS:READTYPE=CCS"]
            sam += [f"m84000_240101_000000_s1/{i}/ccs\t4\t*\t0\t255\t*\t*\t0\t0\t{s}\t{'I' * len(s)}"
                    f"\tRG:Z:{run}" for i, s in enumerate(rs)]
            subprocess.run([samtools, "view", "-b", "-o", bam, "-"], input="\n".join(sam) + "\n",
                           text=True, check=True)
            os.remove(p)
            p = bam
        items.append((label, "hifi", run, p, label))

    hifi("A1", "r1", 8)
    hifi("A1", "r2", 8, fmt="bam")
    hifi("B1", "r1", 8)
    hifi("B1", "r2", 8)
    hifi("B2", "r1", 10, fmt="fasta")
    hifi("C1", "r1", 10)

    def ont(label, run, truth_parts, depth, fmt="tar", sub=0.03, indel=0.03):
        rs = mix([reads(g[t], depth * f, 6000, rng, sub=sub, indel=indel)
                  for t, f in truth_parts], rng)
        if fmt == "tar":
            p = write_ont_tar(os.path.join(outdir, f"{run}_fastq.tar"), rs, run,
                              fail=reads(g["C1"], 0.5, 3000, rng, sub=0.1))
        else:
            p = write_fastq_gz(os.path.join(outdir, f"{run}.ont.fastq.gz"), rs, run)
        items.append((label, "ont", run, p, "+".join(t for t, _ in truth_parts)))

    ont("A1", "ontA1", [("A1", 1.0)], 4)                    # correct
    ont("B1", "ontB1", [("B1", 1.0)], 4, fmt="fastq")       # correct, close relative B2 in panel
    ont("C1", "ontC1_mislabel", [("A1", 1.0)], 4)           # WRONG label: reads are A1
    ont("C1", "ontC1_ok", [("C1", 1.0)], 3)                 # correct
    ont("B2", "ontB2_mix", [("B2", 0.5), ("C1", 0.5)], 4)   # mixture
    ont("A1", "ontX", [("X", 1.0)], 4)                      # genotype not in panel
    ont("B1", "ontB1_tiny", [("B1", 1.0)], 0.02)            # too little data
    return items
