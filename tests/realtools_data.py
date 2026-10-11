"""Tiny synthetic sequencing data for real-tool compatibility tests.

A random 300 kb "genome" with error-free HiFi-like reads (~30x, 15 kb) and
ONT-like reads (~5x, 60-90 kb). Small enough for hifiasm, meryl, Merqury and
BUSCO to run in seconds to minutes on a few threads.
"""
import gzip
import hashlib
import os
import random
import subprocess
import uuid

COMP = str.maketrans("ACGT", "TGCA")


def genome(n=300_000, seed=11):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def sample_reads(g, n_reads, length, seed, jitter=0.2):
    rng = random.Random(seed)
    out = []
    for _ in range(n_reads):
        L = int(length * (1 + rng.uniform(-jitter, jitter)))
        L = min(L, len(g))
        start = rng.randrange(0, len(g) - L + 1)
        s = g[start:start + L]
        if rng.random() < 0.5:
            s = s.translate(COMP)[::-1]
        out.append(s)
    return out


def write_fastq(path, names, seqs, q="I"):
    with gzip.open(path, "wt", compresslevel=1) as fh:
        for n, s in zip(names, seqs):
            fh.write(f"@{n}\n{s}\n+\n{q * len(s)}\n")


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


def build(outdir, samtools=None, movie="m84109_240101_000000_s1",
          bam_movie="m84109_240202_000000_s2"):
    """hifi.fastq.gz (movie A, 500 reads) + hifi_subset.bam (movie B, 100 reads)
    together give ~30x; ont.fastq.gz has 10 reads of ~20 kb and 25 of 60-90 kb."""
    os.makedirs(outdir, exist_ok=True)
    g = genome()
    write_fasta(os.path.join(outdir, "genome.fa"), "chrSynthetic", g)
    hifi = sample_reads(g, 600, 15000, seed=1)
    names = [f"{movie}/{z}/ccs" for z in range(1, 501)]
    write_fastq(os.path.join(outdir, "hifi.fastq.gz"), names, hifi[100:])
    rng = random.Random(3)
    ont = sample_reads(g, 10, 20000, seed=4) + sample_reads(g, 25, 75000, seed=2)
    onames = [f"{uuid.UUID(int=rng.getrandbits(128))} runid=synthrun flow_cell_id=FAKE00000 "
              f"basecall_model_version_id=synthetic" for _ in ont]
    write_fastq(os.path.join(outdir, "ont.fastq.gz"), onames, ont, q="5")
    if samtools:
        # PacBio read-group IDs are the first 8 hex digits of md5("<movie>//<READTYPE>").
        rg = hashlib.md5(f"{bam_movie}//CCS".encode()).hexdigest()[:8]
        sam = ["@HD\tVN:1.6\tSO:unknown\tpb:5.0.0",
               f"@RG\tID:{rg}\tPL:PACBIO\tPU:{bam_movie}\tDS:READTYPE=CCS;BINDINGKIT=101-894-200;"
               "SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000\tPM:REVIO"]
        bam_names = [f"{bam_movie}/{z}/ccs" for z in range(1, 101)]
        for z, (n, s) in enumerate(zip(bam_names, hifi[:100]), 1):
            sam.append(f"{n}\t4\t*\t0\t255\t*\t*\t0\t0\t{s}\t{'I' * len(s)}\tRG:Z:{rg}"
                       f"\tzm:i:{z}\tnp:i:10\trq:f:0.999")
        subprocess.run([samtools, "view", "-b", "-o", os.path.join(outdir, "hifi_subset.bam"), "-"],
                       input="\n".join(sam) + "\n", text=True, check=True)
    return outdir


if __name__ == "__main__":
    import sys
    build(sys.argv[1], samtools=sys.argv[2] if len(sys.argv) > 2 else None)
