"""Synthetic fixture builders for tests (tiny files, no real data)."""
import copy
import gzip
import os
import random
import shutil
import subprocess
import sys
import uuid

PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PIPELINE_DIR, "scripts"))

from g2f.config import load_config  # noqa: E402

SAMTOOLS = shutil.which("samtools") or "/programs/samtools-1.20/bin/samtools"
SEQKIT = shutil.which("seqkit") or "/programs/seqkit-2.13.0/seqkit"
HAVE_SAMTOOLS = os.path.exists(SAMTOOLS)
HAVE_SEQKIT = os.path.exists(SEQKIT)


def write_bytes(path, data):
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def hifi_names(movie, n, start=1):
    return [f"{movie}/{z}/ccs" for z in range(start, start + n)]


def write_fastq(path, names, lens, seed=1, gz=True, qchar="I"):
    rng = random.Random(seed)
    op = gzip.open if gz else open
    with op(path, "wt") as fh:
        for name, n in zip(names, lens):
            fh.write(f"@{name}\n{rand_seq(n, rng)}\n+\n{qchar * n}\n")
    return path


def write_fasta(path, names, lens, seed=1, gz=True):
    rng = random.Random(seed)
    op = gzip.open if gz else open
    with op(path, "wt") as fh:
        for name, n in zip(names, lens):
            s = rand_seq(n, rng)
            fh.write(f">{name}\n" + "\n".join(s[i:i + 60] for i in range(0, n, 60)) + "\n")
    return path


def write_bam(path, names, lens, readtype="CCS", platform="PACBIO", movie="m84109_240101_000000_s1",
              seed=1, extra_rg=""):
    """Write an unaligned BAM via samtools (real BGZF)."""
    rng = random.Random(seed)
    rg = f"@RG\tID:{movie[-6:]}\tPL:{platform}\tPU:{movie}\tDS:READTYPE={readtype};Ipd:Frames{extra_rg}"
    sam = ["@HD\tVN:1.6\tSO:unknown", rg]
    for name, n in zip(names, lens):
        sam.append(f"{name}\t4\t*\t0\t255\t*\t*\t0\t0\t{rand_seq(n, rng)}\t{'I' * n}\tRG:Z:{movie[-6:]}")
    subprocess.run([SAMTOOLS, "view", "-b", "-o", path, "-"], input="\n".join(sam) + "\n",
                   text=True, check=True)
    return path


def ont_names(n, runid="0123456789abcdef0123456789abcdef01234567", seed=1):
    rng = random.Random(seed)
    return [f"{uuid.UUID(int=rng.getrandbits(128))} runid={runid} read={i} ch=1 "
            f"flow_cell_id=PAQ00000 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
            for i in range(n)]


BASE_TOP = {
    "samples": "config/samples.example.yaml",
    "external_assemblies": "config/external_assemblies.example.yaml",
    "science": "config/science.yaml",
    "methods_dir": "config/methods",
    "resources": "config/resources/default.yaml",
    "tools": "config/tools/default.yaml",
    "legacy_dirs": [],
}


def make_cfg(tmp, samples=None, external=None, science_overrides=None):
    """Load the real tracked configuration, pointed at a temporary data/results tree."""
    top = dict(BASE_TOP, data_dir=os.path.join(tmp, "data"), results_dir=os.path.join(tmp, "results"),
               local=os.path.join(tmp, "local.yaml"))
    cfg = load_config(top, PIPELINE_DIR)
    cfg["samples"] = copy.deepcopy(samples or {})
    cfg["local"] = {"external_inputs": external or {}}
    for k, v in (science_overrides or {}).items():
        cfg["science"][k] = v
    return cfg


def mkdirs(tmp, sample, techs=("hifi",)):
    out = {}
    for t in techs:
        d = os.path.join(tmp, "data", sample, t)
        os.makedirs(d, exist_ok=True)
        out[t] = d
    return out
