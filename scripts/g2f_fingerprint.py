#!/usr/bin/env python3
"""Optional genotype fingerprinting (k-mer based sample-identity QC).

Two steps, usable inside the Snakemake workflow (`fingerprint` target) or on
their own with any read files:

  # 1. sketch each read file (needs the sourmash library; one job per file)
  g2f_fingerprint.py sketch --input reads.bam --tech hifi --label LINE_A \\
      --run-id m84109_... --out-sig A.sig.gz --out-json A.json \\
      [--ksizes 21,31 --scaled 2000 --max-bases 15000000000]

  # 2. compare all ONT read sets with all HiFi genotypes (standard library only)
  g2f_fingerprint.py analyze --records sketches/*.json --outdir fingerprint/ \\
      [--ksizes 21,31 --primary-k 21]

--label is the genotype name. HiFi labels define the reference panel. For ONT
it is the EXPECTED genotype being tested; leave it empty ("") if unknown.
The analysis reports evidence only; it never moves or relabels files.
See docs/fingerprinting.md.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _ints(s):
    return [int(x) for x in str(s).split(",") if x.strip()]


def cmd_sketch(a):
    from g2f import fpsketch
    for p in (a.out_sig, a.out_json):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    rec = fpsketch.sketch(a.input, a.out_sig, a.out_json, a.label, a.tech, a.run_id,
                          ksizes=_ints(a.ksizes), scaled=a.scaled, max_bases=a.max_bases,
                          read_fraction=a.read_fraction, min_read_len=a.min_read_len,
                          seed=a.seed, fmt=a.format, samtools=a.samtools, threads=a.threads)
    s = rec["sampling"]
    sys.stderr.write(f"sketched {a.input}: {s['reads_used']} reads, {s['bases_used']} bases "
                     f"in {rec['runtime_s']} s; hashes {rec['hashes']}\n")


def cmd_analyze(a):
    from g2f import fpanalysis
    params = json.loads(a.params) if a.params else {}
    prov = json.loads(a.provenance) if a.provenance else {}
    ks = _ints(a.ksizes)
    primary = a.primary_k or ks[0]
    if primary not in ks:
        sys.exit("--primary-k must be one of --ksizes")
    recs = list(a.records)
    if a.records_from:
        with open(a.records_from) as fh:
            recs += [line.strip() for line in fh if line.strip()]
    fpanalysis.analyze(recs, a.outdir, ks, primary, params, provenance=prov,
                       log=lambda m: (sys.stderr.write(m + "\n"), sys.stderr.flush()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sketch", help="FracMinHash sketch of one read file (streaming)")
    s.add_argument("--input", required=True)
    s.add_argument("--tech", choices=["hifi", "ont"], required=True)
    s.add_argument("--label", required=True, help="genotype (ONT: expected genotype, or '')")
    s.add_argument("--run-id", required=True)
    s.add_argument("--out-sig", required=True)
    s.add_argument("--out-json", required=True)
    s.add_argument("--format", choices=["bam", "fastq", "fasta", "tar"])
    s.add_argument("--ksizes", default="21")
    s.add_argument("--scaled", type=int, default=2000)
    s.add_argument("--max-bases", type=int, default=0, help="0 = whole file")
    s.add_argument("--read-fraction", type=float, default=1.0)
    s.add_argument("--min-read-len", type=int, default=0)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--samtools", default="samtools")
    s.add_argument("--threads", type=int, default=2)
    s.set_defaults(func=cmd_sketch)

    s = sub.add_parser("analyze", help="compare ONT read sets with the HiFi genotype panel")
    s.add_argument("--records", nargs="*", default=[], help="sketch JSON records")
    s.add_argument("--records-from", help="file listing sketch JSON records, one per line")
    s.add_argument("--outdir", required=True)
    s.add_argument("--ksizes", default="21")
    s.add_argument("--primary-k", type=int)
    s.add_argument("--params", help="JSON object of analysis parameters")
    s.add_argument("--provenance", help="JSON object recorded in parameters_and_provenance.json")
    s.set_defaults(func=cmd_analyze)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
