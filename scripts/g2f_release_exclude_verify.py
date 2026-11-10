#!/usr/bin/env python3
"""Independent verification of a g2f_release_exclude.py build (shares no code with the builder).

Re-reads the source and output FASTAs with Python's gzip module and checks:
  X1  contigs: output IDs = source IDs minus the plan, same order; every kept record identical (sequence sha256)
  X2  scaffolds: the same, for the scaffold FASTA
  X3  every excluded record existed in both sources, its contig and scaffold sequences are identical, and its length
      and sha256 (uppercase) match the EXCLUDED_RECORDS table; the table lists exactly the plan
  X4  AGP: output data lines = source data lines minus exactly one whole, + singleton line per excluded record;
      with --relabel-linkage OLD=NEW, gap lines carry NEW in place of OLD and nothing else changes (no OLD remains);
      comment lines = the source comments, or with --agp-header the template rendered here independently
      ({contigs_sha256} recomputed from the decompressed output contig FASTA)
  X5  contig_to_scaffold map and chain: source minus exactly the excluded records; gap BED identical
  X6  totals: bp(out) = bp(source) - excluded bp, for contigs and scaffolds
  X7  .fai of each output agrees with its records; BUILD.json output hashes match the files

  g2f_release_exclude_verify.py --plan plan.tsv --unit U --build DIR --src-contigs C.fa.gz --src-scaffolds S.fa.gz
      --src-agp X.agp --src-chain X.chain --src-map X.tsv --src-gaps-bed X.bed --out verify.json
      [--agp-header TEMPLATE] [--relabel-linkage OLD=NEW]   (the same values given to the builder)
Exit status 1 on any failure.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import sys


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def records(path):
    """Ordered list of (id, length, sha256_exact, sha256_upper)."""
    out = []
    name = he = hu = None
    n = 0
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if name is not None:
                    out.append((name, n, he.hexdigest(), hu.hexdigest()))
                name = line[1:].split()[0].decode()
                he, hu, n = hashlib.sha256(), hashlib.sha256(), 0
            else:
                s = line.rstrip(b"\r\n")
                he.update(s)
                hu.update(s.upper())
                n += len(s)
    if name is not None:
        out.append((name, n, he.hexdigest(), hu.hexdigest()))
    return out


def lines(p):
    with open(p) as fh:
        return fh.readlines()


def chain_blocks(p):
    with open(p) as fh:
        return [b.strip("\n") for b in fh.read().split("\n\n") if b.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("plan", "unit", "build", "src-contigs", "src-scaffolds", "src-agp", "src-chain", "src-map",
              "src-gaps-bed", "out"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--agp-header")
    ap.add_argument("--relabel-linkage")
    a = ap.parse_args(argv)
    checks, fails = {}, []

    def chk(key, ok, msg=""):
        checks.setdefault(key, True)
        if not ok:
            checks[key] = False
            fails.append(f"{key}: {msg}")

    with open(a.plan) as fh:
        plan = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t") if r["unit"] == a.unit]
    drop = [r["record"] for r in plan]
    dset = set(drop)
    chk("X0_plan", len(dset) == len(drop) and len(drop) > 0, "empty or duplicate plan")
    build = json.load(open(os.path.join(a.build, "BUILD.json")))
    outs = {os.path.dirname(k): os.path.join(a.build, k) for k in build["outputs"] if k.endswith(".fa.gz")}
    c_out, s_out = outs["contigs"], outs["scaffolds_rg"]

    excl = {}
    for fa_src, fa_out, key in ((a.src_contigs, c_out, "X1_contigs"), (a.src_scaffolds, s_out, "X2_scaffolds")):
        src, out = records(fa_src), records(fa_out)
        ids_src = [r[0] for r in src]
        chk(key, len(set(ids_src)) == len(ids_src), "duplicate IDs in source")
        exp = [r for r in src if r[0] not in dset]
        chk(key, [r[0] for r in out] == [r[0] for r in exp], "output IDs are not source minus plan, in order")
        chk(key, out == exp, "a kept record differs from its source")
        gone = {r[0]: r for r in src if r[0] in dset}
        chk("X3_excluded", set(gone) == dset, f"{key}: planned records missing in source: {sorted(dset - set(gone))}")
        excl[key] = gone
        chk("X6_totals", sum(r[1] for r in out) == sum(r[1] for r in src) - sum(r[1] for r in gone.values()), key)
        fai = [l.split("\t")[:2] for l in lines(fa_out + ".fai")]
        chk("X7_index", fai == [[r[0], str(r[1])] for r in out], f"{fa_out}.fai disagrees with records")
    chk("X3_excluded", excl["X1_contigs"] == excl["X2_scaffolds"], "excluded contig and scaffold records differ")
    ex_tab = os.path.join(a.build, [k for k in build["outputs"] if k.endswith("EXCLUDED_RECORDS.tsv")][0])
    with open(ex_tab) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    chk("X3_excluded", [r["record"] for r in rows] == drop, "EXCLUDED_RECORDS does not list exactly the plan")
    for r in rows:
        g = excl["X1_contigs"].get(r["record"])
        chk("X3_excluded", g is not None and (int(r["length"]), r["seq_sha256_upper"]) == (g[1], g[3]),
            f"{r['record']}: length/hash differ from source")

    P = lambda src: os.path.join(a.build, "provenance", os.path.basename(src))
    agp_src, agp_out = lines(a.src_agp), lines(P(a.src_agp))
    removed = [l for l in agp_src if not l.startswith("#") and l.split("\t", 1)[0] in dset]
    exp_data = [l for l in agp_src if not l.startswith("#") and l not in removed]
    old = new = None
    if a.relabel_linkage:
        old, new = a.relabel_linkage.split("=")
        relab = []
        for l in exp_data:
            f = l.rstrip("\n").split("\t")
            if f[4] in ("N", "U"):
                f[8] = ";".join(new if t == old else t for t in f[8].split(";"))
            relab.append("\t".join(f) + "\n")
        exp_data = relab
    if a.agp_header:
        dec = hashlib.sha256()
        with gzip.open(c_out, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 24), b""):
                dec.update(b)
        t = open(a.agp_header).read()
        for k_, v_ in (("{contigs_name}", os.path.basename(c_out)), ("{contigs_sha256}", dec.hexdigest()),
                       ("{scaffolds_name}", os.path.basename(s_out)), ("{excluded}", str(len(drop)))):
            t = t.replace(k_, v_)
        exp_head = [l + "\n" for l in t.rstrip("\n").split("\n")]
    else:
        exp_head = [l for l in agp_src if l.startswith("#")]
    chk("X4_agp", [l for l in agp_out if l.startswith("#")] == exp_head, "AGP comment lines differ from the expected header")
    chk("X4_agp", [l for l in agp_out if not l.startswith("#")] == exp_data,
        "AGP data lines are not source minus excluded objects (with the declared relabel only)")
    chk("X4_agp", agp_out[:len(exp_head)] == exp_head, "AGP comments are not at the top")
    if old:
        left = [l for l in agp_out if not l.startswith("#") and l.split("\t")[4] in ("N", "U")
                and old in l.rstrip("\n").split("\t")[8].split(";")]
        chk("X4_agp", not left, f"{len(left)} gap lines still carry {old}")
    chk("X4_agp", len(removed) == len(drop), "an excluded object has more than one AGP line")
    for l in removed:
        f = l.rstrip("\n").split("\t")
        L = excl["X1_contigs"][f[0]][1]
        chk("X4_agp", f[5] == f[0] and f[4] not in ("N", "U") and f[1:3] == ["1", str(L)] and f[6:9] == ["1", str(L), "+"],
            f"{f[0]}: not a whole + singleton")
    m_src, m_out = lines(a.src_map), lines(P(a.src_map))
    chk("X5_tables", m_out == [l for l in m_src if l.split("\t", 1)[0] not in dset]
        and len(m_src) - len(m_out) == len(drop), "contig_to_scaffold map")
    c_src, c_outb = chain_blocks(a.src_chain), chain_blocks(P(a.src_chain))
    chk("X5_tables", c_outb == [b for b in c_src if b.split()[2] not in dset]
        and len(c_src) - len(c_outb) == len(drop), "chain")
    chk("X5_tables", lines(a.src_gaps_bed) == lines(P(a.src_gaps_bed)), "gap BED changed")
    for k, h in build["outputs"].items():
        chk("X7_index", sha256_file(os.path.join(a.build, k)) == h, f"BUILD.json hash mismatch for {k}")

    res = {"schema": "g2f-release-exclude-verify/1", "unit": a.unit, "build": os.path.abspath(a.build),
           "result": "PASS" if not fails else "FAIL", "checks": {k: "PASS" if v else "FAIL" for k, v in checks.items()},
           "failures": fails, "excluded": len(drop),
           "excluded_bp": sum(r[1] for r in excl["X1_contigs"].values()),
           "inputs": {k: sha256_file(getattr(a, k)) for k in ("src_contigs", "src_scaffolds", "src_agp", "src_chain",
                                                               "src_map", "src_gaps_bed", "plan")}}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("unit", "result", "checks", "failures")}))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
