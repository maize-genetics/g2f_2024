#!/usr/bin/env python3
"""Independently verify a compressed FASTA against its uncompressed source, or against an AGP and its components.

Subcommands:
  source  --gz X.fa.gz --source X.fa [--source-sha256 HEX]
          The decompressed stream must be byte-identical to the source, and records (ID, length, sequence hash, order)
          must agree.
  agp     --gz S.fa.gz --agp S.agp --components C.fa[.gz] [--allow-unused]
          Every object must equal its reconstruction from the components; each component is used once and whole.

Both report record count, total bp, unique IDs and, when an index exists next to the file, its agreement with the
records. Results go to stdout as JSON (and --output). Exit status 1 on any failure. Shares no code with builders, so it
can be used as the separate verification pass required for a derived FASTA.
"""
import argparse
import gzip
import hashlib
import json
import os
import sys

BUF = 1 << 24
RC = bytes.maketrans(b"ACGTNacgtnRYKMSWBDHVrykmswbdhv", b"TGCANtgcanYRMKSWVHDByrmkswvhdb")


def opener(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def stream_sha(path):
    h = hashlib.sha256()
    with opener(path) as fh:
        for b in iter(lambda: fh.read(BUF), b""):
            h.update(b)
    return h.hexdigest()


def records(path, keep=False):
    """Yield-free reader: list of (id, length, sha256 of the sequence) and, with keep, the sequences."""
    out, seqs, name, buf = [], {}, None, []

    def close():
        if name is None:
            return
        s = b"".join(buf)
        out.append((name, len(s), hashlib.sha256(s).hexdigest()))
        if keep:
            seqs[name] = s

    with opener(path) as fh:
        for line in fh:
            if line[:1] == b">":
                close()
                name = line[1:].split()[0].decode()
                buf = []
            else:
                buf.append(line.rstrip(b"\r\n"))
    close()
    return out, seqs


def common_checks(path, recs, fails):
    ids = [r[0] for r in recs]
    if len(set(ids)) != len(ids):
        fails.append("duplicate record IDs")
    fai = path + ".fai"
    if os.path.exists(fai):
        with open(fai) as fh:
            idx = [(l.split("\t")[0], int(l.split("\t")[1])) for l in fh]
        if idx != [(r[0], r[1]) for r in recs]:
            fails.append(".fai does not match the records")
    return {"records": len(recs), "bp": sum(r[1] for r in recs)}


def cmd_source(a):
    fails = []
    gz_sha, src_sha = stream_sha(a.gz), stream_sha(a.source)
    if gz_sha != src_sha:
        fails.append(f"decompressed sha256 {gz_sha} != source {src_sha}")
    if a.source_sha256 and a.source_sha256 != src_sha:
        fails.append(f"source sha256 {src_sha} != recorded {a.source_sha256}")
    g, _ = records(a.gz)
    s, _ = records(a.source)
    if g != s:
        fails.append("records differ between compressed file and source (ID, length, sequence or order)")
    info = common_checks(a.gz, g, fails)
    return {"mode": "source", "gz": a.gz, "source": a.source, "decompressed_sha256": gz_sha, **info}, fails


def read_agp(path):
    objs = {}
    with opener(path) as fh:
        for line in fh:
            line = line.decode()
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            objs.setdefault(f[0], []).append(f)
    return objs


def cmd_agp(a):
    fails = []
    g, _ = records(a.gz)
    comp = records(a.components, keep=True)[1]
    used = {}
    objs = read_agp(a.agp)
    want = {}
    for o, rows in objs.items():
        h, n, pos = hashlib.sha256(), 0, 1
        for f in sorted(rows, key=lambda f: int(f[1])):
            b, e = int(f[1]), int(f[2])
            if b != pos:
                fails.append(f"{o}: AGP not contiguous at {pos}")
            pos = e + 1
            if f[4] in ("N", "U"):
                s = b"N" * int(f[5])
            else:
                c = comp.get(f[5])
                if c is None:
                    fails.append(f"{o}: unknown component {f[5]}")
                    continue
                if int(f[6]) != 1 or int(f[7]) != len(c):
                    fails.append(f"{f[5]}: component not used whole")
                s = c[int(f[6]) - 1:int(f[7])]
                if f[8] == "-":
                    s = s.translate(RC)[::-1]
                used[f[5]] = used.get(f[5], 0) + 1
            if len(s) != e - b + 1:
                fails.append(f"{o}:{b}-{e}: span does not match its sequence")
            h.update(s)
            n += len(s)
        want[o] = (n, h.hexdigest())
    if [r[0] for r in g] != list(objs):
        fails.append("record order or IDs differ from the AGP objects")
    for name, ln, sh in g:
        if want.get(name) != (ln, sh):
            fails.append(f"{name}: record differs from its AGP reconstruction")
    if any(v != 1 for v in used.values()):
        fails.append("a component is used more than once")
    if not a.allow_unused and set(used) != set(comp):
        fails.append(f"{len(set(comp) - set(used))} components are not used")
    info = common_checks(a.gz, g, fails)
    return {"mode": "agp", "gz": a.gz, "agp": a.agp, "components": a.components, "components_used": len(used),
            "decompressed_sha256": stream_sha(a.gz), **info}, fails


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("source")
    s.add_argument("--gz", required=True)
    s.add_argument("--source", required=True)
    s.add_argument("--source-sha256")
    s.add_argument("--output")
    g = sub.add_parser("agp")
    g.add_argument("--gz", required=True)
    g.add_argument("--agp", required=True)
    g.add_argument("--components", required=True)
    g.add_argument("--allow-unused", action="store_true")
    g.add_argument("--output")
    a = p.parse_args(argv)
    res, fails = (cmd_source if a.cmd == "source" else cmd_agp)(a)
    res = {"status": "PASS" if not fails else "FAIL", "failures": fails, **res}
    text = json.dumps(res, indent=1)
    if a.output:
        with open(a.output, "w") as fh:
            fh.write(text + "\n")
    print(text)
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
