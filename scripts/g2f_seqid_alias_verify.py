#!/usr/bin/env python3
"""Independent verification of an alias FASTA made by g2f_seqid_alias.py make (shares no code with it).

  g2f_seqid_alias_verify.py --source IN.fa[.gz] --alias OUT.fa.gz --map MAP.tsv [--max-len 13] --out verify.json

Checks, re-reading both FASTAs with Python's gzip module:
  A1 same number of records, same order; record k of the alias FASTA is named map[k].alias and record k of the source
     is map[k].original_id
  A2 every sequence identical (exact bytes of the sequence lines, and uppercase sha256 = map column)
  A3 aliases unique, at most --max-len characters, [A-Za-z0-9_] only; original IDs unique
  A4 map lengths equal the record lengths
Exit status 1 on any failure.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import sys


def recs(path):
    op = gzip.open if path.endswith(".gz") else open
    out, name, he, hu, n = [], None, None, None, 0
    with op(path, "rb") as fh:
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True)
    ap.add_argument("--alias", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--max-len", type=int, default=13)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    fails = []
    src, ali = recs(a.source), recs(a.alias)
    with open(a.map) as fh:
        m = list(csv.DictReader(fh, delimiter="\t"))
    if not (len(src) == len(ali) == len(m)):
        fails.append(f"A1 record counts differ: source {len(src)}, alias {len(ali)}, map {len(m)}")
    for s, x, r in zip(src, ali, m):
        if s[0] != r["original_id"] or x[0] != r["alias"]:
            fails.append(f"A1 {s[0]} / {x[0]} do not match map row {r['original_id']} / {r['alias']}")
        if s[1:] != x[1:] or s[3] != r["seq_sha256_upper"]:
            fails.append(f"A2 sequence of {s[0]} differs")
        if int(r["length"]) != s[1]:
            fails.append(f"A4 length of {s[0]}")
    al = [r["alias"] for r in m]
    if len(set(al)) != len(al) or len({r["original_id"] for r in m}) != len(m):
        fails.append("A3 aliases or original IDs not unique")
    bad = [x for x in al if len(x) > a.max_len or not re.match(r"^[A-Za-z0-9_]+$", x)]
    if bad:
        fails.append(f"A3 {len(bad)} unsafe aliases, e.g. {bad[:3]}")
    res = {"schema": "g2f-seqid-alias-verify/1", "source": a.source, "alias": a.alias, "map": a.map, "records": len(src),
           "max_alias_len": max((len(x) for x in al), default=0), "result": "PASS" if not fails else "FAIL",
           "failures": fails[:50]}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("records", "max_alias_len", "result")}))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
