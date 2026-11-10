#!/usr/bin/env python3
"""Reversible short sequence-ID aliases for tools that rewrite long IDs (e.g. EDTA, which truncates IDs longer than
13 characters and replaces '-' and other punctuation; its output keeps the rewritten IDs).

Subcommands:
  make       --fasta IN.fa[.gz] --code CODE --out-fasta OUT.fa.gz --out-map MAP.tsv [--max-len 13]
             [--chrom-regex REGEX] [--threads 8]
             Alias of each record: CODE_chrNN when the ID matches --chrom-regex (default '_chr(\\d+)$'), otherwise
             CODE_<k><ordinal> (k = 'c' for a file without chromosome records, 'u' for the other records of a file
             that has them), zero-padded to the record count. Aliases must be unique, at most --max-len characters
             and only [A-Za-z0-9_]; anything else is refused. Sequence lines are copied byte for byte; only header lines
             change (any description after the ID is dropped). MAP.tsv: alias, original_id, length, seq_sha256_upper.
  translate  --map MAP.tsv --in FILE --out FILE --format {gff3,bed,rmout,tsv} [--column N] [--to original|alias]
             [--allow-unknown]
             Rewrite sequence IDs (default alias -> original). gff3: column 1 and ##sequence-region; bed / tsv: column 1
             (tsv: --column, 1-based); rmout: RepeatMasker .out query column (5th field), spacing kept. An ID that is
             not in the map stops the run unless --allow-unknown. Attributes are never touched.

Verify a made alias FASTA with g2f_seqid_alias_verify.py (separate code). Every output is written as OUT.partial and
renamed on success; existing outputs are never overwritten.
"""
import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys

SAFE = re.compile(r"^[A-Za-z0-9_]+$")


def die(msg):
    sys.exit(f"g2f_seqid_alias: {msg}")


def refuse_existing(*paths):
    for p in paths:
        if os.path.exists(p) or os.path.exists(p + ".partial"):
            die(f"refusing to overwrite {p}")


def open_text_lines(path, bgzip):
    if path.endswith(".gz"):
        p = subprocess.Popen([bgzip, "-dc", path], stdout=subprocess.PIPE)
        return p, p.stdout
    return None, open(path, "rb")


def headers(path, bgzip):
    p, fh = open_text_lines(path, bgzip)
    ids = [l[1:].split()[0].decode() for l in fh if l.startswith(b">")]
    if p and p.wait() != 0:
        die(f"cannot read {path}")
    return ids


def make_aliases(ids, code, max_len, chrom_regex):
    if not SAFE.match(code):
        die(f"--code {code!r} must be [A-Za-z0-9_] only")
    rx = re.compile(chrom_regex)
    chrom = {i: rx.search(i) for i in ids}
    has_chrom = any(chrom.values())
    other = [i for i in ids if not chrom[i]]
    width = max(4, len(str(len(other))))
    k = "u" if has_chrom else "c"
    ordinal = {i: n for n, i in enumerate(other, 1)}
    al = {}
    for i in ids:
        al[i] = f"{code}_chr{chrom[i].group(1)}" if chrom[i] else f"{code}_{k}{ordinal[i]:0{width}d}"
    bad = [a for a in al.values() if len(a) > max_len or not SAFE.match(a)]
    if bad:
        die(f"{len(bad)} aliases longer than {max_len} or not [A-Za-z0-9_], e.g. {bad[:3]}; use a shorter --code")
    if len(set(al.values())) != len(al) or len(set(ids)) != len(ids):
        die("aliases or source IDs are not unique")
    return al


def cmd_make(a):
    refuse_existing(a.out_fasta, a.out_map)
    ids = headers(a.fasta, a.bgzip)
    al = make_aliases(ids, a.code, a.max_len, a.chrom_regex)
    p, fh = open_text_lines(a.fasta, a.bgzip)
    part = a.out_fasta + ".partial"
    info, cur, h, n = {}, None, None, 0
    with open(part, "wb") as fo:
        enc = subprocess.Popen([a.bgzip, "-@", str(a.threads), "-c"], stdin=subprocess.PIPE, stdout=fo)
        for line in fh:
            if line.startswith(b">"):
                if cur is not None:
                    info[cur] = (n, h.hexdigest())
                cur = line[1:].split()[0].decode()
                h, n = hashlib.sha256(), 0
                enc.stdin.write(f">{al[cur]}\n".encode())
            else:
                s = line.rstrip(b"\r\n")
                h.update(s.upper())
                n += len(s)
                enc.stdin.write(line)
        if cur is not None:
            info[cur] = (n, h.hexdigest())
        enc.stdin.close()
        if enc.wait() != 0 or (p and p.wait() != 0):
            die("bgzip failed")
    with open(a.out_map + ".partial", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["alias", "original_id", "length", "seq_sha256_upper"])
        for i in ids:
            w.writerow([al[i], i, info[i][0], info[i][1]])
    os.rename(part, a.out_fasta)
    os.rename(a.out_map + ".partial", a.out_map)
    print(f"{len(ids)} records; longest alias {max(len(x) for x in al.values())} characters -> {a.out_fasta}")
    return 0


def cmd_translate(a):
    refuse_existing(a.out)
    with open(a.map) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    m = {r["alias"]: r["original_id"] for r in rows} if a.to == "original" else {r["original_id"]: r["alias"] for r in rows}
    if len(m) != len(rows):
        die("map is not one-to-one")
    unknown = set()

    def tr(x):
        if x in m:
            return m[x]
        unknown.add(x)
        return x
    col = (a.column or 1) - 1
    rm = re.compile(r"^(\s*(?:\S+\s+){4})(\S+)(.*)$", re.S)
    n = 0
    with open(a.in_) as fi, open(a.out + ".partial", "w") as fo:
        for l in fi:
            if a.format == "gff3":
                if l.startswith("##sequence-region"):
                    f = l.split()
                    f[1] = tr(f[1])
                    l = " ".join(f) + "\n"
                elif l.strip() and not l.startswith("#"):
                    f = l.split("\t", 1)
                    l = tr(f[0]) + "\t" + f[1]
                    n += 1
            elif a.format in ("bed", "tsv"):
                if l.strip() and not l.startswith(("#", "track", "browser")):
                    f = l.rstrip("\n").split("\t")
                    f[col] = tr(f[col])
                    l = "\t".join(f) + "\n"
                    n += 1
            elif a.format == "rmout":
                mm = rm.match(l)
                if mm and mm.group(1).split()[0].isdigit():
                    old = mm.group(2)
                    new = tr(old)
                    pad = max(0, len(old) - len(new))
                    l = mm.group(1) + new + " " * pad + mm.group(3)
                    n += 1
            fo.write(l)
    if unknown and not a.allow_unknown:
        os.remove(a.out + ".partial")
        die(f"{len(unknown)} sequence IDs not in the map, e.g. {sorted(unknown)[:5]}; nothing written")
    os.rename(a.out + ".partial", a.out)
    print(f"{n} records translated to {a.to} IDs -> {a.out}" + (f"; {len(unknown)} unknown IDs left as is" if unknown else ""))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("make")
    p.add_argument("--fasta", required=True)
    p.add_argument("--code", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True)
    p.add_argument("--max-len", type=int, default=13)
    p.add_argument("--chrom-regex", default=r"_chr(\d+)$")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--bgzip", default=os.environ.get("G2F_BGZIP", "bgzip"))
    p = sub.add_parser("translate")
    p.add_argument("--map", required=True)
    p.add_argument("--in", dest="in_", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--format", required=True, choices=["gff3", "bed", "rmout", "tsv"])
    p.add_argument("--column", type=int)
    p.add_argument("--to", choices=["original", "alias"], default="original")
    p.add_argument("--allow-unknown", action="store_true")
    a = ap.parse_args(argv)
    return cmd_make(a) if a.cmd == "make" else cmd_translate(a)


if __name__ == "__main__":
    sys.exit(main())
