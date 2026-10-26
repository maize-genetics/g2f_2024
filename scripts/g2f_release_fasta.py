#!/usr/bin/env python3
"""Reorder FASTA records for an annotation release, and verify a reordered FASTA independently.

Only the ORDER of complete records changes: each record (its header line and all of its sequence lines, with the
original line wrapping, case and bases) is copied byte for byte from the source FASTA.

Sort key for identifiers of the form <PREFIX>_chr<N>_ctg<M> / <PREFIX>_chrUn_ctg<M> (the G2F v1.1 convention):
numbered chromosomes in natural numeric order (chr1 < chr2 < ... < chr10), then ctg number numerically, then
chrUn records by ctg number. The ctg number is the existing identifier; it is a B73-anchored rank (chromosome
tier) or a length rank (chrUn), not a physical order or a scaffold. Identifiers that do not match the
convention are refused, unless --unmatched last is given (then they keep their source order after chrUn).

  g2f_release_fasta.py sort   --in SRC.fa --out DST.fa [--order-tsv ORDER.tsv]
  g2f_release_fasta.py verify --src SRC.fa --dst DST.fa [--expect-seq-sha256 MAP.tsv:id_col:sha_col] [--json OUT]

`verify` re-reads both files from scratch (it shares no state with `sort`) and checks: identical ID multiset,
each ID exactly once, byte-identical record blocks, identical per-record sequence SHA-256 and length, identical
line layout, total bp, and that DST is in the expected sort order. Exit status 0 only if every check passes.
"""
import argparse
import collections
import csv
import hashlib
import json
import os
import re
import sys

ID_RE = re.compile(r"^(?P<prefix>.+)_chr(?P<chrom>\d+|Un)_ctg(?P<ctg>\d+)$")


def sort_key(rid, unmatched_ok=False):
    m = ID_RE.match(rid)
    if not m:
        if unmatched_ok:
            return (2, 0, 0)
        raise ValueError(f"identifier does not follow <PREFIX>_chr<N|Un>_ctg<M>: {rid!r}")
    chrom = m.group("chrom")
    return (1, 0, int(m.group("ctg"))) if chrom == "Un" else (0, int(chrom), int(m.group("ctg")))


def scan(path):
    """Record blocks in file order: [(id, header_offset, block_end_offset)]. One pass, no sequence kept."""
    recs, cur, off = [], None, 0
    with open(path, "rb") as fh:
        for ln in fh:
            if ln.startswith(b">"):
                if cur:
                    recs.append((cur[0], cur[1], off))
                cur = (ln[1:].split()[0].decode(), off)
            elif cur is None and ln.strip():
                raise ValueError(f"{path}: sequence data before the first header")
            off += len(ln)
    if cur:
        recs.append((cur[0], cur[1], off))
    return recs


def cmd_sort(a):
    recs = scan(a.inp)
    ids = [r[0] for r in recs]
    dup = [k for k, v in collections.Counter(ids).items() if v > 1]
    if dup:
        sys.exit(f"duplicate identifiers in {a.inp}: {dup[:5]}")
    unmatched_ok = a.unmatched == "last"
    order = sorted(range(len(recs)), key=lambda i: (sort_key(recs[i][0], unmatched_ok), i))
    tmp = a.out + ".tmp"
    with open(a.inp, "rb") as fi, open(tmp, "wb") as fo:
        for i in order:
            _, s, e = recs[i]
            fi.seek(s)
            left = e - s
            while left:
                b = fi.read(min(left, 1 << 24))
                if not b:
                    sys.exit(f"{a.inp}: short read")
                fo.write(b)
                left -= len(b)
    os.replace(tmp, a.out)
    if a.order_tsv:
        with open(a.order_tsv, "w") as fh:
            fh.write("new_index\tid\tsource_index\tgroup\n")
            for n, i in enumerate(order, 1):
                m = ID_RE.match(recs[i][0])
                grp = f"chr{int(m.group('chrom'))}" if m and m.group("chrom") != "Un" else ("chrUn" if m else "unmatched")
                fh.write(f"{n}\t{recs[i][0]}\t{i + 1}\t{grp}\n")
    print(f"wrote {a.out}: {len(order)} records reordered")


def digests(path):
    """{id: dict(block_sha, seq_sha, length, layout)} plus file order. Independent single pass."""
    out, order, cur = {}, [], None

    def close():
        if cur is not None:
            out.setdefault(cur["id"], []).append(dict(block_sha=cur["bh"].hexdigest(), seq_sha=cur["sh"].hexdigest(),
                                                     length=cur["n"], layout=tuple(cur["widths"][:2]) + (len(cur["widths"]),)))

    with open(path, "rb") as fh:
        for ln in fh:
            if ln.startswith(b">"):
                close()
                rid = ln[1:].split()[0].decode()
                order.append(rid)
                cur = dict(id=rid, bh=hashlib.sha256(ln), sh=hashlib.sha256(), n=0, widths=[])
            else:
                cur["bh"].update(ln)
                s = ln.rstrip(b"\r\n")
                cur["sh"].update(s)
                cur["n"] += len(s)
                cur["widths"].append(len(s))
    close()
    return out, order


def cmd_verify(a):
    src, src_order = digests(a.src)
    dst, dst_order = digests(a.dst)
    chk = collections.OrderedDict()
    chk["each_src_id_exactly_once_in_src"] = all(len(v) == 1 for v in src.values())
    chk["each_src_id_exactly_once_in_dst"] = all(len(dst.get(k, [])) == 1 for k in src)
    chk["no_extra_ids_in_dst"] = set(dst) == set(src)
    chk["record_count_equal"] = len(src_order) == len(dst_order)
    common = [k for k in src if k in dst and len(src[k]) == len(dst[k]) == 1]
    chk["record_blocks_byte_identical"] = len(common) == len(src) and all(src[k][0]["block_sha"] == dst[k][0]["block_sha"] for k in common)
    chk["sequence_sha256_identical"] = len(common) == len(src) and all(src[k][0]["seq_sha"] == dst[k][0]["seq_sha"] for k in common)
    chk["lengths_identical"] = len(common) == len(src) and all(src[k][0]["length"] == dst[k][0]["length"] for k in common)
    chk["line_layout_identical"] = len(common) == len(src) and all(src[k][0]["layout"] == dst[k][0]["layout"] for k in common)
    chk["total_bp_equal"] = sum(v[0]["length"] for v in src.values()) == sum(v[0]["length"] for v in dst.values())
    exp = sorted(dst_order, key=lambda r: (sort_key(r, a.unmatched == "last"), src_order.index(r)))
    chk["dst_in_expected_sort_order"] = dst_order == exp
    if a.expect_seq_sha256:
        path, idc, shc = a.expect_seq_sha256.rsplit(":", 2)
        with open(path) as fh:
            ref = {r[idc]: r[shc] for r in csv.DictReader(fh, delimiter="\t")}
        chk["sequence_sha256_matches_reference_table"] = set(ref) == set(dst) and all(dst[k][0]["seq_sha"] == ref[k] for k in ref)
    res = dict(src=os.path.abspath(a.src), dst=os.path.abspath(a.dst), records=len(dst_order),
               total_bp=sum(v[0]["length"] for v in dst.values()), checks=chk, passed=all(chk.values()),
               first_ids=dst_order[:3], last_ids=dst_order[-3:])
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=1)
    print(json.dumps(res))
    sys.exit(0 if res["passed"] else 1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("sort")
    s.add_argument("--in", dest="inp", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--order-tsv")
    s.add_argument("--unmatched", choices=["refuse", "last"], default="refuse")
    v = sp.add_parser("verify")
    v.add_argument("--src", required=True)
    v.add_argument("--dst", required=True)
    v.add_argument("--expect-seq-sha256", help="TSV:id_column:sha256_column (sequence SHA-256, line breaks removed)")
    v.add_argument("--unmatched", choices=["refuse", "last"], default="refuse")
    v.add_argument("--json")
    a = p.parse_args()
    cmd_sort(a) if a.cmd == "sort" else cmd_verify(a)


if __name__ == "__main__":
    main()
