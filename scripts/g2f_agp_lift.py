#!/usr/bin/env python3
"""Coordinate conversion between AGP components (e.g. contigs) and AGP objects (e.g. scaffolds or pseudomolecules).

Subcommands:
  gff    convert a GFF3 file:   --to object     (contig -> scaffold)   or   --to component (scaffold -> contig)
  bed    convert a BED3-BED6 file in the same way (BED12 blocks are refused rather than silently mis-converted)
  chain  write a UCSC chain file (component -> object), one gap-free chain per W line
  gaps   write the AGP gaps as BED (0-based, half-open)
  check  verify the AGP against the component and object FASTAs, and test coordinate transforms in both directions
         on random intervals by sequence comparison (and against a chain file with --chain)

Conventions: AGP and GFF3 are 1-based, closed; BED and chain are 0-based, half-open. A component placed '-' is reverse
complemented in the object: start/end are mirrored and the strand flips. GFF3 phase is unchanged (it is defined
relative to the start of each CDS segment in the direction of transcription, and segments move as a unit). AGP
orientations '?', '0' and 'na' are treated as '+' (AGP 2.1: assembled forward).

What cannot be converted is never forced. A feature is written to the --unmapped table when it:
  - lies on a sequence absent from the AGP (NOT_IN_AGP);
  - overlaps a gap (GAP_SPANNING / IN_GAP);
  - crosses from one component into another (MULTI_COMPONENT);
  - covers only part of a component range that the AGP uses partially (OUTSIDE_COMPONENT_RANGE).
GFF3 features are kept or dropped as whole hierarchies (gene with all descendants), so no partial gene model is
written; --split-groups keeps the mappable members instead. IDs and all attributes are copied unchanged.
"""
import argparse
import bisect
import gzip
import hashlib
import json
import random
import sys
from collections import defaultdict

FLIP = {"+": "-", "-": "+"}
RC = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def opener(path, mode="rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


class AGP:
    """Parsed AGP. W rows: (object, obj_beg, obj_end, component, c_beg, c_end, ori); gaps: (object, beg, end)."""

    def __init__(self, path):
        self.w, self.gaps, self.obj_len, self.order = [], [], {}, []
        with opener(path) as fh:
            for n, line in enumerate(fh, 1):
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    raise SystemExit(f"{path}:{n}: AGP line with {len(f)} columns")
                obj, ob, oe, typ = f[0], int(f[1]), int(f[2]), f[4]
                if obj not in self.obj_len:
                    self.order.append(obj)
                self.obj_len[obj] = max(self.obj_len.get(obj, 0), oe)
                if typ in ("N", "U"):
                    self.gaps.append((obj, ob, oe))
                    continue
                ori = f[8] if f[8] in ("+", "-") else "+"
                cb, ce = int(f[6]), int(f[7])
                if ce - cb != oe - ob:
                    raise SystemExit(f"{path}:{n}: component span {ce - cb + 1} != object span {oe - ob + 1}")
                self.w.append((obj, ob, oe, f[5], cb, ce, ori))
        self.by_comp = defaultdict(list)
        self.by_obj = defaultdict(list)
        for r in self.w:
            self.by_comp[r[3]].append(r)
            self.by_obj[r[0]].append(r)
        for v in self.by_obj.values():
            v.sort(key=lambda r: r[1])
        self.obj_starts = {o: [r[1] for r in v] for o, v in self.by_obj.items()}
        self.gaps_by_obj = defaultdict(list)
        for g in self.gaps:
            self.gaps_by_obj[g[0]].append(g)
        self.comp_len = {c: max(r[5] for r in v) for c, v in self.by_comp.items()}

    # ---- transforms (1-based closed intervals) ----
    def to_object(self, comp, s, e, strand="."):
        rows = self.by_comp.get(comp)
        if not rows:
            return None, "NOT_IN_AGP"
        for obj, ob, oe, _, cb, ce, ori in rows:
            if cb <= s and e <= ce:
                if ori == "+":
                    return (obj, ob + s - cb, ob + e - cb, strand), None
                return (obj, ob + ce - e, ob + ce - s, FLIP.get(strand, strand)), None
        return None, "OUTSIDE_COMPONENT_RANGE"

    def to_component(self, obj, s, e, strand="."):
        rows = self.by_obj.get(obj)
        if not rows:
            return None, "NOT_IN_AGP"
        i = bisect.bisect_right(self.obj_starts[obj], s) - 1
        if i >= 0:
            o, ob, oe, comp, cb, ce, ori = rows[i]
            if s <= oe and e <= oe:
                if ori == "+":
                    return (comp, cb + s - ob, cb + e - ob, strand), None
                return (comp, ce - (e - ob), ce - (s - ob), FLIP.get(strand, strand)), None
        in_gap = any(gb <= s and e <= ge for _, gb, ge in self.gaps_by_obj.get(obj, ()))
        if in_gap:
            return None, "IN_GAP"
        if any(gb <= e and s <= ge for _, gb, ge in self.gaps_by_obj.get(obj, ())):
            return None, "GAP_SPANNING"
        return None, "MULTI_COMPONENT"

    def convert(self, to, seq, s, e, strand):
        return self.to_object(seq, s, e, strand) if to == "object" else self.to_component(seq, s, e, strand)

    def target_lengths(self, to):
        return self.obj_len if to == "object" else self.comp_len


# ---------------------------------------------------------------- GFF3
def parse_attrs(col9):
    out = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def cmd_gff(a):
    agp = AGP(a.agp)
    header, feats = [], []
    with opener(a.input) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                print("warning: ##FASTA section not converted (sequences are not carried)", file=sys.stderr)
                break
            if line.startswith("#"):
                if not line.startswith("##sequence-region"):
                    header.append(line)
                continue
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                raise SystemExit(f"GFF3 line with {len(f)} columns: {line[:80]!r}")
            feats.append(f)
    ids = {}
    for i, f in enumerate(feats):
        fid = parse_attrs(f[8]).get("ID")
        if fid:
            ids.setdefault(fid, []).append(i)   # multi-line features (e.g. CDS) share one ID

    def root(i, seen=()):
        par = parse_attrs(feats[i][8]).get("Parent")
        if not par:
            return i
        p = par.split(",")[0]
        if p not in ids or i in seen:
            return i
        return root(ids[p][0], seen + (i,))

    groups = defaultdict(list)
    for i in range(len(feats)):
        groups[root(i)].append(i)
    out_rows, unmapped = [], []
    for g, members in groups.items():
        res = {}
        for i in members:
            f = feats[i]
            res[i] = agp.convert(a.to, f[0], int(f[3]), int(f[4]), f[6])
        failed = [i for i in members if res[i][0] is None]
        targets = {res[i][0][0] for i in members if res[i][0]}
        if not a.split_groups and (failed or len(targets) > 1):
            why = res[failed[0]][1] if failed else "GROUP_SPLIT_ACROSS_SEQUENCES"
            for i in members:
                unmapped.append((feats[i], why if i in failed or not failed else "GROUP_MEMBER_UNMAPPED"))
            continue
        for i in members:
            f = feats[i]
            r, why = res[i]
            if r is None:
                unmapped.append((f, why))
                continue
            depth = 0
            j = i
            while parse_attrs(feats[j][8]).get("Parent") and depth < 20:
                p = parse_attrs(feats[j][8])["Parent"].split(",")[0]
                if p not in ids:
                    break
                j, depth = ids[p][0], depth + 1
            out_rows.append(([r[0], f[1], f[2], str(r[1]), str(r[2]), f[5], r[3], f[7], f[8]], g, depth))
    lens = agp.target_lengths(a.to)
    order = {o: k for k, o in enumerate(agp.order if a.to == "object" else sorted(lens))}
    gstart = {}
    for row, g, _ in out_rows:
        key = (order.get(row[0], 1 << 30), int(row[3]))
        gstart[g] = min(gstart.get(g, key), key)
    out_rows.sort(key=lambda t: (gstart[t[1]], t[1], t[2], int(t[0][3]), int(t[0][4])))
    used = []
    for row, _, _ in out_rows:
        if row[0] not in used:
            used.append(row[0])
    with opener(a.output, "wt") as out:
        if not header or not header[0].startswith("##gff-version"):
            out.write("##gff-version 3\n")
        for h in header:
            out.write(h)
        out.write(f"# coordinates converted by g2f_agp_lift.py --to {a.to} with AGP {a.agp}\n")
        for s in sorted(used, key=lambda x: order.get(x, 1 << 30)):
            out.write(f"##sequence-region {s} 1 {lens[s]}\n")
        for row, _, _ in out_rows:
            out.write("\t".join(row) + "\n")
    write_unmapped(a.unmapped, unmapped, gff=True)
    print(json.dumps({"input_features": len(feats), "written": len(out_rows), "unmapped": len(unmapped),
                      "unmapped_reasons": dict(sorted(_count(w for _, w in unmapped).items()))}), file=sys.stderr)
    return 0


def _count(it):
    c = defaultdict(int)
    for x in it:
        c[x] += 1
    return c


def write_unmapped(path, rows, gff):
    if not path:
        return
    with opener(path, "wt") as fh:
        fh.write("reason\t" + ("seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes" if gff
                               else "chrom\tstart\tend\trest") + "\n")
        for f, why in rows:
            fh.write(why + "\t" + "\t".join(f) + "\n")


# ---------------------------------------------------------------- BED
def cmd_bed(a):
    agp = AGP(a.agp)
    out_rows, unmapped = [], []
    n = 0
    with opener(a.input) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 6:
                raise SystemExit("BED with more than 6 columns (thick/blocks) is not supported; convert BED6 or GFF3")
            n += 1
            strand = f[5] if len(f) >= 6 else "."
            r, why = agp.convert(a.to, f[0], int(f[1]) + 1, int(f[2]), strand)
            if r is None:
                unmapped.append(([f[0], f[1], f[2], "\t".join(f[3:])], why))
                continue
            g = [r[0], str(r[1] - 1), str(r[2])] + f[3:5]
            if len(f) >= 6:
                g.append(r[3])
            out_rows.append(g)
    with opener(a.output, "wt") as out:
        for g in out_rows:
            out.write("\t".join(g) + "\n")
    write_unmapped(a.unmapped, unmapped, gff=False)
    print(json.dumps({"input_records": n, "written": len(out_rows), "unmapped": len(unmapped)}), file=sys.stderr)
    return 0


# ---------------------------------------------------------------- chain / gaps
def chain_lines(agp):
    """One gap-free chain per W line; target = component, query = object (UCSC convention: '-' query coordinates
    are on the reverse strand)."""
    out = []
    for k, (obj, ob, oe, comp, cb, ce, ori) in enumerate(agp.w, 1):
        size = ce - cb + 1
        olen = agp.obj_len[obj]
        if ori == "+":
            qs, qe = ob - 1, oe
        else:
            qs, qe = olen - oe, olen - ob + 1
        out.append(f"chain {size} {comp} {agp.comp_len[comp]} + {cb - 1} {ce} {obj} {olen} {ori} {qs} {qe} {k}\n"
                   f"{size}\n\n")
    return out


def cmd_chain(a):
    agp = AGP(a.agp)
    with opener(a.output, "wt") as out:
        out.writelines(chain_lines(agp))
    return 0


def cmd_gaps(a):
    agp = AGP(a.agp)
    with opener(a.output, "wt") as out:
        for obj, gb, ge in agp.gaps:
            out.write(f"{obj}\t{gb - 1}\t{ge}\n")
    return 0


class Chain:
    """Independent transform from a chain file (component -> object), for cross-checking the AGP transform."""

    def __init__(self, path):
        self.c = defaultdict(list)
        with opener(path) as fh:
            for line in fh:
                if line.startswith("chain"):
                    f = line.split()
                    t, ts, te = f[2], int(f[5]), int(f[6])
                    q, qsize, qstr, qs, qe = f[7], int(f[8]), f[9], int(f[10]), int(f[11])
                    self.c[t].append((ts, te, q, qsize, qstr, qs, qe))

    def to_object(self, comp, s, e):   # 1-based closed in, 1-based closed out
        for ts, te, q, qsize, qstr, qs, qe in self.c.get(comp, ()):
            if ts < s and e <= te:
                if te - ts != qe - qs:
                    return None
                a0, b0 = qs + (s - 1 - ts), qs + (e - ts)   # 0-based half-open on the query strand
                if qstr == "+":
                    return q, a0 + 1, b0
                return q, qsize - b0 + 1, qsize - a0
        return None


# ---------------------------------------------------------------- check
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                if name in seqs:
                    raise SystemExit(f"{path}: duplicate sequence ID {name}")
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def cmd_check(a):
    agp = AGP(a.agp)
    comp = read_fasta(a.component_fasta)
    obj = read_fasta(a.object_fasta)
    fails = []
    # 1. reconstruction and accounting
    used = defaultdict(int)
    for o in agp.order:
        parts = sorted([(r[1], "W", r) for r in agp.by_obj[o]] + [(g[1], "G", g) for g in agp.gaps_by_obj[o]])
        seq, pos = [], 1
        for beg, t, r in parts:
            if beg != pos:
                fails.append(f"{o}: AGP not contiguous at {pos}")
                break
            if t == "W":
                c = comp.get(r[3])
                if c is None:
                    fails.append(f"{o}: component {r[3]} not in component FASTA")
                    break
                s = c[r[4] - 1:r[5]]
                seq.append(s.translate(RC)[::-1] if r[6] == "-" else s)
                used[r[3]] += r[5] - r[4] + 1
                pos = r[2] + 1
            else:
                seq.append("N" * (r[2] - r[1] + 1))
                pos = r[2] + 1
        rebuilt = "".join(seq)
        if o not in obj:
            fails.append(f"{o}: object not in object FASTA")
        elif hashlib.sha256(rebuilt.encode()).digest() != hashlib.sha256(obj[o].encode()).digest():
            fails.append(f"{o}: reconstruction differs from object FASTA")
    extra_obj = set(obj) - set(agp.order)
    if extra_obj:
        fails.append(f"{len(extra_obj)} object FASTA records not in the AGP")
    not_whole = [c for c in comp if used.get(c, 0) != len(comp[c])]
    if not_whole and not a.allow_partial:
        fails.append(f"{len(not_whole)} components not used exactly once and whole (e.g. {not_whole[:3]})")
    # 2. random round trips, sequence-checked
    rng = random.Random(a.seed)
    chain = Chain(a.chain) if a.chain else None
    stats = defaultdict(int)
    rows = agp.w
    for _ in range(a.n if rows else 0):
        o, ob, oe, c, cb, ce, ori = rows[rng.randrange(len(rows))]
        span = ce - cb + 1
        ln = min(span, rng.choice([1, 50, 500, 5000, 20000]))
        mode = rng.random()
        s = cb if mode < 0.15 else (ce - ln + 1 if mode < 0.3 else rng.randint(cb, ce - ln + 1))
        e = s + ln - 1
        r, why = agp.to_object(c, s, e, "+")
        if r is None:
            fails.append(f"to_object failed inside a W range: {c}:{s}-{e} ({why})")
            continue
        cs = comp[c][s - 1:e]
        os_ = obj[r[0]][r[1] - 1:r[2]]
        if (cs.translate(RC)[::-1] if r[3] == "-" else cs) != os_:
            fails.append(f"sequence mismatch {c}:{s}-{e} -> {r}")
        back, why = agp.to_component(r[0], r[1], r[2], r[3])
        if back != (c, s, e, "+"):
            fails.append(f"round trip {c}:{s}-{e} -> {r} -> {back}")
        if chain:
            cr = chain.to_object(c, s, e)
            if cr != (r[0], r[1], r[2]):
                fails.append(f"chain transform {c}:{s}-{e} -> {cr} != AGP {r[:3]}")
        stats["intervals_" + ("minus" if ori == "-" else "plus")] += 1
        stats["at_component_end"] += s == cb or e == ce
    # 3. intervals touching gaps must be refused when converting back to components
    for o, gb, ge in agp.gaps[:a.n]:
        r, why = agp.to_component(o, max(1, gb - 10), min(agp.obj_len[o], ge + 10))
        if r is not None or why != "GAP_SPANNING":
            fails.append(f"gap-spanning interval {o}:{gb - 10}-{ge + 10} not refused ({r}, {why})")
        stats["gap_spanning_refused"] += r is None
    res = {"agp": a.agp, "objects": len(agp.order), "components": len(comp), "w_lines": len(agp.w),
           "gaps": len(agp.gaps), "status": "PASS" if not fails else "FAIL", "failures": fails[:20],
           "n_failures": len(fails), "tests": dict(stats), "chain_checked": bool(chain)}
    out = json.dumps(res, indent=1)
    if a.output:
        with open(a.output, "w") as fh:
            fh.write(out + "\n")
    print(out)
    return 0 if not fails else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("gff", "bed"):
        p = sub.add_parser(name)
        p.add_argument("--agp", required=True)
        p.add_argument("--to", required=True, choices=("object", "component"))
        p.add_argument("-i", "--input", required=True)
        p.add_argument("-o", "--output", required=True)
        p.add_argument("--unmapped", help="TSV of features that could not be converted, with the reason")
        if name == "gff":
            p.add_argument("--split-groups", action="store_true",
                           help="keep mappable members of a hierarchy whose other members cannot be converted")
    for name in ("chain", "gaps"):
        p = sub.add_parser(name)
        p.add_argument("--agp", required=True)
        p.add_argument("-o", "--output", required=True)
    p = sub.add_parser("check")
    p.add_argument("--agp", required=True)
    p.add_argument("--component-fasta", required=True)
    p.add_argument("--object-fasta", required=True)
    p.add_argument("--chain")
    p.add_argument("--n", type=int, default=2000, help="random intervals tested")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--allow-partial", action="store_true", help="components may be used partially / not at all")
    p.add_argument("--output")
    a = ap.parse_args(argv)
    return {"gff": cmd_gff, "bed": cmd_bed, "chain": cmd_chain, "gaps": cmd_gaps, "check": cmd_check}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
