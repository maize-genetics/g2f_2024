#!/usr/bin/env python3
"""Independent verification of an AGP-level scaffold revision (output of g2f_agp_revise.py apply). Imports no builder code;
it re-reads the old and new AGPs, the plan, the contigs and both scaffold FASTAs.

  g2f_agp_revise_verify.py --old-agp OLD.agp --new-agp NEW.agp --plan plan.tsv [--unit U] --contigs C.fa.gz
                           --old-scaffolds OLD.fa.gz --new-scaffolds NEW.fa.gz --out verify.json
                           [--chain NEW.chain] [--gaps-bed NEW.gaps.bed] [--map NEW.contig_to_scaffold.tsv]
                           [--contigs-sha256 HEX] [--intervals 300] [--seed 1]
Checks (each PASS/FAIL, exit 1 on any failure):
  agp_syntax       per object: coordinates contiguous from 1, part numbers 1..n, W length = component span, gaps N/U with
                   length = span, orientation +/-, no object starts or ends with a gap, object length = new FASTA .fai
  agp_reconstructs g2f_fasta_gz_verify.py agp: each new scaffold record equals its AGP reconstruction, every contig used once
                   and whole (sibling script, independent code)
  diff_equals_plan the old->new AGP difference is exactly the plan: flips change only the orientation field; unplaced
                   components (and the gap before each) leave their object and reappear as 1..L + objects; nothing else
  sequence         every record not named in the plan has the same sequence MD5 as in the old scaffold FASTA; flipped spans
                   are the reverse complement of the old span; a truncated object is an identical prefix of the old one
  accounting       component set = contig set, W bp = contig bp, scaffold bp = W bp + gap bp, old - new bp = removed gap bp
  derived_tables   (when given) chain blocks, gap BED and contig_to_scaffold rows agree with the new AGP and .fai
  lift_roundtrip   g2f_agp_lift.py bed (sibling script) maps random intervals contig -> scaffold; sequences agree
                   (reverse complement on - placements); a third of the intervals fall in changed components
  contigs_sha256   (when given) the contig FASTA has the expected sha256
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def agp(p):
    return [l.rstrip("\n").split("\t") for l in open(p) if not l.startswith("#") and l.strip()]


def md5s(path):
    out, name, h = {}, None, None
    with (gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")) as fh:
        for l in fh:
            if l[:1] == b">":
                if name:
                    out[name] = h.hexdigest()
                name, h = l[1:].split()[0].decode(), hashlib.md5()
            else:
                h.update(l.strip().upper())
    if name:
        out[name] = h.hexdigest()
    return out


def fx(fa, region):
    s = subprocess.run(["samtools", "faidx", fa, region], capture_output=True, text=True, check=True).stdout
    return "".join(s.splitlines()[1:]).upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--old-agp", "--new-agp", "--plan", "--contigs", "--old-scaffolds", "--new-scaffolds", "--out"):
        ap.add_argument(k, required=True)
    for k in ("--unit", "--chain", "--gaps-bed", "--map", "--contigs-sha256"):
        ap.add_argument(k)
    ap.add_argument("--intervals", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"refusing to overwrite {a.out}")
    for fa in (a.contigs, a.old_scaffolds, a.new_scaffolds):
        for ix in [fa + ".fai"] + ([fa + ".gzi"] if fa.endswith(".gz") else []):
            if not os.path.exists(ix):
                sys.exit(f"missing index {ix} (samtools faidx would write it next to the input; create it yourself)")
    res, fails = {}, []

    def chk(name, ok, detail):
        res[name] = {"pass": bool(ok), "detail": detail}
        if not ok:
            fails.append(name)

    plan = list(csv.DictReader(open(a.plan), delimiter="\t"))
    if a.unit is not None:
        plan = [p for p in plan if p.get("unit") == a.unit]
    flips = {(p["object"], p["component"]): p for p in plan if p["action"] == "REVERSE_ORIENTATION"}
    unpl = {p["component"]: p for p in plan if p["action"] == "UNPLACE_TERMINAL_EXTENSION"}
    old, new = agp(a.old_agp), agp(a.new_agp)
    fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.new_scaffolds + ".fai")}
    cfai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.contigs + ".fai")}
    if a.contigs_sha256:
        chk("contigs_sha256", sha(a.contigs) == a.contigs_sha256, "contig FASTA sha256 vs expected")
    # agp_syntax
    by = collections.OrderedDict()
    for x in new:
        by.setdefault(x[0], []).append(x)
    bad = []
    for o, xs in by.items():
        pos = 1
        for i, x in enumerate(xs, 1):
            s, e = int(x[1]), int(x[2])
            if s != pos or int(x[3]) != i or e < s:
                bad.append((o, i, "coordinates/part number"))
            if x[4] in ("N", "U"):
                if int(x[5]) != e - s + 1:
                    bad.append((o, i, "gap length"))
            elif x[8] not in "+-" or int(x[7]) - int(x[6]) != e - s:
                bad.append((o, i, "component span/orientation"))
            pos = e + 1
        if xs[0][4] in ("N", "U") or xs[-1][4] in ("N", "U"):
            bad.append((o, "object starts or ends with a gap"))
        if fai.get(o) != pos - 1:
            bad.append((o, "length differs from .fai"))
    chk("agp_syntax", not bad and list(by) == list(fai), f"{len(by)} objects; problems {bad[:5]}")
    # agp_reconstructs
    j = a.out + ".fasta_gz_verify.json"
    p = subprocess.run([sys.executable, "-I", os.path.join(HERE, "g2f_fasta_gz_verify.py"), "agp", "--gz", a.new_scaffolds,
                        "--agp", a.new_agp, "--components", a.contigs, "--output", j], capture_output=True, text=True)
    chk("agp_reconstructs", p.returncode == 0, f"g2f_fasta_gz_verify.py agp exit {p.returncode}; {j}")
    # diff_equals_plan
    key = lambda x: tuple(x[:5])
    ok_, nk = {key(x): x for x in old}, {key(x): x for x in new}
    oldpart = {(x[0], int(x[3])): x for x in old}
    unexpected = []
    for k, x in ok_.items():
        y = nk.get(k)
        if y == x:
            continue
        if y and x[4] not in ("N", "U") and (x[0], x[5]) in flips and y[:8] == x[:8] and y[8] == flips[(x[0], x[5])]["new_orientation"] != x[8]:
            continue
        if y is None and x[4] not in ("N", "U") and x[5] in unpl and x[0] == unpl[x[5]]["object"]:
            continue
        if y is None and x[4] in ("N", "U"):
            nxt = oldpart.get((x[0], int(x[3]) + 1))
            if nxt and nxt[5] in unpl:
                continue
        unexpected.append(("old", x))
    for k, y in nk.items():
        if k in ok_:
            continue
        if y[0] in unpl and y[1:] == ["1", y[7], "1", y[4], y[0], "1", y[7], "+"]:
            continue
        unexpected.append(("new", y))
    nflip = sum(1 for k, x in ok_.items() if x[4] not in ("N", "U") and (x[0], x[5]) in flips and nk.get(k, x)[8] != x[8])
    chk("diff_equals_plan", not unexpected and nflip == len(flips) and all(c in by for c in unpl),
        f"flips {nflip}/{len(flips)}; unplaced {len(unpl)}; unexpected lines {unexpected[:4]}")
    # sequence
    o5, n5 = md5s(a.old_scaffolds), md5s(a.new_scaffolds)
    named = {p_["object"] for p_ in plan} | set(unpl)
    diff = [o for o in n5 if o not in named and n5[o] != o5.get(o)]
    det, ok5 = [], not diff
    for (o, c), p_ in flips.items():
        x = next(z for z in old if z[0] == o and z[4] not in ("N", "U") and z[5] == c)
        s, e = int(x[1]), int(x[2])
        good = fx(a.new_scaffolds, f"{o}:{s}-{e}") == fx(a.old_scaffolds, f"{o}:{s}-{e}").translate(RC)[::-1]
        rest = [(1, s - 1), (e + 1, fai[o])]
        same_rest = all(fx(a.new_scaffolds, f"{o}:{u}-{v}") == fx(a.old_scaffolds, f"{o}:{u}-{v}") for u, v in rest if v >= u)
        det.append(f"{o}:{s}-{e} reverse complement {good}; rest of object identical {same_rest}")
        ok5 &= good and same_rest
    for o in {p_["object"] for p_ in unpl.values()}:
        L = fai[o]
        same = hashlib.md5(fx(a.new_scaffolds, f"{o}:1-{L}").encode()).hexdigest() == hashlib.md5(fx(a.old_scaffolds, f"{o}:1-{L}").encode()).hexdigest()
        det.append(f"{o}: identical prefix of the old object {same}")
        ok5 &= same
    for c in unpl:
        same = n5.get(c) == hashlib.md5(fx(a.contigs, c).encode()).hexdigest()
        det.append(f"{c}: new record equals the contig {same}")
        ok5 &= same
    chk("sequence", ok5, f"unchanged records differing {diff[:5]}; " + "; ".join(det))
    # accounting
    comps = [x[5] for x in new if x[4] not in ("N", "U")]
    wbp = sum(int(x[7]) - int(x[6]) + 1 for x in new if x[4] not in ("N", "U"))
    gbp = sum(int(x[5]) for x in new if x[4] in ("N", "U"))
    ogbp = sum(int(x[5]) for x in old if x[4] in ("N", "U"))
    oldbp = sum(int(l.split("\t")[1]) for l in open(a.old_scaffolds + ".fai"))
    acc = (sorted(comps) == sorted(cfai) and len(set(comps)) == len(comps) and wbp == sum(cfai.values())
           and sum(fai.values()) == wbp + gbp and oldbp - sum(fai.values()) == ogbp - gbp)
    chk("accounting", acc, f"contigs {len(cfai)} used {len(comps)}; W bp {wbp}; scaffold bp {sum(fai.values())} "
        f"(old {oldbp}); gap bp removed {ogbp - gbp}")
    # derived tables
    w = {x[5]: x for x in new if x[4] not in ("N", "U")}
    probs = []
    if a.chain:
        cl = [l.split() for l in open(a.chain) if l.startswith("chain")]
        for c in cl:
            x = w[c[2]]
            S, E, SL = int(x[1]), int(x[2]), fai[x[0]]
            exp = (S - 1, E) if x[8] == "+" else (SL - E, SL - S + 1)
            if c[7] != x[0] or c[9] != x[8] or (int(c[10]), int(c[11])) != exp or int(c[3]) != cfai[c[2]] or int(c[8]) != SL:
                probs.append(("chain", c[2]))
        if len(cl) != len(w):
            probs.append(("chain count", len(cl)))
    if a.gaps_bed:
        g = [l.split() for l in open(a.gaps_bed)]
        if [(x[0], int(x[1]) + 1, int(x[2])) for x in g] != [(x[0], int(x[1]), int(x[2])) for x in new if x[4] in ("N", "U")]:
            probs.append(("gaps.bed",))
    if a.map:
        for r in csv.DictReader(open(a.map), delimiter="\t"):
            x = w[r["contig"]]
            if (r["scaffold"], int(r["scaffold_start"]), int(r["scaffold_end"]), r["orientation"]) != (x[0], int(x[1]), int(x[2]), x[8]):
                probs.append(("map", r["contig"]))
    if a.chain or a.gaps_bed or a.map:
        chk("derived_tables", not probs, f"problems {probs[:5]}")
    # lift round trip
    random.seed(a.seed)
    focus = [p_["component"] for p_ in plan]
    bed = []
    for i in range(a.intervals):
        c = random.choice(focus) if (focus and i < a.intervals // 3) else random.choice(sorted(cfai))
        L = cfai[c]
        n = min(L, random.choice([50, 500, 5000]))
        s = random.choice([0, L - n, random.randint(0, L - n)])
        bed.append((c, s, s + n, f"f{i}", 0, random.choice("+-")))
    bi, bo, bu = a.out + ".lift_in.bed", a.out + ".lift_out.bed", a.out + ".lift_unmapped.tsv"
    with open(bi, "w") as fh:
        fh.writelines("\t".join(map(str, b)) + "\n" for b in bed)
    subprocess.run([sys.executable, "-I", os.path.join(HERE, "g2f_agp_lift.py"), "bed", "--agp", a.new_agp, "--to", "object",
                    "-i", bi, "-o", bo, "--unmapped", bu], check=True, capture_output=True)
    lifted = {l.split("\t")[3]: l.rstrip("\n").split("\t") for l in open(bo) if not l.startswith(("#", "track"))}
    bad8 = []
    for c, s, e, nm, _, st in bed:
        y = lifted.get(nm)
        if not y:
            bad8.append((nm, "unmapped"))
            continue
        q, t = fx(a.contigs, f"{c}:{s + 1}-{e}"), fx(a.new_scaffolds, f"{y[0]}:{int(y[1]) + 1}-{y[2]}")
        o = w[c][8]
        if (t != q if o == "+" else t != q.translate(RC)[::-1]) or ((y[5] != st) != (o == "-")):
            bad8.append((nm, c, y[0]))
    chk("lift_roundtrip", not bad8, f"{len(bed)} intervals; failures {bad8[:3]}")
    out = {"schema": "g2f-agp-revise-verify/1", "pass": not fails, "fails": fails, "checks": res,
           "inputs": {k: getattr(a, k) for k in ("old_agp", "new_agp", "plan", "unit", "contigs", "old_scaffolds", "new_scaffolds")},
           "new_agp_sha256": sha(a.new_agp), "new_scaffolds_sha256": sha(a.new_scaffolds)}
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(("PASS" if not fails else "FAIL " + ",".join(fails)) + f"  ({a.out})")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
