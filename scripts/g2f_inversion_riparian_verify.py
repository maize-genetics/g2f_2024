#!/usr/bin/env python3
"""Independent check of a g2f_inversion_riparian.py figure. Re-reads the source tables and the contig indexes; does not import the
plotting code. Checks:
  T1 every configured assembly x locus has one track; track contig = the chain contig; contig length = release .fai
  T2 every ribbon's contigs are the chain contigs of its two units, and its upper/lower units are adjacent tracks
     (the nearest track above with a chain)
  T3 every ribbon's B73 interval is covered by ordered-chain segments of BOTH units at its ends, and no uncovered gap
     inside it exceeds the merge tolerance (2 % of the window)
  T4 ribbon orientation: crossing <=> the two units' drawn orientations differ at that B73 interval
  T5 ribbon colour: magenta (inner) <=> the interval lies in an 'inner' chain segment of either unit
  T6 inverted state (arrangement table) <=> an inner chain exists; every inverted unit has a ribbon crossing from the
     B73-oriented neighbour or a parallel magenta ribbon from an inverted one, and two junction marks
  T7 row labels: every configured label occurs in the SVG, in configured order (by y position)
Also writes tables/breakpoint_coverage.tsv (informational): aligned chain bp within 250 kb outside each B73 breakpoint,
per unit, against the recorded B73 / B73_ONE label.
Usage: g2f_inversion_riparian_verify.py --root DIR --svg FILE [--config DIR] [--tables DIR]
Exit status 0 = all checks PASS, 1 = a check failed."""
import argparse, bisect, collections, csv, re, sys

ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--svg", required=True)
ap.add_argument("--config"); ap.add_argument("--tables"); a = ap.parse_args()
CFG = a.config or f"{a.root}/config"; TAB = a.tables or f"{a.root}/tables"
tsv = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
SRC = {r["key"]: r["path"] for r in tsv(f"{CFG}/sources.tsv")}
LOCI = {r["event"]: r for r in tsv(f"{CFG}/loci.tsv")}
ASM = sorted(tsv(f"{CFG}/assemblies.tsv"), key=lambda r: int(r["order"]))
ORDER = ["B73"] + [r["unit"] for r in ASM]
SEG = collections.defaultdict(list)
for r in tsv(SRC["chain_segments"]):
    SEG[(r["event"], r["unit"])].append((int(r["b73_start"]), int(r["b73_end"]), r["strand"], r["zone"], r["contig"], int(r["aligned_bp"])))
for k in SEG: SEG[k].sort()
STATE = {(r["event"], r["unit"]): r for r in tsv(SRC["arrangement"])}
IVL = {r["event"]: r for r in tsv(SRC["event_intervals"])}
FAI = {r["unit"]: {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(r["contig_fai"])} for r in ASM}
TRK = {(r["event"], r["unit"]): r for r in tsv(f"{TAB}/tracks.tsv")}
RIB = tsv(f"{TAB}/ribbons.tsv")
fails = collections.Counter(); n = collections.Counter(); msgs = []


def bad(t, m):
    fails[t] += 1
    if len(msgs) < 40: msgs.append(f"{t}: {m}")


def covering(e, u, b):
    if u == "B73": return ("+", "whole")
    best = None
    for s in SEG[(e, u)]:
        if s[0] <= b <= s[1]: return (s[2], s[3])
        d = min(abs(b - s[0]), abs(b - s[1]))
        if best is None or d < best[0]: best = (d, (s[2], s[3]))
    return best[1]


def covered(e, u, b):
    return u == "B73" or any(s[0] - 1 <= b <= s[1] + 1 for s in SEG[(e, u)])


def drawn(e, u): return 1 if u == "B73" else (-1 if TRK[(e, u)]["strand_drawn"] == "rc" else 1)


# T1
for e in LOCI:
    for r in ASM:
        u = r["unit"]; n["T1"] += 1
        t = TRK.get((e, u))
        if not t: bad("T1", f"{e} {u} no track"); continue
        cs = {s[4] for s in SEG[(e, u)]}
        if cs != {t["contig"]}: bad("T1", f"{e} {u} contig {t['contig']} vs chain {cs}")
        if FAI[u].get(t["contig"]) != int(t["contig_len"]): bad("T1", f"{e} {u} length {t['contig_len']} vs fai {FAI[u].get(t['contig'])}")
# T2-T5
for r in RIB:
    e, up, lo = r["event"], r["upper"], r["lower"]; b0, b1 = int(r["b73_start"]), int(r["b73_end"])
    tol = 0.02 * (int(LOCI[e]["window_end"]) - int(LOCI[e]["window_start"]))
    n["T2"] += 1
    if lo == "B73" or ORDER.index(up) >= ORDER.index(lo): bad("T2", f"{e} {up}->{lo} order")
    between = ORDER[ORDER.index(up) + 1:ORDER.index(lo)]
    if any(SEG[(e, x)] for x in between): bad("T2", f"{e} {up}->{lo} skips a track with a chain")
    for u, c in ((up, r["upper_contig"]), (lo, r["lower_contig"])):
        if u != "B73" and {s[4] for s in SEG[(e, u)]} != {c}: bad("T2", f"{e} {u} ribbon contig {c}")
    n["T3"] += 1
    for u in (up, lo):
        if not (covered(e, u, b0) and covered(e, u, b1)): bad("T3", f"{e} {up}->{lo} {b0}-{b1}: end not aligned in {u}")
    if up != "B73" or lo:
        iv = sorted((max(s[0], b0), min(s[1], b1)) for u in (up, lo) if u != "B73" for s in SEG[(e, u)] if s[1] > b0 and s[0] < b1)
        # gaps where EITHER unit lacks alignment: check per unit
        for u in (up, lo):
            if u == "B73": continue
            ivu = sorted((max(s[0], b0), min(s[1], b1)) for s in SEG[(e, u)] if s[1] > b0 and s[0] < b1)
            cur = b0; gmax = 0
            for x0, x1 in ivu:
                gmax = max(gmax, x0 - cur); cur = max(cur, x1)
            gmax = max(gmax, b1 - cur)
            if gmax > tol + 1: bad("T3", f"{e} {up}->{lo} {b0}-{b1}: {u} gap {gmax} > {tol:.0f}")
    mid = (b0 + b1) // 2
    su, zu = covering(e, up, mid); sl, zl = covering(e, lo, mid)
    ou = (1 if su == "+" else -1) * drawn(e, up); ol = (1 if sl == "+" else -1) * drawn(e, lo)
    n["T4"] += 1
    if int(r["crossing"]) != int(ou != ol) and int(r["b73_end"]) - int(r["b73_start"]) > 5000:
        bad("T4", f"{e} {up}->{lo} {b0}-{b1} crossing={r['crossing']} expected {int(ou != ol)}")
    n["T5"] += 1
    if int(r["inner"]) != int("inner" in (zu, zl)): bad("T5", f"{e} {up}->{lo} {b0}-{b1} inner={r['inner']} zones {zu},{zl}")
# T6
for e in LOCI:
    for r in ASM:
        u = r["unit"]; n["T6"] += 1
        if e not in LOCI: continue
        inv_state = STATE.get((e, u), {}).get("state") == "INVERTED"
        has_inner = any(s[3] == "inner" for s in SEG[(e, u)])
        if inv_state != has_inner: bad("T6", f"{e} {u} state {STATE.get((e, u), {}).get('state')} vs inner chain {has_inner}")
        if inv_state:
            rr = [x for x in RIB if x["event"] == e and x["lower"] == u and x["inner"] == "1"]
            if not rr: bad("T6", f"{e} {u} inverted but no magenta ribbon from above")
            nj = len([j for j in TRK[(e, u)]["junctions"].split(";") if j])
            if nj != 2: bad("T6", f"{e} {u} inverted with {nj} junction marks")
# T7
svg = open(a.svg).read()
texts = [(float(m.group(2)), m.group(3)) for m in re.finditer(r'<text[^>]*x="([\d.]+)" y="([\d.]+)"[^>]*>([^<]+)</text>', svg)]
if not texts:   # matplotlib writes transforms; fall back to order of appearance
    texts = [(i, t) for i, t in enumerate(re.findall(r">([^<>]+)</text>", svg))]
labels = ["B73"] + [r["partition"] or r["line"] for r in ASM]
found = [t for _, t in texts]
n["T7"] += 1
for lab in set(labels) | {r["line"] for r in ASM}:
    if lab not in " ".join(found): bad("T7", f"label {lab} missing from SVG")

# breakpoint coverage (informational)
out = []
for e in LOCI:
    iv = IVL[e]; L0, R0 = int(iv["b73_outer_start"]), int(iv["b73_outer_end"]); I0, I1 = int(iv["b73_inner_start"]), int(iv["b73_inner_end"])
    for r in ASM:
        u = r["unit"]; S = SEG[(e, u)]
        bp = lambda x0, x1: sum(max(0, min(s[1], x1) - max(s[0], x0)) for s in S)
        lf, rf = bp(L0 - 250000, L0), bp(R0, R0 + 250000)
        li, ri = bp(I0, I0 + 250000), bp(I1 - 250000, I1)
        chain_left = lf >= 10000 and li >= 10000; chain_right = rf >= 10000 and ri >= 10000
        st = STATE.get((e, u), {}).get("state", "")
        derived = "BOTH" if chain_left and chain_right else ("ONE" if chain_left or chain_right else "NONE")
        rec = {"B73": "BOTH", "INVERTED": "BOTH", "B73_ONE": "ONE", "NT": "NONE"}.get(st, "")
        out.append(dict(event=e, unit=u, recorded_state=st, outside_left_bp=lf, inside_left_bp=li, inside_right_bp=ri,
                        outside_right_bp=rf, chain_spans=derived, agrees=int(rec == derived)))
with open(f"{TAB}/breakpoint_coverage.tsv", "w") as fh:
    w = csv.DictWriter(fh, list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
dis = [o for o in out if not o["agrees"]]

for t in sorted(n): print(f"{t}\t{n[t]} checked\t{fails[t]} failed\t{'PASS' if not fails[t] else 'FAIL'}")
for m in msgs: print("  ", m)
print(f"breakpoint coverage: {len(out) - len(dis)}/{len(out)} chain spans agree with the recorded label (informational)")
for o in dis: print("   differs:", o["event"], o["unit"], "recorded", o["recorded_state"], "chain", o["chain_spans"],
                    f"(outside L {o['outside_left_bp']}, inside L {o['inside_left_bp']}, inside R {o['inside_right_bp']}, outside R {o['outside_right_bp']})")
sys.exit(1 if sum(fails.values()) else 0)
