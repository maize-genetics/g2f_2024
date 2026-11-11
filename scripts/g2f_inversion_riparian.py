#!/usr/bin/env python3
"""Vertical, GENESPACE-style riparian figure of validated inversions relative to the B73 v5
reference (reference chromosomes chr1..chr10; the reference row is labelled "B73"), with every assembly drawn. Requires matplotlib.

Everything is driven by three configuration tables (no assembly or locus is named in this code):
  config/sources.tsv     input tables (alignment chains, breakpoint intervals, states, junction reads, B73 .fai)
  config/loci.tsv        windows to draw, in order; consecutive windows on one chromosome share a heading and a break
  config/assemblies.tsv  one row per track, in display order: unit id (as in the chain table), line, partition,
                         source, contig .fai (contig lengths, for contig ends)

Drawing rules (identical for every assembly):
  * Tracks are contigs, placed in B73 display coordinates by a slope +-1 map fitted to the flank chain (weighted median
    offset), so lengths are true. A contig that runs opposite to B73 is drawn reverse-complemented ("rc").
  * Track fill: grey = contig sequence; coloured = sequence in the ordered alignment chain (chromosome colour = B73
    orientation, magenta = the reversed interval). A dark bar marks a contig end inside the window.
  * Ribbons join adjacent tracks. B73 -> first track uses the chain segments directly. Between two assemblies a ribbon
    joins the parts of both contigs whose chain segments cover the SAME B73 interval (B73-anchored; there is no direct
    assembly-assembly alignment). Positions inside a segment are linearly interpolated. Consecutive pieces are merged only
    when the gap is <= MERGE_FRAC of the window in B73 and on both tracks, with unchanged orientation; larger gaps stay
    open. Nothing is drawn where either track lacks ordered alignment, so no correspondence crosses an unaligned gap.
  * Ribbon colour: magenta if the piece lies in the reversed (inner) zone of either track; chromosome colour otherwise.
    Crossing ribbons mark an orientation difference between the two tracks, so colour is not the only cue.
  * Junction triangles: filled = at least one single HiFi read crosses the junction; open = covered only by
    overlapping reads (no single crossing read). Counts are kept in tables/tracks.tsv, not printed.
  * Partitions of one line (hap1/hap2) are bracketed under the line name; they are not chromosome-scale phase.
Outputs: PREFIX.svg/.pdf/.png, tables/ribbons.tsv, tables/tracks.tsv, tables/guides.tsv and tables/tags.tsv when
         applicable (inputs for scripts/g2f_inversion_riparian_verify.py).
Usage: g2f_inversion_riparian.py --root DIR --out PREFIX [--config DIR] [--tables DIR] [--legend TSV] [--suggest-order]"""
import argparse, bisect, collections, csv, itertools, logging, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True, help="figure directory holding config/ and tables/")
ap.add_argument("--out", help="output prefix (svg, pdf, png)")
ap.add_argument("--config", help="config directory (default ROOT/config)")
ap.add_argument("--tables", help="output table directory (default ROOT/tables)")
ap.add_argument("--width", type=float, default=183.0, help="figure width in mm")
ap.add_argument("--legend", help="legend table (default CONFIG/legend.tsv)")
ap.add_argument("--junction-marks", choices=["show", "hide"], default="show",
                help="'hide' draws no junction triangles or numbers (records stay in tables/tracks.tsv)")
ap.add_argument("--headings", choices=["show", "hide"], default="show", help="'hide' omits chromosome headings")
ap.add_argument("--ideogram", choices=["genome", "focus"], default="genome",
                help="'focus' draws only the chromosomes of the displayed loci, enlarged, with connectors that end on "
                     "the B73 inversion interval edges; also drops the '//' window break")
ap.add_argument("--flag-tag", choices=["show", "hide"], default="show", help="'hide' omits the dagger flag tag")
ap.add_argument("--tag-anchor", choices=["window", "contig"], default="window",
                help="'contig' places the rc tag a fixed distance right of the drawn contig end, not the window edge")
ap.add_argument("--suggest-order", action="store_true", help="print a row order that minimises arrangement changes")
a = ap.parse_args()

import matplotlib
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Polygon
from matplotlib.path import Path

MM = 1 / 25.4
matplotlib.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 5, "svg.fonttype": "none",
                            "pdf.fonttype": 42, "svg.hashsalt": "g2f"})
INK, INK2, GRID, CONTIG = "#1a1a1a", "#5a5a5a", "#c8c8c8", "#e3e3e3"
COL = {"chr3": "#5b7fa6", "chr4": "#4f9a8f", "chr6": "#c08a3e"}       # per-chromosome colours; others DEFAULT_COL
DEFAULT_COL = "#7f7f7f"
INV = "#b0357e"
HAIR = 0.25
MERGE_FRAC = 0.02


def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))


cfg = a.config or os.path.join(a.root, "config")
LEGF = a.legend or f"{cfg}/legend.tsv"
TAB = a.tables or os.path.join(a.root, "tables")
SRC = {r["key"]: r["path"] for r in tsv(f"{cfg}/sources.tsv")}
LOCI = sorted(tsv(f"{cfg}/loci.tsv"), key=lambda r: int(r["order"]))
ASM = sorted(tsv(f"{cfg}/assemblies.tsv"), key=lambda r: int(r["order"]))
if len({r["unit"] for r in ASM}) != len(ASM): sys.exit("assemblies.tsv: duplicate unit")

# ---------- inputs ----------
SEG = collections.defaultdict(list)          # (event, unit) -> segments
for r in tsv(SRC["chain_segments"]):
    SEG[(r["event"], r["unit"])].append(dict(contig=r["contig"], zone=r["zone"], strand=r["strand"],
                                             b0=int(r["b73_start"]), b1=int(r["b73_end"]),
                                             c0=int(r["contig_start"]), c1=int(r["contig_end"]), bp=int(r["aligned_bp"])))
IVL = {r["event"]: r for r in tsv(SRC["event_intervals"])}
STATE = {(r["event"], r["unit"]): r for r in tsv(SRC["arrangement"])}
B73LEN = {}
for l in open(SRC["b73_fai"]):
    f = l.split("\t")
    if f[0] in [f"chr{i}" for i in range(1, 11)]: B73LEN[f[0]] = int(f[1])
CLEN = {}                                    # unit -> {contig: length}
for r in ASM:
    p = r.get("contig_fai", "")
    CLEN[r["unit"]] = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p)} if p and os.path.exists(p) else {}

# junctions: normalised to (event, unit) -> [dict(side, contig_lo, contig_hi, reads, source)]
JR = collections.defaultdict(list)
shared = collections.defaultdict(list)
if SRC.get("junctions_shared"):
    for r in tsv(SRC["junctions_shared"]):
        ev = r["junction"].split("_")[1]
        lo, hi = (int(x) for x in r["interval"].split("-"))
        shared[(ev, r["unit"])].append(dict(side="left" if r["junction"].endswith("_L") else "right", contig_lo=lo,
                                            contig_hi=hi, reads=int(r["full_span"]), source="junctions_shared"))
if SRC.get("junctions_chain"):
    for r in tsv(SRC["junctions_chain"]):
        if (r["event"], r["unit"]) in shared: continue          # the shared-sequence measurement takes precedence
        lo, hi = (int(x) for x in r["contig_junction"].split("-"))
        JR[(r["event"], r["unit"])].append(dict(side=r["b73_breakpoint"], contig_lo=lo, contig_hi=hi,
                                                reads=int(r["single_reads_crossing"] or 0), source="junctions_chain"))
for k, v in shared.items(): JR[k] = v
# audited junction-spanning read counts replace the sources above when configured
# (event, unit, side, contig lo/hi, status assessed|unassessed, reads, label, same_reads_as)
if SRC.get("junction_counts"):
    JR = collections.defaultdict(list)
    for r in tsv(SRC["junction_counts"]):
        if r["status"] != "assessed":
            JR[(r["event"], r["unit"])].append(dict(side=r["side"], contig_lo=None, contig_hi=None, reads=None, label="n.a.",
                                                    same="", source="junction_counts")); continue
        JR[(r["event"], r["unit"])].append(dict(side=r["side"], contig_lo=int(r["lo"]), contig_hi=int(r["hi"]), reads=int(r["reads"]),
                                                label=r["label"], same=r.get("same_reads_as", ""), source=r["definition"],
                                                hidden=r.get("display", "show") == "hide"))
# optional per-track interval features (e.g. a high-copy repeat block), contig coordinates
FEAT = collections.defaultdict(list)
if os.path.exists(f"{cfg}/track_features.tsv"):
    for r in tsv(f"{cfg}/track_features.tsv"):
        FEAT[(r["event"], r["unit"])].append((int(r["contig_start"]), int(r["contig_end"]), r["color"]))


def is_inv(e, u):
    s = STATE.get((e, u), {}).get("state", "")
    return s == "INVERTED"


# ---------- optional: suggest an order ----------
if a.suggest_order:
    lines = collections.OrderedDict()
    for r in ASM:
        lines.setdefault(r["line"], []).append(r["unit"])
    pat = {ln: tuple(int(any(is_inv(L["event"], u) for u in us)) for L in LOCI) for ln, us in lines.items()}
    groups = collections.OrderedDict()
    for ln, p in pat.items(): groups.setdefault(p, []).append(ln)
    keys = list(groups)
    ham = lambda x, y: sum(i != j for i, j in zip(x, y))
    zero = tuple(0 for _ in LOCI)
    if len(keys) <= 8:
        best = min(itertools.permutations(keys), key=lambda o: sum(ham(x, y) for x, y in zip((zero,) + o, o)))
    else:                                                          # greedy nearest neighbour
        best, cur, rest = [], zero, list(keys)
        while rest:
            nxt = min(rest, key=lambda k: ham(cur, k)); best.append(nxt); rest.remove(nxt); cur = nxt
    cost = sum(ham(x, y) for x, y in zip((zero,) + tuple(best), best))
    print(f"# pattern over {[L['event'] for L in LOCI]} (1 = inverted in any partition); changes along the order = {cost}")
    for p in best: print("".join(map(str, p)), ", ".join(groups[p]))
    sys.exit(0)
if not a.out: sys.exit("--out required")


# ---------- geometry helpers ----------
def cmap(e, u):
    """slope +-1 map contig -> B73 coordinate, fitted on the flank chain (weighted median offset)."""
    S = [x for x in SEG[(e, u)] if x["zone"] != "inner"] or SEG[(e, u)]
    plus = sum(x["bp"] for x in S if x["strand"] == "+"); minus = sum(x["bp"] for x in S if x["strand"] == "-")
    sgn = -1 if minus > plus else 1
    offs = sorted(((x["b0"] + x["b1"]) / 2 - sgn * (x["c0"] + x["c1"]) / 2, x["bp"]) for x in S)
    tot = sum(w for _, w in offs); acc = 0
    for o, w in offs:
        acc += w
        if acc >= tot / 2: return sgn, o


def interp(s, b):
    """contig coordinate aligned to B73 coordinate b inside segment s (linear)."""
    t = (b - s["b0"]) / max(s["b1"] - s["b0"], 1)
    return s["c0"] + t * (s["c1"] - s["c0"]) if s["strand"] == "+" else s["c1"] - t * (s["c1"] - s["c0"])


class Track:
    """one row in one window: B73 (segments=None) or a contig."""
    def __init__(self, e, unit):
        self.e, self.unit = e, unit
        if unit == "B73":
            self.segs, self.sgn, self.off, self.contig = None, 1, 0, None
        else:
            self.segs = sorted(SEG.get((e, unit), []), key=lambda s: s["b0"])
            self.contig = self.segs[0]["contig"] if self.segs else None
            if self.segs:
                if len({s["contig"] for s in self.segs}) != 1: sys.exit(f"{e} {unit}: more than one contig in chain")
                self.sgn, self.off = cmap(e, unit)
            self.b0s = [s["b0"] for s in self.segs]

    def disp(self, c):           # contig coordinate -> B73 display coordinate
        return self.off + self.sgn * c

    def pieces(self, lo, hi):    # [(b_lo, b_hi, x_at_lo, x_at_hi, inner)] for the B73 interval lo..hi
        if self.segs is None: return [(lo, hi, lo, hi, False)]
        out = []
        i = max(bisect.bisect_right(self.b0s, lo) - 1, 0)
        i = max(i - 50, 0)                                            # segments may overlap slightly in B73
        for s in self.segs[i:]:
            if s["b0"] >= hi: break
            x0, x1 = max(lo, s["b0"]), min(hi, s["b1"])
            if x1 > x0: out.append((x0, x1, self.disp(interp(s, x0)), self.disp(interp(s, x1)), s["zone"] == "inner"))
        return out


def correspond(A, B, w0, w1):
    res = []
    for (q0, q1, a0, a1, ainv) in A.pieces(w0, w1):
        for (r0, r1, b0, b1, binv) in B.pieces(q0, q1):
            # positions on A at r0, r1 by interpolation within A's piece
            fa = lambda b: a0 + (b - q0) / max(q1 - q0, 1) * (a1 - a0)
            res.append([r0, r1, fa(r0), fa(r1), b0, b1, ainv or binv])
    res.sort(key=lambda p: p[0])
    return res


def merge(ps, tol):
    out = []
    for p in ps:
        if out:
            o = out[-1]
            same_dir = ((p[3] - p[2]) * (o[3] - o[2]) >= 0) and ((p[5] - p[4]) * (o[5] - o[4]) >= 0)
            if (p[6] == o[6] and same_dir and p[0] - o[1] <= tol and abs(p[2] - o[3]) <= tol and abs(p[4] - o[5]) <= tol):
                o[1], o[3], o[5] = p[1], p[3], p[5]
                o[7] += 1
                continue
        out.append(p + [1])
    return out


# ---------- layout (mm) ----------
W = a.width
LAB_R = 24.0                                   # right edge of the label column
COL_X0, COL_X1 = 26.5, W - 4.0
GAP_SAME, GAP_DIFF = 6.0, 7.0                  # between windows of one chromosome / different chromosomes
nL = len(LOCI)
gaps = [GAP_SAME if LOCI[i]["chrom"] == LOCI[i + 1]["chrom"] else GAP_DIFF for i in range(nL - 1)]
CW = (COL_X1 - COL_X0 - sum(gaps)) / nL
XS = [COL_X0 + i * CW + sum(gaps[:i]) for i in range(nL)]
TH = 1.6                                       # track height
PITCH = 11.5                                   # track-to-track distance
Y_OV = 8.0                                     # overview bars
NOTE_H = 2.4 if any(L.get("note") for L in LOCI) else 0.0   # optional per-window note line under the interval label
if a.headings == "show":
    Y_HEAD = 23.0                              # chromosome headings
    Y_GUIDE = Y_HEAD - 3.2                     # lower end of the zoom guide lines
    Y_IVL = 26.6                               # interval labels (text bottom)
    Y_AX = 30.6 + NOTE_H                       # axis line
elif a.ideogram == "genome":                   # compact layout: no headings; interval label just above the ticks
    Y_GUIDE = 16.5
    Y_IVL = 19.8
    Y_AX = 23.1 + NOTE_H
else:                                          # focus ideogram: larger bars, axis pulled up
    Y_OV = 6.6                                 # top of the enlarged chromosome bars
    Y_AX = 21.2 + NOTE_H
    Y_IVL = Y_AX - 3.0
    Y_GUIDE = None                             # connectors end on the B73 track (set per window)
Y_B73 = Y_AX + 0.6
ASM_ORDER = {r["unit"]: i for i, r in enumerate(ASM)}
ROWS = [dict(unit="B73", line="B73", partition="", source="reference")] + ASM
Y_ROW = [Y_B73 + i * PITCH for i in range(len(ROWS))]
Y_LEG = Y_ROW[-1] + TH + 8.0
_nleg = len({r["row"] for r in tsv(LEGF)}) if os.path.exists(LEGF) else 0
H = Y_LEG + _nleg * 3.1 + 2.0

fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
RIB, TRK, GUIDES = [], [], []
TRK_TAGX = {}


def ribbon(x_top0, x_top1, yt, x_bot0, x_bot1, yb, col, alpha, clip):
    m = (yt + yb) / 2
    v = [(x_top0, yt), (x_top0, m), (x_bot0, m), (x_bot0, yb), (x_bot1, yb), (x_bot1, m), (x_top1, m), (x_top1, yt), (x_top0, yt)]
    cd = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    p = PathPatch(Path(v, cd), fc=col, ec="none", alpha=alpha, lw=0, zorder=1); ax.add_patch(p); p.set_clip_path(clip)


def mbl(x, step):
    d = 0 if step >= 1e6 else 1
    return f"{x / 1e6:.{d}f}"


# ---------- B73 overview ----------
OV_K = 2.3                                     # focus: enlargement of bar height, marker and overhang (genome ideogram = 1)
def ivl_edges(iv):                             # outer edges of the drawn B73 interval (inner block + breakpoint intervals)
    return (min(int(iv["b73_inner_start"]), int(iv["b73_outer_start"])), max(int(iv["b73_inner_end"]), int(iv["b73_outer_end"])))
if a.ideogram == "focus":
    fchr = list(dict.fromkeys(L["chrom"] for L in LOCI))
    span = {c: [min(x for L, x in zip(LOCI, XS) if L["chrom"] == c), max(x + CW for L, x in zip(LOCI, XS) if L["chrom"] == c)] for c in fchr}
    GAPF = 10.0                                # minimum gap between chromosome bars
    # one common Mb scale (true relative lengths). Each bar is centred over its panels, then pushed inside the
    # column and apart; the scale is the largest at which this works.
    def place(sc):
        xs, prev = {}, COL_X0 - GAPF
        for c in fchr:
            w = B73LEN[c] * sc; x = (span[c][0] + span[c][1]) / 2 - w / 2
            x = max(x, prev + GAPF, COL_X0); x = min(x, COL_X1 - w); xs[c] = x; prev = x + w
        ok = all(xs[c] >= COL_X0 - 1e-9 for c in fchr) and all(xs[b] - (xs[a_] + B73LEN[a_] * sc) >= GAPF - 1e-9 for a_, b in zip(fchr, fchr[1:]))
        return ok, xs
    lo_, hi_ = 0.0, (COL_X1 - COL_X0) / sum(B73LEN[c] for c in fchr)
    for _ in range(60):
        mid = (lo_ + hi_) / 2
        if place(mid)[0]: lo_ = mid
        else: hi_ = mid
    ov_scale = lo_; ov_x = place(ov_scale)[1]
    BH = 1.3 * OV_K
    for c in fchr:
        ax.add_patch(Rectangle((ov_x[c], Y_OV), B73LEN[c] * ov_scale, BH, fc="#d6d6d6", ec="none", zorder=2))
        ax.text(ov_x[c] + B73LEN[c] * ov_scale / 2, Y_OV - 0.35 * OV_K - 0.6, f"Chromosome {c[3:]}", fontsize=7, ha="center",
                va="bottom", color=INK, fontweight="bold")
    ax.text(ov_x[fchr[0]] - 1.5, Y_OV + BH / 2, "B73 v5", fontsize=5.5, ha="right", va="center", color=INK)
    OVMARK = {}
    for L, x0 in zip(LOCI, XS):
        iv = IVL[L["event"]]; c = L["chrom"]; col = COL.get(c, DEFAULT_COL)
        m0 = ov_x[c] + int(iv["b73_inner_start"]) * ov_scale; m1 = ov_x[c] + int(iv["b73_inner_end"]) * ov_scale
        if m1 - m0 < 0.6 * OV_K: mc = (m0 + m1) / 2; m0, m1 = mc - 0.3 * OV_K, mc + 0.3 * OV_K
        ax.add_patch(Rectangle((m0, Y_OV - 0.35 * OV_K), m1 - m0, 2.0 * OV_K, fc=col, ec="none", zorder=3))
        OVMARK[L["event"]] = (m0, m1, Y_OV - 0.35 * OV_K + 2.0 * OV_K)
chrs = [f"chr{i}" for i in range(1, 11)] if a.ideogram == "genome" else []
OV_GAP = 2.2
OV_FS = 5 if a.headings == "show" else 5.5    # same size as the row labels
if a.ideogram == "genome":
    ov_scale = (COL_X1 - COL_X0 - OV_GAP * 9) / sum(B73LEN[c] for c in chrs)
    ov_x = {}
x = COL_X0
for c in chrs:
    ov_x[c] = x
    ax.add_patch(Rectangle((x, Y_OV), B73LEN[c] * ov_scale, 1.3, fc="#d6d6d6", ec="none", zorder=2))
    shown = c in {L["chrom"] for L in LOCI}
    ax.text(x + B73LEN[c] * ov_scale / 2, Y_OV - 0.7, c[3:], fontsize=OV_FS, ha="center", va="bottom", color=INK,
            fontweight="bold" if (shown and a.headings == "hide") else "normal")
    x += B73LEN[c] * ov_scale + OV_GAP
if a.ideogram == "genome": ax.text(LAB_R, Y_OV + 0.65, "B73 v5", fontsize=OV_FS, ha="right", va="center", color=INK)
for L, x0 in (zip(LOCI, XS) if a.ideogram == "genome" else []):
    iv = IVL[L["event"]]; c = L["chrom"]; col = COL.get(c, DEFAULT_COL)
    m0 = ov_x[c] + int(iv["b73_inner_start"]) * ov_scale; m1 = ov_x[c] + int(iv["b73_inner_end"]) * ov_scale
    if m1 - m0 < 0.6: mc = (m0 + m1) / 2; m0, m1 = mc - 0.3, mc + 0.3
    ax.add_patch(Rectangle((m0, Y_OV - 0.35), m1 - m0, 2.0, fc=col, ec="none", zorder=3))
    for xa, xb in ((m0, x0), (m1, x0 + CW)):
        ax.plot([xa, xb], [Y_OV + 1.75, Y_GUIDE], color=col, lw=0.35, ls=(0, (2.2, 1.6)), zorder=1)
    if a.headings == "show":                   # two exterior lines only (no horizontal closing line)
        ax.plot([x0, x0 + CW], [Y_GUIDE] * 2, color=col, lw=0.35, ls=(0, (2.2, 1.6)), zorder=1)

# ---------- chromosome headings ----------
groups = []
for L, x0 in zip(LOCI, XS):
    if groups and groups[-1][0] == L["chrom"]: groups[-1][2] = x0 + CW
    else: groups.append([L["chrom"], x0, x0 + CW, L["heading"]])
for c, xa, xb, head in (groups if a.headings == "show" else []):
    ax.text((xa + xb) / 2, Y_HEAD, head, fontsize=6, ha="center", va="bottom", fontweight="bold", color=INK)

# ---------- windows ----------
for li, (L, x0) in enumerate(zip(LOCI, XS)):
    e, ch = L["event"], L["chrom"]; s, t = int(L["window_start"]), int(L["window_end"]); step = int(L["tick_step"])
    col = COL.get(ch, DEFAULT_COL); sc = CW / (t - s); X = lambda b: x0 + (b - s) * sc
    clip = Rectangle((x0, 0), CW, H, transform=ax.transData)
    iv = IVL[e]; bi0, bi1 = int(iv["b73_inner_start"]), int(iv["b73_inner_end"])
    bo0, bo1 = int(iv["b73_outer_start"]), int(iv["b73_outer_end"])
    if a.ideogram == "genome":
        ax.text(x0 + CW / 2, Y_IVL, f"{bi0 / 1e6:.3f}–{bi1 / 1e6:.3f} Mb", fontsize=5, ha="center", va="bottom", color=INK)
    if L.get("note"):
        ax.text(x0 + CW / 2, Y_IVL + NOTE_H, L["note"], fontsize=4.5, ha="center", va="bottom", color=INK2, style="italic")
    # axis
    ax.plot([x0, x0 + CW], [Y_AX] * 2, color=INK, lw=HAIR)
    for b in range((s // step + 1) * step, t + 1, step):
        if b == s: continue
        ax.plot([X(b)] * 2, [Y_AX, Y_AX - 0.5], color=INK, lw=HAIR)
        tk = ax.text(X(b), Y_AX - 0.65, mbl(b, step), fontsize=5, ha="center", va="bottom", color=INK, zorder=4)
        if a.ideogram == "focus":                  # connectors pass behind tick numbers
            tk.set_bbox(dict(fc="white", ec="none", pad=0.4))
    if a.ideogram == "focus":                  # connectors: marker ends -> outer edges of the drawn B73 interval
        m0, m1, ym = OVMARK[e]; q0, q1 = ivl_edges(iv)
        for xa, xb in ((m0, X(q0)), (m1, X(q1))):
            ax.plot([xa, xb], [ym, Y_B73], color=col, lw=0.35, ls=(0, (2.2, 1.6)), zorder=1, solid_capstyle="butt")
            GUIDES.append(dict(event=e, x_marker=round(xa, 4), y_marker=round(ym, 4), x_end=round(xb, 4), y_end=round(Y_B73, 4)))
        # interval label: centred between the two connectors at its baseline height, checked for clearance
        xl = lambda y: m0 + (X(q0) - m0) * (y - ym) / (Y_B73 - ym); xr = lambda y: m1 + (X(q1) - m1) * (y - ym) / (Y_B73 - ym)
        lt = ax.text(0, Y_IVL, f"{bi0 / 1e6:.3f}–{bi1 / 1e6:.3f} Mb", fontsize=5, ha="center", va="bottom", color=INK, zorder=4)
        bb = lt.get_window_extent(renderer=fig.canvas.get_renderer()); inv_ = ax.transData.inverted()
        (bx0, by0), (bx1, by1) = inv_.transform((bb.x0, bb.y0)), inv_.transform((bb.x1, bb.y1))
        tw, ttop = abs(bx1 - bx0), min(by0, by1)
        cx = (xl(Y_IVL) + xr(Y_IVL)) / 2; lt.set_x(cx)
        room = min(xr(y) for y in (ttop, Y_IVL)) - max(xl(y) for y in (ttop, Y_IVL))
        GUIDES[-1]["label_clearance_mm"] = GUIDES[-2]["label_clearance_mm"] = round(room - tw, 3)
        if room - tw < 0.5:                        # too narrow: connectors pass behind a white-backed label
            lt.set_bbox(dict(fc="white", ec="none", pad=0.6))
            GUIDES[-1]["label_backed"] = GUIDES[-2]["label_backed"] = 1
    if a.ideogram == "genome" and li + 1 < nL and LOCI[li + 1]["chrom"] == ch:   # break between windows of one chromosome
        xm = x0 + CW + gaps[li] / 2
        for yy in (Y_AX, Y_B73 + TH / 2):
            ax.text(xm, yy, "//", fontsize=6, ha="center", va="center", color=INK2)
    tracks = [Track(e, r["unit"]) for r in ROWS]
    # B73 track
    ax.add_patch(Rectangle((x0, Y_B73), CW, TH, fc=col, ec="none", alpha=0.55, zorder=2))
    ax.add_patch(Rectangle((X(bi0), Y_B73), X(bi1) - X(bi0), TH, fc=INV, ec="none", alpha=0.85, zorder=3))
    for u0, u1 in ((bo0, bi0), (bi1, bo1)):                          # breakpoint intervals (uncertainty)
        if u1 > u0:
            ax.add_patch(Rectangle((X(u0), Y_B73), X(u1) - X(u0), TH, fc="white", ec="none", zorder=3.4))
            ax.add_patch(Rectangle((X(u0), Y_B73), X(u1) - X(u0), TH, fc=INV, ec="none", alpha=0.35, zorder=3.5))
    # contig tracks
    for k, (r, T) in enumerate(zip(ROWS, tracks)):
        if k == 0: continue
        y = Y_ROW[k]; u = r["unit"]
        if not T.segs:
            ax.text(x0 + CW / 2, y + TH / 2, "no contig in window", fontsize=4.5, ha="center", va="center", color=INK2)
            TRK.append(dict(event=e, unit=u, contig="", note="no chain")); continue
        Lc = CLEN.get(u, {}).get(T.contig)
        ends = [T.disp(0), T.disp(Lc)] if Lc else [-1e18, 1e18]
        lo, hi = max(min(ends), s), min(max(ends), t)
        ax.add_patch(Rectangle((X(lo), y), X(hi) - X(lo), TH, fc=CONTIG, ec="none", zorder=2))
        for sg in T.segs:
            q0, q1 = sorted((T.disp(sg["c0"]), T.disp(sg["c1"])))
            if q1 < s or q0 > t: continue
            q0, q1 = max(q0, s), min(q1, t)
            ax.add_patch(Rectangle((X(q0), y), max(X(q1) - X(q0), 0.06), TH, fc=INV if sg["zone"] == "inner" else col,
                                   ec="none", alpha=0.85 if sg["zone"] == "inner" else 0.6, lw=0, zorder=3))
        cap = []
        for ce in ends:
            if s < ce < t:
                ax.add_patch(Rectangle((X(ce) - 0.15, y - 0.45), 0.3, TH + 0.9, fc=INK, ec="none", zorder=5)); cap.append(int(ce))
        for f0, f1, fc in FEAT.get((e, u), []):
            q0, q1 = sorted((T.disp(f0), T.disp(f1)))
            ax.add_patch(Rectangle((X(max(q0, s)), y - 0.75), X(min(q1, t)) - X(max(q0, s)), 0.55, fc=fc, ec="none", zorder=4))
        tags = []
        if T.sgn < 0: tags.append("rc")
        if STATE.get((e, u), {}).get("flag") and a.flag_tag == "show": tags.append("†")
        tag_x = (X(hi) if a.tag_anchor == "contig" else x0 + CW) + 0.5
        if a.tag_anchor == "contig" and hi < t: tag_x += 0.15      # clear the half-width of a contig-end bar
        if tags: TRK_TAGX[(e, u)] = round(tag_x, 4)
        if tags: ax.text(tag_x, y + TH / 2, " ".join(tags), fontsize=4.5, ha="left", va="center", color=INK2)
        for j in JR.get((e, u), []):
            if j.get("hidden") or a.junction_marks == "hide": continue
            if j["contig_lo"] is None:                       # unassessed: text only, never drawn as 0
                ax.text(x0 + (0.02 if j["side"] == "left" else 0.98) * CW, y + TH + 1.5, "n.a.", fontsize=4.5,
                        ha="center", va="top", color=INK2); continue
            xm = X(T.disp((j["contig_lo"] + j["contig_hi"]) / 2))
            if x0 <= xm <= x0 + CW:
                ax.add_patch(Polygon([(xm - 0.5, y + TH + 1.05), (xm + 0.5, y + TH + 1.05), (xm, y + TH + 0.15)], closed=True,
                                     fc=INK if (j["reads"] or 0) > 0 else "white", ec=INK, lw=0.3, zorder=6))
                if "label" in j:
                    later = j.get("same") and ASM_ORDER.get(j["same"], 1e9) < ASM_ORDER.get(u, 0)
                    ax.text(xm, y + TH + 1.25, "=" if later else j["label"], fontsize=4.5, ha="center", va="top", color=INK, zorder=6)
        TRK.append(dict(event=e, unit=u, contig=T.contig, contig_len=Lc or "", strand_drawn="rc" if T.sgn < 0 else "fwd",
                        offset=round(T.off), state=STATE.get((e, u), {}).get("state", ""),
                        flag=STATE.get((e, u), {}).get("flag", ""), contig_ends_in_window=";".join(map(str, cap)),
                        junctions=";".join(f"{j['side']}:{j['contig_lo']}-{j['contig_hi']}:{'hidden' if j.get('hidden') else j.get('label', j['reads'])}:{j['source']}"
                                           for j in JR.get((e, u), []))))
    # ribbons: each track to the nearest track above that has a chain
    tol = MERGE_FRAC * (t - s)
    for k in range(1, len(ROWS)):
        B = tracks[k]
        if not B.segs: continue
        ka = max(i for i in range(k) if i == 0 or tracks[i].segs)
        A = tracks[ka]
        ps = merge(correspond(A, B, s, t), tol)
        for p in ps:
            b_lo, b_hi, a_lo, a_hi, c_lo, c_hi, inner, n = p
            if b_hi - b_lo <= 0: continue
            ribbon(X(a_lo), X(a_hi), Y_ROW[ka] + TH, X(c_lo), X(c_hi), Y_ROW[k], INV if inner else col,
                   0.45 if inner else 0.22, clip)
            RIB.append(dict(event=e, upper=ROWS[ka]["unit"], lower=ROWS[k]["unit"], b73_start=b_lo, b73_end=b_hi,
                            upper_contig=A.contig or "B73", upper_disp_start=round(a_lo), upper_disp_end=round(a_hi),
                            lower_contig=B.contig, lower_disp_start=round(c_lo), lower_disp_end=round(c_hi),
                            crossing=int((a_hi - a_lo) * (c_hi - c_lo) < 0), inner=int(inner), pieces=n))

# ---------- row labels ----------
by_line = collections.OrderedDict()
for k, r in enumerate(ROWS): by_line.setdefault(r["line"], []).append(k)
for ln, ks in by_line.items():
    if len(ks) == 1 and not ROWS[ks[0]]["partition"]:
        ax.text(LAB_R, Y_ROW[ks[0]] + TH / 2, ln, fontsize=5.5, ha="right", va="center", color=INK,
                fontweight="bold" if ln == "B73" else "normal")
    else:
        for k in ks: ax.text(LAB_R, Y_ROW[k] + TH / 2, ROWS[k]["partition"], fontsize=5, ha="right", va="center", color=INK2)
        ya, yb = Y_ROW[ks[0]] + TH / 2, Y_ROW[ks[-1]] + TH / 2; xb = LAB_R - 5.4
        ax.plot([xb + 0.6, xb, xb, xb + 0.6], [ya, ya, yb, yb], color=INK2, lw=0.3)
        ax.text(xb - 0.6, (ya + yb) / 2, ln, fontsize=5.5, ha="right", va="center", color=INK)

# ---------- legend (config/legend.tsv: row, kind, text; items of one row are laid out left to right) ----------
lx, ly, dy = COL_X0, Y_LEG, 3.1
REN = fig.canvas.get_renderer()
def textw(tobj):
    bb = tobj.get_window_extent(renderer=REN); inv = ax.transData.inverted()
    return abs(inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0])
def sym(kind, x, y):
    if kind == "aligned":
        cs = [COL[c] for c in dict.fromkeys(L["chrom"] for L in LOCI) if c in COL] or [DEFAULT_COL]
        for i, c in enumerate(cs): ax.add_patch(Rectangle((x + i * 1.8, y - 0.7), 1.6, 1.4, fc=c, alpha=0.6, ec="none"))
        return 1.8 * len(cs)
    if kind == "reversed": ax.add_patch(Rectangle((x, y - 0.7), 2.4, 1.4, fc=INV, alpha=0.85, ec="none")); return 2.4
    if kind == "bp_interval": ax.add_patch(Rectangle((x, y - 0.7), 2.4, 1.4, fc=INV, alpha=0.35, ec="none")); return 2.4
    if kind == "contig": ax.add_patch(Rectangle((x, y - 0.7), 2.4, 1.4, fc=CONTIG, ec="none")); return 2.4
    if kind == "contig_end": ax.add_patch(Rectangle((x + 0.6, y - 1.0), 0.3, 2.0, fc=INK, ec="none")); return 1.5
    if kind == "triangle":
        ax.add_patch(Polygon([(x, y + 0.5), (x + 1.0, y + 0.5), (x + 0.5, y - 0.4)], closed=True, fc=INK, ec=INK, lw=0.3)); return 1.0
    if kind == "triangle_open":
        ax.add_patch(Polygon([(x, y + 0.5), (x + 1.0, y + 0.5), (x + 0.5, y - 0.4)], closed=True, fc="white", ec=INK, lw=0.3)); return 1.0
    if kind == "rc": ax.text(x, y, "rc", va="center", color=INK2, fontsize=4.5); return 1.6
    if kind.startswith("feature:"): ax.add_patch(Rectangle((x, y - 0.3), 2.4, 0.55, fc=kind.split(":", 1)[1], ec="none")); return 2.4
    return 0.0
LEG = tsv(LEGF) if os.path.exists(LEGF) else []
rows_ = collections.OrderedDict()
for r in LEG: rows_.setdefault(int(r["row"]), []).append(r)
for i, (rn, items) in enumerate(rows_.items()):
    y = ly + i * dy; x = lx
    for r in items:
        w_ = sym(r["kind"], x, y)
        t_ = ax.text(x + w_ + (0.8 if w_ else 0), y, r["text"], va="center", color=INK2 if r["kind"] == "note" else INK)
        x += w_ + (0.8 if w_ else 0) + textw(t_) + 4.0
fig.savefig(a.out + ".svg"); fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=300)
os.makedirs(TAB, exist_ok=True)
TAGS = [dict(event=k[0], unit=k[1], tag_x_mm=v) for k, v in TRK_TAGX.items()]
for name, rows in (("ribbons", RIB), ("tracks", TRK)) + ((("guides", GUIDES),) if GUIDES else ()) + ((("tags", TAGS),) if TAGS else ()):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(f"{TAB}/{name}.tsv", "w") as fh:
        w_ = csv.DictWriter(fh, keys, delimiter="\t", restval=""); w_.writeheader(); w_.writerows(rows)
print(f"wrote {a.out}.svg/.pdf/.png ({W:.0f} x {H:.0f} mm); {len(RIB)} ribbons, {len(TRK)} tracks")
