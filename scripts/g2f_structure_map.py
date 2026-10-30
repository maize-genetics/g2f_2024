#!/usr/bin/env python3
"""Structural-confidence map of assemblies: AnchorWave synteny classes combined with read evidence.

Takes a finished g2f_riparian.py run (sequences.tsv, ribbons.tsv, RUN.json) and, optionally, read evidence, and
classifies every anchored record into regions:

  consistent        drawn ribbon on the record's chromosome, in the record's majority orientation
  inversion         drawn ribbon on the record's chromosome in the opposite orientation (inverted vs. B73)
  off_chromosome    drawn ribbon on another chromosome (terminal = within --terminal-bp of a record end)
  ambiguous         ribbons of different chromosomes overlapping on the record, a record whose assignment differs
                    from its dominant chromosome, or a record with < --min-dominant-frac of anchors on one chromosome
  redundant         consistent ribbon whose reference interval is >= --redundant-frac covered by a ribbon of a
                    record with more anchors in the same unit (extra allele, duplicate or paralogue: a candidate)
  no_anchor         gaps >= --gap-bp between ribbons on anchored records, and whole unanchored records >= --gap-bp

Evidence (optional, TSV; coordinates 1-based, on the delivered records):
  unit, record, start, end, kind, call, source, detail
    kind  junction  a read test of continuity at a position (call: SUPPORTED / CONTINUOUS = reads span it;
                    CONTRADICTED = reads support another arrangement; BREAK = no read spans it;
                    UNRESOLVED / WEAK / anything else = unresolved)
          edit      a structural edit made at that position (call: e.g. broken / split / retained)
Depth windows (optional, TSV): unit, path (+ optional coordinate_map) pointing at g2f_asmqc
  regions/window_classes.tsv[.gz] (contig, start, end, depth_ratio, class); each region gets its median depth ratio
  and the fraction of windows in each depth class.

Synteny never decides alone. Each inversion / off-chromosome region gets `read_evidence` from the junction calls at
its boundaries and an `interpretation`:
  read-supported      reads span every tested boundary  -> reference difference (structural variant), not an error
  read-contradicted   reads contradict a boundary       -> assembly error candidate
  not-spanned         no read crosses a boundary (BREAK)  -> candidate problem; absence of support, not contradiction
  unresolved          otherwise; recurrence of the same reference breakpoints in other units is reported as context
                      (`recurrence`), never as proof

Outputs (new directory): flagged_regions.tsv (all non-consistent regions + junction/edit rows), regions.tsv (every
region incl. consistent), unit_summary.tsv, RUN.json, METHODS.md, figures/structure_<group>.{pdf,svg,png},
figures/structure_by_unit.pdf.

Usage:
  g2f_structure_map.py run --riparian RIPDIR --out OUTDIR [--evidence EV.tsv] [--windows WIN.tsv]
"""
import argparse
import csv
import datetime
import gzip
import json
import os
import sys
import textwrap
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2f_riparian as rip  # noqa: E402

CODE_VERSION = 1
SUPPORTED = {"SUPPORTED", "CONTINUOUS"}
CONTRADICTED = {"CONTRADICTED"}  # reads positively support another arrangement (e.g. split reads to a partner locus)
NOT_SPANNED = {"BREAK"}  # no read crosses some point: absence of support, not contradiction
RANK = ["read-contradicted", "not-spanned", "unresolved", "read-supported"]
# region colours (status-like meaning, always paired with a legend label and a marker shape)
CLASS_COL = {"consistent": "#c9d8c5", "inversion": "#6250d6", "off_chromosome": "#c0262d", "ambiguous": "#eda100",
             "redundant": "#eb6834", "no_anchor": "#ffffff"}
CALL_COL = {"read-supported": "#1f7a3a", "read-contradicted": "#c0262d", "not-spanned": "#d4540f",
            "unresolved": "#b07800", "untested": "#8a8984"}
INK, INK2, MUTED = rip.INK, rip.INK2, rip.MUTED


def die(msg):
    sys.exit(f"g2f_structure_map: error: {msg}")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_run(ripdir):
    for f in ("sequences.tsv", "ribbons.tsv", "RUN.json"):
        if not os.path.isfile(os.path.join(ripdir, f)):
            die(f"{ripdir}: missing {f} (not a g2f_riparian run?)")
    with open(os.path.join(ripdir, "RUN.json")) as fh:
        run = json.load(fh)
    seqs = rip.read_tsv(os.path.join(ripdir, "sequences.tsv"))
    ribs = rip.read_tsv(os.path.join(ripdir, "ribbons.tsv"))
    for r in seqs:
        for k in ("length_bp", "gene_anchors", "anchors_on_dominant"):
            r[k] = int(r[k]) if r.get(k) not in (None, "") else 0
        r["frac_plus_on_dominant"] = float(r["frac_plus_on_dominant"]) if r.get("frac_plus_on_dominant") else 1.0
        r["frac_anchors_on_dominant"] = (float(r["frac_anchors_on_dominant"]) if r.get("frac_anchors_on_dominant")
                                         else 0.0)
    for r in ribs:
        for k in ("ref_start", "ref_end", "query_start", "query_end", "gene_anchors"):
            r[k] = int(r[k])
    return run, seqs, ribs


def read_evidence(path, known):
    ev = defaultdict(list)
    for r in rip.read_tsv(path):
        for k in ("unit", "record", "start", "end", "kind", "call"):
            if not r.get(k):
                die(f"{path}: row lacks {k}: {r}")
        if (r["unit"], r["record"]) not in known:
            die(f"{path}: {r['unit']} {r['record']} is not a record of the riparian run")
        if r["kind"] not in ("junction", "edit"):
            die(f"{path}: kind must be junction or edit, not {r['kind']}")
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        if r["end"] < r["start"] or r["end"] > known[(r["unit"], r["record"])]:
            die(f"{path}: bad interval {r['record']}:{r['start']}-{r['end']}")
        ev[(r["unit"], r["record"])].append(r)
    return ev


def read_windows(path, units_by_id, lengths):
    """unit -> record -> sorted [(start, end, depth_ratio, class)] on delivered-record coordinates (1-based)."""
    out = {}
    for spec in rip.read_tsv(path):
        u = spec["unit"]
        if u not in units_by_id:
            die(f"{path}: unknown unit {u}")
        cmap = rip.read_coordinate_map(spec["coordinate_map"]) if spec.get("coordinate_map") else None
        per = defaultdict(list)
        with opener(spec["path"]) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {k: hdr.index(k) for k in ("contig", "start", "end", "depth_ratio", "class")}
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                c, s, e = f[ix["contig"]], int(f[ix["start"]]) + 1, int(f[ix["end"]])
                if cmap and c in cmap:
                    mid = (s + e) // 2
                    for ps, pe, rec in cmap[c]:
                        if ps <= mid <= pe:
                            c, s, e = rec, max(s, ps) - ps + 1, min(e, pe) - ps + 1
                            break
                if (u, c) not in lengths:
                    continue  # windows on records the run does not know (should not happen for the same FASTA)
                try:
                    d = float(f[ix["depth_ratio"]])
                except ValueError:
                    d = float("nan")
                per[c].append((s, e, d, f[ix["class"]]))
        for c in per:
            per[c].sort()
        out[u] = per
    return out


def window_stats(wins, s, e):
    if not wins:
        return {}
    import bisect
    i = bisect.bisect_left(wins, (s - 20_000,))
    sel = []
    while i < len(wins) and wins[i][0] <= e:
        if wins[i][1] >= s:
            sel.append(wins[i])
        i += 1
    if not sel:
        return {}
    ds = sorted(d for _, _, d, _ in sel if d == d)
    cls = defaultdict(int)
    for w in sel:
        cls[w[3]] += 1
    n = len(sel)
    return dict(depth_median=round(ds[len(ds) // 2], 3) if ds else "",
                depth_classes=";".join(f"{k}:{v / n:.2f}" for k, v in sorted(cls.items(), key=lambda kv: -kv[1])),
                frac_low=round((cls.get("low", 0) + cls.get("none", 0)) / n, 3),
                frac_half=round(cls.get("half", 0) / n, 3), frac_dup=round((cls.get("dup", 0) + cls.get("multi", 0))
                                                                           / n, 3))


def union(iv):
    return rip.union_len(iv)


def classify_unit(unit, seqs, ribs, P):
    """Regions for one unit. Returns (regions, boundaries) where boundaries are the synteny breakpoints between
    consecutive drawn ribbons on a record."""
    rows = {r["record"]: r for r in seqs}
    drawn = defaultdict(list)
    for r in ribs:
        if r["drawn"] == "yes":
            drawn[r["record"]].append(r)
    # redundancy at ribbon level: reference interval covered by a record with more anchors on that chromosome
    by_chr = defaultdict(list)
    for rec, rb in drawn.items():
        for r in rb:
            by_chr[r["ref_chr"]].append((rows[rec]["anchors_on_dominant"], rec, r["ref_start"], r["ref_end"]))
    regions, bounds = [], []
    for rec, row in rows.items():
        L = row["length_bp"]
        rb = sorted(drawn.get(rec, []), key=lambda r: r["query_start"])
        if not rb:
            if row["gene_anchors"] == 0 and L >= P["gap_bp"]:
                regions.append(dict(record=rec, start=1, end=L, cls="no_anchor", detail="record without gene anchors"))
            continue
        home = row.get("assigned_chr") or row["dominant_chr"]
        major = "+" if row["frac_plus_on_dominant"] >= 0.5 else "-"
        rec_flags = []
        if row.get("assignment_agrees") == "no":
            rec_flags.append(f"assigned {row['assigned_chr']} but dominant {row['dominant_chr']}")
        if row["frac_anchors_on_dominant"] < P["min_dominant_frac"]:
            rec_flags.append(f"only {row['frac_anchors_on_dominant']:.2f} of anchors on {row['dominant_chr']}")
        for r in rb:
            qs, qe = r["query_start"], r["query_end"]
            # a boundary is a junction only if synteny continues on the other side
            base = dict(record=rec, start=qs, end=qe, ref_chr=r["ref_chr"], ref_start=r["ref_start"],
                        ref_end=r["ref_end"], strand=r["strand"], gene_anchors=r["gene_anchors"],
                        _left=any(o["query_end"] < qs for o in rb), _right=any(o["query_start"] > qe for o in rb))
            if r["ref_chr"] != home:
                term = qs <= P["terminal_bp"] or qe >= L - P["terminal_bp"]
                regions.append(dict(base, cls="off_chromosome", subtype="terminal" if term else "internal",
                                    detail=f"{r['ref_chr']} block on a {home} record"))
            elif r["strand"] != major:
                regions.append(dict(base, cls="inversion", detail=f"{r['ref_chr']} block inverted vs. rest of record"))
            else:
                best, best_f = "", 0.0
                own = r["ref_end"] - r["ref_start"] + 1
                for n_other, orec, s, e in by_chr[r["ref_chr"]]:
                    if orec == rec or n_other <= row["anchors_on_dominant"]:
                        continue
                    ov = max(0, min(e, r["ref_end"]) - max(s, r["ref_start"]) + 1)
                    if ov / own > best_f:
                        best, best_f = orec, ov / own
                if best_f >= P["redundant_frac"]:
                    regions.append(dict(base, cls="redundant", partner=best,
                                        detail=f"{best_f * 100:.0f}% of {r['ref_chr']}:{rip.mb(r['ref_start'])}-"
                                               f"{rip.mb(r['ref_end'])} Mb also on {best}"))
                else:
                    regions.append(dict(base, cls="consistent", detail=""))
        # overlapping ribbons of different chromosomes = competing assignment of the same sequence
        for i, a in enumerate(rb):
            for b in rb[i + 1:]:
                if b["query_start"] > a["query_end"]:
                    break
                if b["ref_chr"] != a["ref_chr"]:
                    s, e = b["query_start"], min(a["query_end"], b["query_end"])
                    rec_flags_txt = f"{a['ref_chr']} and {b['ref_chr']} ribbons overlap on the record"
                    regions.append(dict(record=rec, start=s, end=e, cls="ambiguous", detail=rec_flags_txt))
        if rec_flags:
            regions.append(dict(record=rec, start=1, end=L, cls="ambiguous", detail="; ".join(rec_flags)))
        # gaps without anchors on the record (including the ends)
        cov, pos = sorted((r["query_start"], r["query_end"]) for r in rb), 1
        for s, e in cov + [(L + 1, L + 1)]:
            if s - pos >= P["gap_bp"]:
                regions.append(dict(record=rec, start=pos, end=s - 1, cls="no_anchor",
                                    detail="no drawn anchor block"))
            pos = max(pos, e + 1)
        # synteny breakpoints between consecutive ribbons (chromosome or strand change)
        for a, b in zip(rb, rb[1:]):
            if a["ref_chr"] != b["ref_chr"] or a["strand"] != b["strand"]:
                bounds.append(dict(record=rec, start=min(a["query_end"], b["query_start"]),
                                   end=max(a["query_end"], b["query_start"]),
                                   type="inter_chromosomal" if a["ref_chr"] != b["ref_chr"] else "strand_switch",
                                   left=f"{a['ref_chr']}:{a['ref_start']}-{a['ref_end']}{a['strand']}",
                                   right=f"{b['ref_chr']}:{b['ref_start']}-{b['ref_end']}{b['strand']}"))
    for r in regions:
        r.setdefault("subtype", "")
    return regions, bounds


def call_class(call):
    c = call.upper()
    if c in SUPPORTED:
        return "read-supported"
    if c in CONTRADICTED:
        return "read-contradicted"
    if c in NOT_SPANNED:
        return "not-spanned"
    return "unresolved"


def near(ev_rows, s, e, slop):
    return [x for x in ev_rows if x["start"] <= e + slop and x["end"] >= s - slop]


def attach_evidence(unit, regions, bounds, ev, P):
    """Read evidence for each flagged region from junction calls at its boundaries."""
    for r in regions:
        if r["cls"] not in ("inversion", "off_chromosome"):
            continue
        rows = ev.get((unit, r["record"]), [])
        tested = []
        # only boundaries with synteny on both sides are junctions (a block running to a record end has one)
        for pos, internal in ((r["start"], r.get("_left")), (r["end"], r.get("_right"))):
            if not internal:
                continue
            hits = [x for x in near(rows, pos, pos, P["slop_bp"]) if x["kind"] == "junction"]
            edits = [x for x in near(rows, pos, pos, P["slop_bp"]) if x["kind"] == "edit"]
            note = ("; edit: " + "; ".join(x["call"] for x in edits)) if edits else ""
            if hits:
                cls = sorted({call_class(x["call"]) for x in hits}, key=RANK.index)[0]
                tested.append((cls, "; ".join(f"{x['call']} [{x.get('source', '')}] {x.get('detail', '')[:160]}"
                                              for x in hits) + note))
            elif edits:
                tested.append(("edited", note[2:]))
            else:
                tested.append(("untested", "no read test at this boundary"))
        classes = [t[0] for t in tested]
        if not tested:
            r["read_evidence"] = "no internal junction"
        elif "read-contradicted" in classes:
            r["read_evidence"] = "read-contradicted"
        elif "not-spanned" in classes:
            r["read_evidence"] = "not-spanned"
        elif all(c == "read-supported" for c in classes):
            r["read_evidence"] = "read-supported"
        elif any(c == "unresolved" for c in classes) or "read-supported" in classes:
            r["read_evidence"] = "unresolved"
        else:
            r["read_evidence"] = "untested" if all(c == "untested" for c in classes) else "edited"
        r["evidence_detail"] = " | ".join(f"{c}: {d}" for c, d in tested)


def recurrence(all_regions, P):
    """Inversion regions sharing reference boundaries (both within tol) in other units/lines."""
    inv = [(u, line, r) for u, line, rs in all_regions for r in rs if r["cls"] == "inversion"]
    for u, line, r in inv:
        others = sorted({l2 for u2, l2, r2 in inv if l2 != line and r2["ref_chr"] == r["ref_chr"]
                         and abs(r2["ref_start"] - r["ref_start"]) <= P["recur_tol_bp"]
                         and abs(r2["ref_end"] - r["ref_end"]) <= P["recur_tol_bp"]})
        r["recurrence"] = f"{len(others) + 1} lines" + (f" (also {','.join(others)})" if others else "")


def interpret(r):
    ev = r.get("read_evidence", "")
    rec = r.get("recurrence", "")
    multi = rec and not rec.startswith("1 line")
    if r["cls"] == "inversion":
        if ev == "read-supported":
            return "reference difference: inversion vs. B73 spanned by reads (structural variant, not an error)"
        if ev == "read-contradicted":
            return "assembly error candidate: reads contradict a boundary"
        if ev == "not-spanned":
            return ("candidate local assembly problem: no read spans a boundary (absence of support, not "
                    "contradiction)" + ("; same B73 breakpoints in other lines suggest the inversion itself is real"
                                         if multi else ""))
        if multi:
            return ("unresolved by reads; same B73 breakpoints in other lines suggest a real structural variant "
                    "(context, not proof)")
        return "unresolved: inversion vs. B73 without a conclusive read test"
    if r["cls"] == "off_chromosome":
        if ev == "read-supported":
            return ("read-supported join to another chromosome: candidate translocation / reference difference "
                    "(not an error by reads)")
        if ev == "read-contradicted":
            return "assembly error candidate: reads contradict the join"
        if ev == "not-spanned":
            return "candidate misjoin: no read spans the join (absence of support, not contradiction)"
        return "unresolved join to another chromosome (candidate misjoin; reads do not decide)"
    if r["cls"] == "redundant":
        return ("redundant-copy candidate (same reference genes on two records of this FASTA): allele, duplicate "
                "or paralogue; see depth")
    if r["cls"] == "ambiguous":
        return "competing chromosome signal; review"
    if r["cls"] == "no_anchor":
        return "no gene-anchored synteny (structure not assessable by AnchorWave)"
    return "chromosome-consistent"


def confidence(r):
    if r["cls"] != "consistent":
        return ""
    if r.get("frac_low", 0) and r["frac_low"] > 0.05:
        return "consistent; >5% low-depth windows"
    return "high" if r.get("depth_median", "") != "" else "consistent (no depth data)"


# ------------------------------------------------------------------------------------------------------------ plots
def plot(out, units, seqs_by_unit, regions_by_unit, ev, groups, chroms, P):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle, Patch
    from matplotlib.lines import Line2D
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir)
    col = dict(zip(chroms, rip.PALETTE))
    order_key = "assembly_order" if any(r.get("assembly_order") for s in seqs_by_unit.values() for r in s) \
        else "plot_order"
    group_key = "assembly_group" if order_key == "assembly_order" else "dominant_chr"
    MARK = {"read-supported": "^", "read-contradicted": "X", "not-spanned": "D", "unresolved": "o", "edit": "|"}

    def layout(u):
        rows = [r for r in seqs_by_unit[u] if r.get(order_key)]
        rows.sort(key=lambda r: int(r[order_key]))
        return rows

    def legend(fig, y):
        h = [Patch(fc=CLASS_COL[k], ec=MUTED if k in ("no_anchor", "consistent") else "none", lw=0.6,
                   label=k.replace("_", " ")) for k in CLASS_COL]
        h += [Line2D([], [], ls="none", marker=MARK[k], color=CALL_COL.get(k, INK), ms=7,
                     label={"read-supported": "junction: reads span", "read-contradicted": "junction: reads contradict",
                            "not-spanned": "junction: no read spans it",
                            "unresolved": "junction: unresolved / weak", "edit": "release edit (J1-J7 split/break)"}[k])
              for k in MARK]
        fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, y), ncol=10, fontsize=8, frameon=False,
                   handlelength=1.4, columnspacing=1.0)

    def draw_track(ax, y, u, rows, xmax, gapw, ggap, label, detail):
        regs = defaultdict(list)
        for r in regions_by_unit[u]:
            regs[r["record"]].append(r)
        x, prev, gx = 0.0, None, {}
        for row in rows:
            g = row.get(group_key) or "unassigned"
            if prev is not None:
                x += ggap if g != prev else gapw
            rec, L = row["record"], row["length_bp"]
            ax.add_patch(Rectangle((x, y - 0.3), L, 0.6, fc=CLASS_COL["no_anchor"], ec=MUTED, lw=0.3, zorder=2))
            for r in sorted(regs[rec], key=lambda r: ["no_anchor", "consistent", "redundant", "ambiguous",
                                                      "inversion", "off_chromosome"].index(r["cls"])):
                if r["cls"] == "no_anchor":
                    continue
                hh = 0.6 if r["cls"] not in ("ambiguous",) else 0.2
                yy = y - 0.3 if r["cls"] != "ambiguous" else y + 0.3
                ax.add_patch(Rectangle((x + r["start"], yy), max(r["end"] - r["start"], xmax * 4e-4), hh,
                                       fc=CLASS_COL[r["cls"]], ec="none", zorder=3 + (r["cls"] in
                                                                                       ("inversion", "off_chromosome"))))
            for e in ev.get((u, rec), []):
                k = "edit" if e["kind"] == "edit" else call_class(e["call"])
                ax.plot(x + (e["start"] + e["end"]) / 2, y + 0.55, ls="none", marker=MARK[k],
                        color=CALL_COL.get(k, INK), ms=5 if not detail else 6, zorder=6,
                        mec="white" if k != "edit" else None, mew=0.4 if k != "edit" else 1.4)
            gx.setdefault(g, [x, x + L])[1] = x + L
            x += L
            prev = g
        for g, (a, b) in gx.items():
            ax.plot([a, b], [y - 0.45, y - 0.45], color=col.get(g, MUTED), lw=2, solid_capstyle="butt")
            if True:
                ax.text((a + b) / 2, y - 0.55, g, ha="center", va="top", fontsize=6.5, color=INK2, style="italic")
        ax.text(-xmax * 0.005, y, label, ha="right", va="center", fontsize=8, color=INK)
        return x

    def width(u, gapw, ggap):
        rows, x, prev = layout(u), 0.0, None
        for row in rows:
            g = row.get(group_key) or "unassigned"
            if prev is not None:
                x += ggap if g != prev else gapw
            x += row["length_bp"]
            prev = g
        return x

    gapw, ggap = 1.5e6, 12e6
    xmax = max(width(u["unit"], gapw, ggap) for u in units)
    sub = textwrap.fill(
        "Records in " + ("assembly order, grouped by their own chromosome assignment (not a physical order)"
                         if order_key == "assembly_order" else "B73-synteny order") +
        "; bar = record, coloured by AnchorWave synteny class (white = no gene-anchored block). Thin yellow strip "
        "above a bar = ambiguous chromosome signal. Markers above bars = read tests of continuity at junctions "
        "(HiFi/ONT) and release edits. Synteny classes are candidates; read evidence decides between reference "
        "difference and assembly error.", 250)
    written = []
    for g, members in groups.items():
        fig = plt.figure(figsize=(20, 0.62 * len(members) + 1.9))
        H = 0.62 * len(members) + 1.9
        ax = fig.add_axes([0.07, 0.02, 0.92, 1 - 1.55 / H])
        for i, u in enumerate(members):
            draw_track(ax, -i * 1.6, u["unit"], layout(u["unit"]), xmax, gapw, ggap, u["unit"], detail=False)
        ax.set_xlim(0, xmax * 1.01)
        ax.set_ylim(-(len(members) - 1) * 1.6 - 1.0, 1.0)
        ax.axis("off")
        fig.suptitle(f"Structural-confidence map — {g}", x=0.01, ha="left", y=1 - 0.15 / H, fontsize=13, color=INK)
        fig.text(0.01, 1 - 0.42 / H, sub, fontsize=8, color=INK2, va="top")
        legend(fig, 1 - 0.95 / H)
        for ext in ("pdf", "svg", "png"):
            p = os.path.join(figdir, f"structure_{g}.{ext}")
            fig.savefig(p, dpi=150 if ext == "png" else None)
            written.append(p)
        plt.close(fig)
    # per unit: one row per chromosome group, with labels
    pdf_path = os.path.join(figdir, "structure_by_unit.pdf")
    with PdfPages(pdf_path) as pdf:
        for u in units:
            rows = layout(u["unit"])
            byg = defaultdict(list)
            for r in rows:
                byg[r.get(group_key) or "unassigned"].append(r)
            gs = [c for c in chroms if c in byg] + (["unassigned"] if "unassigned" in byg else [])
            gx = max(sum(r["length_bp"] for r in byg[g]) + gapw * (len(byg[g]) - 1) for g in gs)
            H = 0.75 * len(gs) + 2.0
            fig = plt.figure(figsize=(20, H))
            ax = fig.add_axes([0.06, 0.02, 0.93, 1 - 1.6 / H])
            for i, g in enumerate(gs):
                y = -i * 1.9
                draw_track(ax, y, u["unit"], byg[g], gx, gapw, ggap, g, detail=True)
                x = 0.0
                for r in byg[g]:
                    if r["length_bp"] >= 5e6:
                        ax.text(x + r["length_bp"] / 2, y, r["record"], ha="center", va="center", fontsize=5.5,
                                color=INK, zorder=7)
                    x += r["length_bp"] + gapw
            ax.set_xlim(0, gx * 1.01)
            ax.set_ylim(-(len(gs) - 1) * 1.9 - 1.0, 1.0)
            ax.axis("off")
            fig.suptitle(f"Structural-confidence map — {u['unit']} ({u.get('label', '')})", x=0.01, ha="left",
                         y=1 - 0.15 / H, fontsize=12, color=INK)
            fig.text(0.01, 1 - 0.4 / H, sub, fontsize=7.5, color=INK2, va="top")
            legend(fig, 1 - 1.0 / H)
            pdf.savefig(fig)
            plt.close(fig)
    written.append(pdf_path)
    return written


REGION_COLS = ["line", "unit", "record", "start", "end", "length_bp", "region_class", "subtype", "ref_chr",
               "ref_start", "ref_end", "strand", "gene_anchors", "partner", "detail", "read_evidence",
               "evidence_detail", "recurrence", "depth_median", "depth_classes", "frac_low", "frac_half", "frac_dup",
               "interpretation", "confidence"]

METHODS = """# Structural-confidence map: methods ({stamp})

Generated by `g2f_structure_map.py` (code version {cv}) from the g2f_riparian run `{rip}`. Parameters: RUN.json.

**Regions.** Each drawn AnchorWave ribbon of an anchored record (>= the riparian run's minimum anchors) becomes one
region. The record's home chromosome is its assigned chromosome when the run has an assignment, else its dominant
chromosome.
- `consistent`: home chromosome, the record's majority strand, not redundant.
- `inversion`: home chromosome, opposite strand.
- `off_chromosome`: another chromosome; `terminal` if within {terminal_bp:,} bp of a record end.
- `redundant`: a consistent ribbon whose B73 interval is >= {redundant_pct:.0f} % covered by a ribbon of a record
  with more anchors on that chromosome in the same unit.
- `ambiguous`: ribbons of two chromosomes overlapping on the record; assignment != dominant chromosome; or
  < {min_dominant_frac} of anchors on the dominant chromosome.
- `no_anchor`: >= {gap_bp:,} bp of a record without a drawn ribbon, and unanchored records >= {gap_bp:,} bp.

**Read evidence.** For inversion and off-chromosome regions, each internal boundary (not a record end) is matched to
evidence rows within {slop_bp:,} bp; boundaries without synteny on the far side (block runs to a record end) are not
junctions. SUPPORTED / CONTINUOUS = reads span it; CONTRADICTED = reads support another arrangement; BREAK = no read
crosses some point (`not-spanned`: absence of support, e.g. a repeat longer than the reads, not contradiction); other
calls = unresolved; no row = untested. The region is read-supported only if every tested boundary is supported;
otherwise the worst boundary decides (contradicted > not-spanned > unresolved).

**Recurrence.** Inversion regions of other lines whose B73 start and end both lie within {recur_tol_bp:,} bp. Shared
breakpoints in independently assembled lines make an assembly artefact unlikely, but this is context, not proof.

**Depth.** When asmqc windows are given, each region reports the median depth ratio (observed / expected) and the
fraction of windows per depth class. For haplotype FASTAs mapped one haplotype at a time, `half` is expected in
heterozygous sequence.

**Confidence.** `high` = consistent region with depth data and <= 5 % low/zero-depth windows. AnchorWave collinearity
alone is never called validated.
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "run"])
    ap.add_argument("--riparian", required=True, help="a finished g2f_riparian.py run directory")
    ap.add_argument("--evidence", help="junction / edit evidence TSV (see header)")
    ap.add_argument("--windows", help="TSV: unit, path (asmqc window_classes), optional coordinate_map")
    ap.add_argument("--out")
    ap.add_argument("--terminal-bp", type=int, default=2_000_000)
    ap.add_argument("--redundant-frac", type=float, default=0.5)
    ap.add_argument("--min-dominant-frac", type=float, default=0.8)
    ap.add_argument("--gap-bp", type=int, default=1_000_000)
    ap.add_argument("--slop-bp", type=int, default=100_000)
    ap.add_argument("--recur-tol-bp", type=int, default=500_000)
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    P = dict(terminal_bp=a.terminal_bp, redundant_frac=a.redundant_frac, min_dominant_frac=a.min_dominant_frac,
             gap_bp=a.gap_bp, slop_bp=a.slop_bp, recur_tol_bp=a.recur_tol_bp)
    run, seqs, ribs = read_run(a.riparian)
    units = run["units"]
    chroms = list(run["reference_chromosomes"])
    seqs_by_unit, ribs_by_unit = defaultdict(list), defaultdict(list)
    for r in seqs:
        seqs_by_unit[r["unit"]].append(r)
    for r in ribs:
        ribs_by_unit[r["unit"]].append(r)
    lengths = {(r["unit"], r["record"]): r["length_bp"] for r in seqs}
    ev = read_evidence(a.evidence, lengths) if a.evidence else {}
    wins = read_windows(a.windows, {u["unit"]: u for u in units}, lengths) if a.windows else {}
    all_regions, regions_by_unit = [], {}
    for u in units:
        regs, bounds = classify_unit(u["unit"], seqs_by_unit[u["unit"]], ribs_by_unit[u["unit"]], P)
        for r in regs:
            r["_len"] = lengths[(u["unit"], r["record"])]
        attach_evidence(u["unit"], regs, bounds, ev, P)
        for r in regs:
            r.update(window_stats(wins.get(u["unit"], {}).get(r["record"]), r["start"], r["end"]))
        regions_by_unit[u["unit"]] = regs
        all_regions.append((u["unit"], u["line"], regs))
    recurrence(all_regions, P)
    out_rows = []
    for u in units:
        for r in regions_by_unit[u["unit"]]:
            r["interpretation"] = interpret(r)
            r["confidence"] = confidence(r)
            out_rows.append(dict(r, line=u["line"], unit=u["unit"], region_class=r["cls"],
                                 length_bp=r["end"] - r["start"] + 1))
        # junction / edit evidence rows that no flagged region claimed are listed too
        for (uu, rec), rows in ev.items():
            if uu != u["unit"]:
                continue
            for e in rows:
                k = "edit" if e["kind"] == "edit" else "junction"
                out_rows.append(dict(line=u["line"], unit=uu, record=rec, start=e["start"], end=e["end"],
                                     length_bp=e["end"] - e["start"] + 1, region_class=k,
                                     subtype=e.get("type", ""), detail=e.get("detail", ""),
                                     read_evidence="edit" if k == "edit" else call_class(e["call"]),
                                     evidence_detail=f"{e['call']} [{e.get('source', '')}]",
                                     interpretation=(f"release edit: {e['call']}" if k == "edit" else
                                                     {"read-supported": "reads span this junction",
                                                      "read-contradicted": "reads contradict this junction",
                                                      "not-spanned": "no read spans this junction",
                                                      "unresolved": "read test inconclusive"}[call_class(e["call"])])))
    if a.command == "check":
        n = defaultdict(int)
        for r in out_rows:
            n[r["region_class"]] += 1
        print(f"check OK: {dict(n)} (nothing written)", file=sys.stderr)
        return 0
    if not a.out:
        die("run needs --out")
    if os.path.exists(a.out) or os.path.exists(a.out + ".partial"):
        die(f"{a.out} (or .partial) exists; choose a new output directory")
    start = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    partial = a.out + ".partial"
    os.makedirs(partial)
    rip.write_tsv(os.path.join(partial, "regions.tsv"), out_rows, REGION_COLS)
    rip.write_tsv(os.path.join(partial, "flagged_regions.tsv"),
                  [r for r in out_rows if r["region_class"] != "consistent"], REGION_COLS)
    summ = []
    for u in units:
        rs = [r for r in out_rows if r["unit"] == u["unit"]]
        s = dict(line=u["line"], unit=u["unit"])
        for c in ("consistent", "inversion", "off_chromosome", "ambiguous", "redundant", "no_anchor"):
            sel = [r for r in rs if r["region_class"] == c]
            s[f"{c}_n"] = len(sel)
            s[f"{c}_Mb"] = round(sum(r["length_bp"] for r in sel) / 1e6, 2)
        for c in ("read-supported", "read-contradicted", "not-spanned", "unresolved", "untested"):
            s[f"flagged_{c}"] = sum(r.get("read_evidence") == c for r in rs
                                    if r["region_class"] in ("inversion", "off_chromosome"))
        s["junction_tests"] = sum(r["region_class"] == "junction" for r in rs)
        s["junctions_contradicted"] = sum(r["region_class"] == "junction" and
                                          r["read_evidence"] == "read-contradicted" for r in rs)
        s["junctions_not_spanned"] = sum(r["region_class"] == "junction" and
                                         r["read_evidence"] == "not-spanned" for r in rs)
        s["edits"] = sum(r["region_class"] == "edit" for r in rs)
        s["consistent_high_confidence_Mb"] = round(sum(r["length_bp"] for r in rs if r.get("confidence") == "high")
                                                   / 1e6, 2)
        summ.append(s)
    rip.write_tsv(os.path.join(partial, "unit_summary.tsv"), summ, list(summ[0]))
    groups = {}
    for u in units:
        groups.setdefault(u["group"], []).append(u)
    figs = [] if a.no_plots else plot(partial, units, seqs_by_unit, regions_by_unit, ev, groups, chroms, P)
    with open(os.path.join(partial, "METHODS.md"), "w") as fh:
        fh.write(METHODS.format(stamp=start, cv=CODE_VERSION, rip=os.path.abspath(a.riparian),
                                redundant_pct=P["redundant_frac"] * 100, **P))
    inputs = {k: (os.path.abspath(p), rip.sha256(p)) for k, p in
              (("riparian_sequences", os.path.join(a.riparian, "sequences.tsv")),
               ("riparian_ribbons", os.path.join(a.riparian, "ribbons.tsv")),
               ("evidence", a.evidence), ("windows", a.windows)) if p}
    with open(os.path.join(partial, "RUN.json"), "w") as fh:
        json.dump(dict(schema="g2f-structure-map/1", code_version=CODE_VERSION, start=start,
                       end=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
                       command=[sys.executable] + sys.argv, cwd=os.getcwd(), script_sha256=rip.sha256(__file__),
                       pipeline_git=rip.git_state(), params=P, inputs=inputs,
                       outputs=[os.path.relpath(p, partial) for p in figs] +
                       ["regions.tsv", "flagged_regions.tsv", "unit_summary.tsv", "METHODS.md"]), fh, indent=1)
    os.rename(partial, a.out)
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
