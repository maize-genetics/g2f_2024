#!/usr/bin/env python3
"""Riparian (GENESPACE-style) synteny overview of assemblies against a reference, from AnchorWave gene anchors.

Reads existing AnchorWave `proali` anchor files (no alignment is run), lifts anchors onto the delivered FASTA records
(optionally through a coordinate map for split/broken releases), aggregates anchors into block-level ribbons, assigns
each record a dominant reference chromosome, and writes:

  sequences.tsv         one row per FASTA record (plotted or not): dominant chromosome, anchors, spans, orientation,
                        candidate inter-chromosomal flag, redundancy flag, plot order
  ribbons.tsv           one row per aggregated anchor block (drawn and undrawn)
  unit_summary.tsv      one row per unit (assembly FASTA): counts per dominant chromosome, flags
  INPUTS.tsv, RUN.json  every input with size and sha256 (anchors, maps, .fai), parameters, Git commit, command
  METHODS.md            how ordering, aggregation and assignment were done, with the parameter values used
  figures/              riparian_<group>.{pdf,svg,png} (one panel per line), riparian_by_line.pdf (one page per line),
                        lines/<line>.{svg,png}

Units table (TSV, one row per assembly FASTA; several rows may share a `line`, e.g. hap1 + hap2):
  line, unit, group, label, anchors, fai          required
  coordinate_map   optional: TSV mapping delivered records to the contigs AnchorWave saw
                   (columns record|candidate_record, original_contig, original_start, original_end; 1-based)
  name_map         optional: TSV whose first two columns map record id -> display name (header line skipped)
  assignment       optional: TSV giving each record's chromosome label in the assembly's own release
                   (columns named by assignment_record_col / assignment_chr_col; defaults record / chromosome)
  any other columns (fasta, sha256, evidence ...) are copied into INPUTS.tsv untouched

Two record orders (--order):
  synteny   (default) GENESPACE-like: dominant reference chromosome, then median anchor position; records whose
            anchors are mostly '-' are drawn reverse-complemented
  assembly  records grouped by the `assignment` table, in FASTA order within a group, in their released
            orientation; a whole group is mirrored only when >= --group-flip-frac of its anchors are '-'. This is a
            chromosome-grouped contig view: the order inside a group is not a physical order
  both      both figure sets, plus figures/compare_<line> (synteny above, assembly below) for --compare lines

Interpretation: a record with ribbons to several chromosomes is a *candidate* for review, not a proven misjoin
(reference differences, duplicated genes and repeats give the same pattern).

Usage:
  g2f_riparian.py check --units units.tsv --ref-fai ref.fa.fai
  g2f_riparian.py run   --units units.tsv --ref-fai ref.fa.fai --out OUTDIR [--order both --compare LINE,LINE]
"""
import copy
import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
import textwrap
from collections import defaultdict

CODE_VERSION = 2
REQUIRED = ("line", "unit", "group", "label", "anchors", "fai")
# 10 reference chromosomes: the 8 validated categorical slots, then blue-600 and blue-300 (adjacent pairs pass the
# CVD / normal-vision checks; identity is also given by position on the reference track and by text labels).
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948", "#184f95",
           "#6da7ec"]
FLAG = "#c0262d"
INK, INK2, MUTED, BAR = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d3"


def die(msg):
    sys.exit(f"g2f_riparian: error: {msg}")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def read_tsv(path):
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def read_fai(path):
    out = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            out[f[0]] = int(f[1])
    return out


def read_units(path):
    units = read_tsv(path)
    if not units:
        die(f"{path}: no units")
    for u in units:
        miss = [k for k in REQUIRED if not u.get(k)]
        if miss:
            die(f"{path}: unit {u.get('unit', '?')} lacks {', '.join(miss)}")
        for k in ("anchors", "fai", "coordinate_map", "name_map", "assignment"):
            if u.get(k) and not os.path.isfile(u[k]):
                die(f"{u['unit']}: {k} not found: {u[k]}")
    if len({u["unit"] for u in units}) != len(units):
        die("unit ids must be unique")
    return units


def ref_chroms(fai_path, pattern):
    rx = re.compile(pattern)
    fai = read_fai(fai_path)
    chroms = sorted((c for c in fai if rx.fullmatch(c)), key=natural)
    if not chroms:
        die(f"no reference sequence matches {pattern!r} in {fai_path}")
    return chroms, {c: fai[c] for c in chroms}


def read_anchors(path, chroms):
    """Gene anchors (interanchor rows skipped) as dicts; returns (anchors, n_other_ref)."""
    out, other = [], 0
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("refChr"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[7] == "interanchor":
                continue
            if f[0] not in chroms:
                other += 1
                continue
            out.append(dict(rc=f[0], rs=int(f[1]), re=int(f[2]), q=f[3], qs=int(f[4]), qe=int(f[5]), strand=f[6],
                            gene=f[7], block=f[8]))
    return out, other


def read_coordinate_map(path):
    pieces = defaultdict(list)  # original contig -> [(start, end, record)]
    for r in read_tsv(path):
        rec = r.get("record") or r.get("candidate_record")
        if not rec:
            die(f"{path}: needs a record or candidate_record column")
        pieces[r["original_contig"]].append((int(r["original_start"]), int(r["original_end"]), rec))
    for c, p in pieces.items():
        p.sort()
        pos = 1
        for s, e, rec in p:
            if s != pos:
                die(f"{path}: {c} pieces are not contiguous at {rec} ({s} != {pos})")
            pos = e + 1
    return pieces


def lift(anchors, fai, cmap):
    """Move anchors onto delivered records. A piece is chosen by the anchor midpoint; anchors that straddle a split
    point are counted and kept on the piece holding their midpoint."""
    lifted, straddle, missing = [], 0, set()
    for a in anchors:
        if cmap and a["q"] in cmap:
            mid = (a["qs"] + a["qe"]) // 2
            for s, e, rec in cmap[a["q"]]:
                if s <= mid <= e:
                    break
            else:
                die(f"anchor at {a['q']}:{mid} outside the coordinate map")
            if a["qs"] < s or a["qe"] > e:
                straddle += 1
            a = dict(a, q=rec, qs=max(a["qs"], s) - s + 1, qe=min(a["qe"], e) - s + 1)
        if a["q"] not in fai:
            missing.add(a["q"])
            continue
        if a["qe"] > fai[a["q"]]:
            die(f"anchor {a['gene']} at {a['q']}:{a['qe']} beyond record length {fai[a['q']]}")
        lifted.append(a)
    if missing:
        die(f"{len(missing)} anchored contigs are not in the FASTA index (e.g. {sorted(missing)[:3]}); "
            "a coordinate map is missing or the anchors belong to another assembly")
    return lifted, straddle


def segment(anchors, max_gap, min_n):
    """Aggregate one record's anchors (any order) into ribbons. Consecutive anchors along the record join a ribbon
    while reference chromosome and strand agree and both reference and query gaps stay <= max_gap. Ribbons with fewer
    than min_n anchors are then set aside and neighbouring compatible ribbons are merged across them, so a single
    out-of-place gene does not cut a long block in two."""
    def runs(items):
        out = []
        for a in items:
            cur = out[-1] if out else None
            if (cur and cur["rc"] == a["rc"] and cur["strand"] == a["strand"]
                    and a["qs"] - cur["qe"] <= max_gap
                    and min(abs(a["rs"] - cur["last_rs"]), abs(a["re"] - cur["last_rs"])) <= max_gap):
                cur["items"].append(a)
                cur["qe"] = max(cur["qe"], a["qe"])
                cur["last_rs"] = a["rs"]
            else:
                out.append(dict(rc=a["rc"], strand=a["strand"], qs=a["qs"], qe=a["qe"], last_rs=a["rs"], items=[a]))
        return out
    items = sorted(anchors, key=lambda a: (a["qs"], a["qe"]))
    first = runs(items)
    big = [a for r in first if len(r["items"]) >= min_n for a in r["items"]]
    small = [a for r in first if len(r["items"]) < min_n for a in r["items"]]
    merged = runs(big)
    for r in merged:
        r["drawn"] = len(r["items"]) >= min_n
    leftover = runs(small)
    for r in leftover:
        r["drawn"] = False
    ribbons = sorted(merged + leftover, key=lambda r: r["qs"])
    for r in ribbons:
        r["n"] = len(r["items"])
        r["rs"] = min(a["rs"] for a in r["items"])
        r["re"] = max(a["re"] for a in r["items"])
        r["qs"] = min(a["qs"] for a in r["items"])
        r["qe"] = max(a["qe"] for a in r["items"])
    return ribbons


def union_len(iv):
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted(iv):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                tot += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return tot + (cur_e - cur_s + 1 if cur_e is not None else 0)


def mb(x):
    return f"{x / 1e6:.2f}"


def summarise(unit, fai, anchors, names, chroms, P):
    by_rec = defaultdict(list)
    for a in anchors:
        by_rec[a["q"]].append(a)
    rows, ribbons = [], []
    for rec in sorted(fai, key=natural):
        A = by_rec.get(rec, [])
        rb = segment(A, P["max_gap"], P["min_ribbon_anchors"]) if A else []
        for r in rb:
            r["record"] = rec
        drawn = [r for r in rb if r["drawn"]]
        per = defaultdict(lambda: dict(n=0, plus=0, rib_q=[], rib_r=[], drawn_n=0))
        for a in A:
            per[a["rc"]]["n"] += 1
            per[a["rc"]]["plus"] += a["strand"] == "+"
        for r in drawn:
            per[r["rc"]]["rib_q"].append((r["qs"], r["qe"]))
            per[r["rc"]]["rib_r"].append((r["rs"], r["re"]))
            per[r["rc"]]["drawn_n"] += r["n"]
        row = dict(line=unit["line"], unit=unit["unit"], record=rec, display_name=names.get(rec, rec),
                   length_bp=fai[rec], gene_anchors=len(A), anchors_in_drawn_ribbons=sum(r["n"] for r in drawn),
                   ribbons_drawn=len(drawn), plotted="yes" if drawn else "no")
        if not A:
            row.update(dominant_chr="", status="no anchors")
            rows.append(row)
            continue
        # dominant chromosome: most gene anchors; tie -> larger drawn query span
        dom = max(per, key=lambda c: (per[c]["n"], union_len(per[c]["rib_q"]), -chroms.index(c)))
        d = per[dom]
        frac_plus = d["plus"] / d["n"]
        orient = "forward" if frac_plus >= 0.9 else ("reverse" if frac_plus <= 0.1 else "mixed")
        major = "+" if frac_plus >= 0.5 else "-"
        minority = [r for r in drawn if r["rc"] == dom and r["strand"] != major]
        dom_anchor_pos = sorted((a["rs"] + a["re"]) // 2 for a in A if a["rc"] == dom)
        substantial = [c for c in per if per[c]["n"] >= P["sec_min_anchors"]
                       and union_len(per[c]["rib_q"]) >= P["sec_min_span"]]
        # drawn blocks on another chromosome: terminal ones (near a record end) are the signature of a join
        off = [r for r in drawn if r["rc"] != dom]
        for r in off:
            r["terminal"] = r["qs"] <= P["terminal_bp"] or r["qe"] >= fai[rec] - P["terminal_bp"]
        off_txt = "; ".join(f"{r['rc']} {r['strand']} {r['n']} anchors, record {mb(r['qs'])}-{mb(r['qe'])} Mb "
                            f"({'terminal' if r['terminal'] else 'internal'}) -> {r['rc']}:{mb(r['rs'])}-{mb(r['re'])} Mb"
                            for r in off)
        tier = "yes" if len(substantial) >= 2 else ("minor" if off else "no")
        secondary = sorted((c for c in per if c != dom and per[c]["n"] > 0), key=lambda c: -per[c]["n"])
        sec_txt = "; ".join(
            f"{c}: {per[c]['n']} anchors" + (
                f", record {mb(min(s for s, _ in per[c]['rib_q']))}-{mb(max(e for _, e in per[c]['rib_q']))} Mb"
                f" -> {c}:{mb(min(s for s, _ in per[c]['rib_r']))}-{mb(max(e for _, e in per[c]['rib_r']))} Mb"
                if per[c]["rib_q"] else " (no ribbon >= threshold)") for c in secondary)
        row.update(
            dominant_chr=dom,
            anchors_on_dominant=d["n"],
            frac_anchors_on_dominant=round(d["n"] / len(A), 4),
            dominant_ref_start_bp=min(s for s, _ in d["rib_r"]) if d["rib_r"] else "",
            dominant_ref_end_bp=max(e for _, e in d["rib_r"]) if d["rib_r"] else "",
            dominant_ref_span_bp=(max(e for _, e in d["rib_r"]) - min(s for s, _ in d["rib_r"]) + 1) if d["rib_r"] else 0,
            dominant_ref_covered_bp=union_len(d["rib_r"]),
            dominant_query_covered_bp=union_len(d["rib_q"]),
            frac_record_in_dominant_ribbons=round(union_len(d["rib_q"]) / fai[rec], 4),
            median_anchor_pos_on_dominant=dom_anchor_pos[len(dom_anchor_pos) // 2],
            n_chr_touched=len(per),
            n_chr_substantial=len(substantial),
            chr_substantial=",".join(sorted(substantial, key=natural)),
            candidate_interchromosomal=tier,
            off_chr_blocks=len(off),
            off_chr_terminal_blocks=sum(r["terminal"] for r in off),
            off_chr_block_detail=off_txt,
            orientation_vs_ref=orient,
            frac_plus_on_dominant=round(frac_plus, 4),
            display_reversed="yes" if major == "-" else "no",
            minority_strand_ribbons=len(minority),
            minority_strand_query_bp=union_len([(r["qs"], r["qe"]) for r in minority]),
            other_chromosomes=sec_txt,
            status="",
        )
        rows.append(row)
        ribbons.extend(rb)
    plotted = [r for r in rows if r["plotted"] == "yes"]
    # redundancy: >= redundant_frac of a record's dominant reference interval also covered by a record with more
    # anchors on that chromosome (possible extra allele / duplicated copy; flag only)
    cover = defaultdict(list)
    for r in ribbons:
        if r["drawn"]:
            cover[(r["record"], r["rc"])].append((r["rs"], r["re"]))
    for row in plotted:
        dom, mine = row["dominant_chr"], cover.get((row["record"], row["dominant_chr"]), [])
        own = union_len(mine)
        best, best_frac = "", 0.0
        for other in plotted:
            if other is row or other["dominant_chr"] != dom or other["anchors_on_dominant"] <= row["anchors_on_dominant"]:
                continue
            ov = own + union_len(cover[(other["record"], dom)]) - union_len(mine + cover[(other["record"], dom)])
            if own and ov / own > best_frac:
                best, best_frac = other["record"], ov / own
        row["ref_overlap_with_larger"] = f"{best} ({best_frac * 100:.0f}%)" if best else ""
        row["redundant_candidate"] = "yes" if best_frac >= P["redundant_frac"] else "no"
    # status label
    for row in rows:
        if row["status"]:
            continue
        if row["plotted"] == "no":
            row["status"] = "anchors below ribbon threshold"
        elif row["candidate_interchromosomal"] == "yes":
            row["status"] = "candidate inter-chromosomal"
        elif row["candidate_interchromosomal"] == "minor":
            row["status"] = ("minor block on another chromosome (record end)" if row["off_chr_terminal_blocks"]
                             else "minor block on another chromosome (internal)")
        elif row["redundant_candidate"] == "yes":
            row["status"] = "overlaps a larger record (redundancy candidate)"
        elif row["minority_strand_ribbons"]:
            row["status"] = "one chromosome; contains inverted block(s)"
        else:
            row["status"] = "one chromosome, collinear"
    # plot order: dominant chromosome, then median anchor position on it, then longer first
    plotted.sort(key=lambda r: (chroms.index(r["dominant_chr"]), r["median_anchor_pos_on_dominant"], -r["length_bp"]))
    for i, row in enumerate(plotted, 1):
        row["plot_order"] = i
    return rows, ribbons


def read_assignment(unit, chroms):
    """record -> chromosome label from the unit's assignment table; labels that are not reference chromosomes (e.g.
    'unplaced') map to ''."""
    rc, cc = unit.get("assignment_record_col") or "record", unit.get("assignment_chr_col") or "chromosome"
    rows = read_tsv(unit["assignment"])
    if rows and (rc not in rows[0] or cc not in rows[0]):
        die(f"{unit['unit']}: assignment table lacks column {rc!r} or {cc!r}")
    out = {}
    for r in rows:
        if r[rc] in out and out[r[rc]] != r[cc]:
            die(f"{unit['unit']}: record {r[rc]} assigned twice")
        out[r[rc]] = r[cc] if r[cc] in chroms else ""
    return out


def assembly_order(rows, fai_order, assign, chroms, flip_frac):
    """Assembly-native layout: group records by their assigned chromosome (unassigned last), FASTA order inside a
    group, released orientation. A group is mirrored as a whole (order and every record reversed) only when
    >= flip_frac of the anchors that records assigned to it have on that chromosome are on '-'."""
    groups = defaultdict(list)
    for row in rows:
        a = assign.get(row["record"])
        row["assigned_chr"] = a if a is not None else ""
        row["assignment_agrees"] = ("" if not a or not row["dominant_chr"] else
                                    "yes" if a == row["dominant_chr"] else "no")
        if a or row["plotted"] == "yes":
            groups[a or "unassigned"].append(row)
    order = [c for c in chroms if c in groups] + (["unassigned"] if "unassigned" in groups else [])
    i = 0
    for g in order:
        members = sorted(groups[g], key=lambda r: fai_order[r["record"]])
        n = sum(r["anchors_on_dominant"] for r in members if r["dominant_chr"] == g)
        minus = sum(r["anchors_on_dominant"] * (1 - r["frac_plus_on_dominant"]) for r in members
                    if r["dominant_chr"] == g)
        rev = g != "unassigned" and n > 0 and minus / n >= flip_frac
        for r in (members[::-1] if rev else members):
            i += 1
            r["assembly_order"] = i
            r["assembly_group"] = g
            r["group_frac_minus"] = round(minus / n, 4) if n else ""
            r["group_reversed"] = "yes" if rev else "no"


def view(data, mode):
    """Per-figure copy of the record rows with the fields the drawing code reads set for one order."""
    out = {}
    for unit, d in data.items():
        rows = []
        for r in d["rows"]:
            r = dict(r)
            if mode == "synteny":
                r["shown"] = r["plotted"]
                r["_group"] = r["dominant_chr"]
            else:
                r["shown"] = "yes" if r.get("assembly_order") else "no"
                r["plot_order"] = r.get("assembly_order")
                r["_group"] = r.get("assembly_group", "")
                r["display_reversed"] = r.get("group_reversed", "no")
            rows.append(r)
        out[unit] = dict(d, rows=rows)
    return out


COLS = ["line", "unit", "plot_order", "record", "display_name", "length_bp", "plotted", "status", "dominant_chr",
        "gene_anchors", "anchors_on_dominant", "frac_anchors_on_dominant", "anchors_in_drawn_ribbons", "ribbons_drawn",
        "dominant_ref_start_bp", "dominant_ref_end_bp", "dominant_ref_span_bp", "dominant_ref_covered_bp",
        "dominant_query_covered_bp", "frac_record_in_dominant_ribbons", "median_anchor_pos_on_dominant",
        "n_chr_touched", "n_chr_substantial", "chr_substantial", "candidate_interchromosomal", "off_chr_blocks",
        "off_chr_terminal_blocks", "off_chr_block_detail", "orientation_vs_ref",
        "frac_plus_on_dominant", "display_reversed", "minority_strand_ribbons", "minority_strand_query_bp", "redundant_candidate",
        "ref_overlap_with_larger", "other_chromosomes"]
ASM_COLS = ["assigned_chr", "assignment_agrees", "assembly_order", "assembly_group", "group_reversed",
            "group_frac_minus"]


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


def unit_summary(unit, rows, chroms, straddle, other_ref):
    pl = [r for r in rows if r["plotted"] == "yes"]
    out = dict(line=unit["line"], unit=unit["unit"], label=unit["label"], records_total=len(rows),
               bp_total=sum(r["length_bp"] for r in rows), records_plotted=len(pl),
               bp_plotted=sum(r["length_bp"] for r in pl),
               records_anchors_below_threshold=sum(r["status"] == "anchors below ribbon threshold" for r in rows),
               records_no_anchors=sum(r["gene_anchors"] == 0 for r in rows),
               bp_no_anchors=sum(r["length_bp"] for r in rows if r["gene_anchors"] == 0),
               gene_anchors=sum(r["gene_anchors"] for r in rows),
               candidate_interchromosomal=sum(r.get("candidate_interchromosomal") == "yes" for r in pl),
               candidate_interchromosomal_records=",".join(r["record"] for r in pl
                                                           if r["candidate_interchromosomal"] == "yes"),
               minor_interchromosomal=sum(r.get("candidate_interchromosomal") == "minor" for r in pl),
               minor_interchromosomal_terminal=sum(r.get("candidate_interchromosomal") == "minor"
                                                   and r["off_chr_terminal_blocks"] > 0 for r in pl),
               minor_interchromosomal_records=",".join(r["record"] for r in pl
                                                       if r["candidate_interchromosomal"] == "minor"),
               records_touching_gt1_chr=sum(r.get("n_chr_touched", 0) > 1 for r in pl),
               records_with_inverted_blocks=sum(r.get("minority_strand_ribbons", 0) > 0 for r in pl),
               redundant_candidates=sum(r.get("redundant_candidate") == "yes" for r in pl),
               anchors_straddling_split=straddle, anchors_on_other_ref_seqs=other_ref)
    for c in chroms:
        out[f"n_{c}"] = sum(r["dominant_chr"] == c for r in pl)
    for c in chroms:
        out[f"bp_{c}"] = sum(r["length_bp"] for r in pl if r["dominant_chr"] == c)
    return out


# ------------------------------------------------------------------------------------------------------------ plots
def band(ax, x1a, x1b, y1, x2a, x2b, y2, **kw):
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    ym = (y1 + y2) / 2
    verts = [(x1a, y1), (x1a, ym), (x2a, ym), (x2a, y2), (x2b, y2), (x2b, ym), (x1b, ym), (x1b, y1), (x1a, y1)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), **kw))


def layout_query(rows, gap, group_gap):
    pos, x, prev = {}, 0, None
    for r in sorted(rows, key=lambda r: r["plot_order"]):
        if prev is not None:
            x += group_gap if r["_group"] != prev else gap
        pos[r["record"]] = x
        x += r["length_bp"]
        prev = r["_group"]
    return pos, x


def draw_line(ax, line_units, data, chroms, clen, xmax, P, detail, mode="synteny"):
    """One panel: reference track in the middle (two units: unit 1 above, unit 2 below) or on top (one unit)."""
    from matplotlib.patches import FancyBboxPatch, Rectangle
    col = dict(zip(chroms, PALETTE))
    H = 0.10  # bar height
    ref_gap = P["ref_gap_frac"] * sum(clen.values())
    ref_w = sum(clen.values()) + ref_gap * (len(chroms) - 1)
    ref_x0 = (xmax - ref_w) / 2
    rpos, x = {}, ref_x0
    for c in chroms:
        rpos[c] = x
        x += clen[c] + ref_gap
    if len(line_units) == 1:
        tracks = [(line_units[0], 0.0, -1)]
        yref = 1.0
    else:
        tracks = [(line_units[0], 2.0, +1), (line_units[1], 0.0, -1)]
        yref = 1.0
    xpad = xmax * 0.004
    fs = 7 if not detail else 8
    # reference bars and labels
    for c in chroms:
        ax.add_patch(FancyBboxPatch((rpos[c], yref - H / 2), clen[c], H, boxstyle=f"round,pad=0,rounding_size={H / 3}",
                                    mutation_aspect=1 / (xmax / 12), fc=col[c], ec="none", zorder=3))
        # label on the side without ribbons (above for one unit; inside the gap above the bar for two)
        ax.text(rpos[c] + clen[c] / 2, yref + H / 2 + 0.03, c, ha="center", va="bottom", fontsize=fs, color=INK,
                fontweight="bold", zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.8) if len(line_units) == 2 else None)
    ax.text(ref_x0 - xpad * 2, yref, "B73 v5", ha="right", va="center", fontsize=fs + 1, color=INK)
    for unit, yq, side in tracks:
        rows = [r for r in data[unit["unit"]]["rows"] if r["shown"] == "yes"]
        rib = data[unit["unit"]]["ribbons"]
        qgap = P["query_gap_frac"] * xmax
        qpos, qw = layout_query(rows, qgap, ref_gap)
        q0 = (xmax - qw) / 2
        info = {r["record"]: r for r in rows}
        # ribbons: from the reference bar edge facing the query to the query bar edge facing the reference
        y_r = yref + H / 2 if side > 0 else yref - H / 2
        y_q = yq - H / 2 if side > 0 else yq + H / 2
        # records whose dominant-chromosome anchors are mostly on '-' are drawn reverse-complemented (as GENESPACE
        # does), so a twisted ribbon means opposite orientation to the rest of the record, not hifiasm's arbitrary
        # contig orientation
        def qx(row, a, b):
            if row["display_reversed"] == "yes":
                a, b = row["length_bp"] - b + 1, row["length_bp"] - a + 1
            return q0 + qpos[row["record"]] + a, q0 + qpos[row["record"]] + b

        for r in rib:
            if not r["drawn"] or r["record"] not in info:
                continue
            row = info[r["record"]]
            # off = not the chromosome the record is grouped under (its dominant one, or its released assignment)
            off = r["rc"] != (row["_group"] if row["_group"] in col else row["dominant_chr"])
            xa, xb = rpos[r["rc"]] + r["rs"], rpos[r["rc"]] + r["re"]
            qa, qb = qx(row, r["qs"], r["qe"])
            if (r["strand"] == "-") != (row["display_reversed"] == "yes"):
                qa, qb = qb, qa
            band(ax, xa, xb, y_r, qa, qb, y_q, fc=col[r["rc"]], alpha=0.85 if off else 0.45,
                 ec=INK if off else col[r["rc"]], lw=0.6 if off else 0.15, zorder=6 if off else 2)
        # query bars
        groups = defaultdict(list)
        for row in rows:
            x0 = q0 + qpos[row["record"]]
            tier = row.get("candidate_interchromosomal") or "no"
            flagged = tier != "no"
            ax.add_patch(Rectangle((x0, yq - H / 2), row["length_bp"], H, fc=BAR, ec=FLAG if flagged else "none",
                                   lw=1.2 if flagged else 0, ls=(0, (2, 1.2)) if tier == "minor" else "-",
                                   zorder=3 + flagged))
            groups[row["_group"]].append((x0, x0 + row["length_bp"]))
            if not row["dominant_chr"]:
                continue
            # minority-strand (inverted relative to the rest of the record) blocks, on the side away from the ribbons
            major = "+" if row["frac_plus_on_dominant"] >= 0.5 else "-"
            ys = yq + side * (H / 2 + 0.035)
            for r in rib:
                if r["record"] == row["record"] and r["drawn"] and r["rc"] == row["dominant_chr"] and r["strand"] != major:
                    a, b = qx(row, r["qs"], r["qe"])
                    ax.add_patch(Rectangle((a, ys - 0.025), max(b - a, xmax * 1.5e-3), 0.05,
                                           fc=INK, ec="none", zorder=5))
            lab = None
            if flagged:
                lab = row["display_name"]
            elif detail and row["length_bp"] >= P["label_min_bp"]:
                lab = row["display_name"]
            if lab:
                ax.text(x0 + row["length_bp"] / 2, yq + side * (H / 2 + 0.09), lab, rotation=90,
                        ha="center", va="bottom" if side > 0 else "top", fontsize=fs - 2 if not flagged else fs - 1,
                        color=FLAG if flagged else INK2, zorder=7, fontweight="bold" if flagged else "normal")
        # group brackets on the ribbon side of the bar (dominant chromosome, or released assignment)
        for c, iv in groups.items():
            a, b = min(s for s, _ in iv), max(e for _, e in iv)
            yb = yq - side * (H / 2 + 0.03)
            ax.plot([a, b], [yb, yb], color=col.get(c, MUTED), lw=2.2, solid_capstyle="butt", zorder=4)
            if mode == "assembly":  # the group label, since the order no longer follows the reference
                ax.text((a + b) / 2, yb - side * 0.02, c, ha="center", va="top" if side > 0 else "bottom",
                        fontsize=fs - 1, color=INK, zorder=8, style="italic",
                        bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
        allrows = data[unit["unit"]]["rows"]
        nskip = sum(r["shown"] == "no" for r in allrows)
        bskip = sum(r["length_bp"] for r in allrows if r["shown"] == "no")
        nflag = sum(r.get("candidate_interchromosomal") == "yes" for r in rows)
        nminor = sum(r.get("candidate_interchromosomal") == "minor" for r in rows)
        if mode == "assembly":
            nun = sum(r["_group"] == "unassigned" for r in rows)
            ndis = sum(r.get("assignment_agrees") == "no" for r in rows)
            txt = (f"{len(rows)} records shown ({nun} unassigned)\n{nskip} not shown ({bskip / 1e6:.1f} Mb)\n"
                   f"assignment != dominant: {ndis}")
        else:
            txt = (f"{len(rows)} records plotted\n{nskip} not plotted ({bskip / 1e6:.1f} Mb)\n"
                   f"inter-chr: {nflag} candidate, {nminor} minor")
        ax.text(q0 - xpad * 2, yq, unit["track_label"], ha="right", va="center", fontsize=fs + 1, color=INK)
        ax.text(q0 + qw + xpad * 2, yq, txt, ha="left", va="center", fontsize=fs - 1, color=INK2, linespacing=1.15)
    ylo = -1.0 if detail else -0.85
    yhi = (3.0 if detail else 2.85) if len(line_units) == 2 else 1.35
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(-xmax * 0.06, xmax * 1.09)
    ax.axis("off")


def legend(fig, chroms, y, mode="synteny"):
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.lines import Line2D
    col = dict(zip(chroms, PALETTE))
    h = [Patch(fc=col[c], ec="none", alpha=0.75, label=c) for c in chroms]
    h += [Patch(fc=BAR, ec=FLAG, lw=1.2, label="candidate inter-chromosomal record"),
          Patch(fc=BAR, ec=FLAG, lw=1.2, ls=(0, (2, 1.2)), label="minor block on another chromosome"),
          Patch(fc="#999999", ec=INK, lw=0.6, label="ribbon to a non-dominant chromosome"),
          Rectangle((0, 0), 1, 1, fc=INK, ec="none", label="block inverted vs. rest of record"),
          Line2D([], [], color=MUTED, lw=2.2, label="bracket = dominant-chromosome group" if mode == "synteny"
                 else "bracket = assigned-chromosome group")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, y), ncol=8, fontsize=8, frameon=False,
               handlelength=1.4, columnspacing=1.2)


SUBTITLE = {
    "synteny": ("Ribbons: AnchorWave gene-anchor blocks (>= {min_ribbon_anchors} anchors) coloured by B73 chromosome. "
                "Records drawn in their majority orientation, so a twist = inverted relative to the rest of the record. "
                "Records ordered by dominant B73 chromosome, then median anchor position. "
                "Visual QC only: multi-chromosome links are candidates, not confirmed errors."),
    "assembly": ("CHROMOSOME-GROUPED CONTIGS, NOT SCAFFOLDS: records grouped by the assembly's own chromosome "
                 "assignment, in FASTA order inside a group; order and spacing inside a group carry no physical "
                 "meaning. Released orientation kept (a group is mirrored as a whole only if >= {group_flip_pct:.0f}% "
                 "of its anchors are reverse), so a twisted ribbon = that segment runs opposite to B73 as drawn. "
                 "Ribbons coloured by B73 chromosome; outlined ribbons leave the record's assigned chromosome."),
}
TITLE = {"synteny": "Riparian synteny vs B73 v5", "assembly": "Chromosome-grouped assembly contigs vs B73 v5"}
PREFIX = {"synteny": "riparian", "assembly": "asmorder"}


def plot_all(out, units, data, chroms, clen, P, modes=("synteny",), compare=()):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    figdir = os.path.join(out, "figures")
    lines, groups = {}, {}
    for u in units:
        lines.setdefault(u["line"], []).append(u)
        groups.setdefault(u["group"], [])
        if u["line"] not in groups[u["group"]]:
            groups[u["group"]].append(u["line"])
    views = {m: view(data, m) for m in modes}
    # one common bp scale for every panel, figure and order
    ref_w = sum(clen.values()) * (1 + P["ref_gap_frac"] * (len(chroms) - 1) / len(chroms))
    xmax = ref_w
    for v in views.values():
        for u in units:
            rows = [r for r in v[u["unit"]]["rows"] if r["shown"] == "yes"]
            _, w = layout_query(rows, 0, P["ref_gap_frac"] * sum(clen.values()))
            xmax = max(xmax, w / (1 - P["query_gap_frac"] * len(rows)))
    written = []

    def panel_h(line):
        return 2.9 if len(lines[line]) == 2 else 1.8

    for mode in modes:
        d, pre = views[mode], PREFIX[mode]
        sub = SUBTITLE[mode].format(**P)
        os.makedirs(os.path.join(figdir, "lines" if mode == "synteny" else "lines_asmorder"))
        for g, glines in groups.items():
            hs = [panel_h(l) for l in glines]
            fig = plt.figure(figsize=(20, sum(hs) + 1.3))
            top = 1 - (1.15 if mode == "synteny" else 1.45) / (sum(hs) + 1.3)
            gs = fig.add_gridspec(len(glines), 1, height_ratios=hs, top=top, bottom=0.01, left=0.01, right=0.99,
                                  hspace=0.12)
            for i, line in enumerate(glines):
                ax = fig.add_subplot(gs[i])
                draw_line(ax, lines[line], d, chroms, clen, xmax, P, detail=False, mode=mode)
                ax.set_title(f"{line} — {lines[line][0]['label']}", loc="left", fontsize=10, color=INK, pad=2)
            fig.suptitle(f"{TITLE[mode]} — {g}", x=0.01, ha="left", y=1 - 0.15 / (sum(hs) + 1.3),
                         fontsize=13, color=INK)
            if mode == "synteny":  # unchanged from code version 1
                fig.text(0.01, 1 - 0.45 / (sum(hs) + 1.3), sub, fontsize=8, color=INK2)
                legend(fig, chroms, 1 - 0.55 / (sum(hs) + 1.3))
            else:
                fig.text(0.01, 1 - 0.45 / (sum(hs) + 1.3), textwrap.fill(sub, 260), fontsize=8, color=INK2,
                         va="center")
                legend(fig, chroms, 1 - 0.66 / (sum(hs) + 1.3), mode)
            for ext in ("pdf", "svg", "png"):
                p = os.path.join(figdir, f"{pre}_{g}.{ext}")
                fig.savefig(p, dpi=150 if ext == "png" else None)
                written.append(p)
            plt.close(fig)
        pdf_path = os.path.join(figdir, f"{pre}_by_line.pdf")
        with PdfPages(pdf_path) as pdf:
            for g, glines in groups.items():
                for line in glines:
                    h = panel_h(line) * 1.6 + 1.3
                    fig = plt.figure(figsize=(20, h))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 1 - 1.1 / h])
                    draw_line(ax, lines[line], d, chroms, clen, xmax, P, detail=True, mode=mode)
                    fig.suptitle(f"{line} vs B73 v5 — {lines[line][0]['label']}", x=0.01, ha="left",
                                 y=1 - 0.15 / h, fontsize=13, color=INK)
                    fig.text(0.01, 1 - 0.3 / h, textwrap.fill(sub + f" Labels: records >= "
                             f"{P['label_min_bp'] / 1e6:g} Mb and all flagged records.", 240), fontsize=8,
                             color=INK2, va="top")
                    legend(fig, chroms, 1 - 0.55 / h, mode)
                    pdf.savefig(fig)
                    for ext in ("svg", "png"):
                        p = os.path.join(figdir, "lines" if mode == "synteny" else "lines_asmorder", f"{line}.{ext}")
                        fig.savefig(p, dpi=150 if ext == "png" else None)
                        written.append(p)
                    plt.close(fig)
        written.append(pdf_path)
    # side-by-side comparison of the two orders for selected lines
    if compare and len(modes) == 2:
        os.makedirs(os.path.join(figdir, "compare"))
        for line in compare:
            one = len(lines[line]) == 1  # labels hang below a single track: leave room for them
            hs = [panel_h(line) * (1.3 if one else 1.0)] * 2
            H = sum(hs) + 2.0 + (1.6 if one else 0.4)
            fig = plt.figure(figsize=(20, H))
            gs = fig.add_gridspec(2, 1, height_ratios=hs, top=1 - 1.75 / H, bottom=(1.0 if one else 0.3) / H,
                                  left=0.01, right=0.99, hspace=0.75 if one else 0.2)
            for i, mode in enumerate(modes):
                ax = fig.add_subplot(gs[i])
                draw_line(ax, lines[line], views[mode], chroms, clen, xmax, P, detail=True, mode=mode)
                ax.set_title(f"{'A' if i == 0 else 'B'}. {TITLE[mode]}: " +
                             ("records sorted and flipped by B73 synteny" if mode == "synteny" else
                              "records in assembly order and released orientation"),
                             loc="left", fontsize=10, color=INK, pad=2)
            fig.suptitle(f"{line} — the same anchors in two orders ({lines[line][0]['label']})", x=0.01,
                         ha="left", y=1 - 0.15 / H, fontsize=13, color=INK)
            fig.text(0.01, 1 - 0.42 / H, textwrap.fill(
                "A sorts and flips every record by its B73 anchors, which makes any assembly look collinear by "
                "construction. B: " + SUBTITLE["assembly"].format(**P), 260), fontsize=8, color=INK2, va="top")
            legend(fig, chroms, 1 - 0.95 / H)
            for ext in ("pdf", "svg", "png"):
                p = os.path.join(figdir, "compare", f"compare_{line}.{ext}")
                fig.savefig(p, dpi=150 if ext == "png" else None)
                written.append(p)
            plt.close(fig)
    return written


# ------------------------------------------------------------------------------------------------------------- main
def git_state():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", here, "status", "--porcelain", "--", os.path.abspath(__file__)],
                               capture_output=True, text=True).stdout.strip()
        return dict(commit=commit, script_modified_or_untracked=bool(dirty))
    except OSError:
        return dict(commit=None, script_modified_or_untracked=None)


def load(units, chroms, P, quiet=False):
    data = {}
    for u in units:
        fai = read_fai(u["fai"])
        cmap = read_coordinate_map(u["coordinate_map"]) if u.get("coordinate_map") else None
        names = {}
        if u.get("name_map"):
            with open(u["name_map"]) as fh:
                next(fh)
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) >= 2:
                        names[f[0]] = f[1]
        raw, other = read_anchors(u["anchors"], set(chroms))
        if not raw:
            die(f"{u['unit']}: no gene anchors on the selected reference chromosomes")
        anchors, straddle = lift(raw, fai, cmap)
        if cmap:  # every delivered record from a mapped contig must exist in the index with the mapped length
            for c, pieces in cmap.items():
                for s, e, rec in pieces:
                    if fai.get(rec) != e - s + 1:
                        die(f"{u['unit']}: {rec} length {fai.get(rec)} != coordinate map {e - s + 1}")
        rows, ribbons = summarise(u, fai, anchors, names, chroms, P)
        if u.get("assignment"):
            assign = read_assignment(u, set(chroms))
            unknown = sorted(set(assign) - set(fai))
            if unknown:
                die(f"{u['unit']}: {len(unknown)} assigned records are not in the FASTA index (e.g. {unknown[:3]})")
            assembly_order(rows, {k: i for i, k in enumerate(fai)}, assign, chroms, P["group_flip_frac"])
        elif P["order"] != "synteny":
            die(f"{u['unit']}: --order {P['order']} needs an assignment table for every unit")
        data[u["unit"]] = dict(rows=rows, ribbons=ribbons, straddle=straddle, other_ref=other)
        if not quiet:
            pl = sum(r["plotted"] == "yes" for r in rows)
            fl = sum(r.get("candidate_interchromosomal") == "yes" for r in rows)
            mi = sum(r.get("candidate_interchromosomal") == "minor" for r in rows)
            print(f"  {u['unit']}: {len(raw)} gene anchors, {len(rows)} records, {pl} plotted, inter-chromosomal "
                  f"{fl} candidate / {mi} minor, {straddle} anchors straddle a split", file=sys.stderr)
    return data


def params(a):
    return dict(min_ribbon_anchors=a.min_ribbon_anchors, max_gap=a.max_gap_bp, sec_min_anchors=a.sec_min_anchors,
                sec_min_span=a.sec_min_span_bp, terminal_bp=a.terminal_bp, redundant_frac=a.redundant_frac, label_min_bp=a.label_min_bp,
                ref_gap_frac=0.006, query_gap_frac=0.0012, ref_pattern=a.ref_pattern, order=a.order,
                group_flip_frac=a.group_flip_frac, group_flip_pct=a.group_flip_frac * 100)


METHODS = """# Riparian synteny overview: methods ({stamp})

Generated by `g2f_riparian.py` (code version {cv}). Parameters are those recorded in RUN.json.

**Input.** Existing AnchorWave `proali` anchor files (one per unit; see INPUTS.tsv). Only gene-anchor rows are used
(`interanchor` rows are skipped). Anchors to reference sequences other than `{ref_pattern}` are counted
(`anchors_on_other_ref_seqs`) and not plotted. No alignment was run.

**Lifting to delivered records.** Where a unit has a coordinate map (records split or broken from the contigs that
AnchorWave saw), each anchor is moved to the piece that holds its midpoint, with coordinates re-based to the piece.
The map must tile each original contig without gaps, and each piece length must equal the FASTA index length.
Anchors that cross a split point are counted (`anchors_straddling_split`).

**Ribbons (aggregation).** For each record, anchors are taken in record order. Consecutive anchors join one ribbon
while the B73 chromosome and strand agree, and both the record gap and the B73 gap are <= {max_gap:,} bp. Ribbons
with fewer than {min_ribbon_anchors} anchors are set aside (kept in ribbons.tsv as `drawn=no`). Compatible neighbours
are then re-joined across them, so one out-of-place gene does not cut a block. A drawn ribbon spans the minimum to
maximum coordinate of its anchors on both sides.

**Display orientation.** hifiasm contig orientation is arbitrary. As in GENESPACE, a record whose dominant-chromosome
anchors are mostly on `-` is drawn reverse-complemented (`display_reversed = yes`; table coordinates are never
flipped). A twisted ribbon therefore means **opposite orientation to the rest of its record**: an inversion relative
to B73, or a block of another chromosome joined in reverse. A black mark beside the bar shows the same inverted
dominant-chromosome blocks.

**Dominant chromosome.** The B73 chromosome with the most gene anchors on the record. A tie goes to the larger drawn
record span.
- `frac_anchors_on_dominant` = anchors on it / all anchors on the record.
- `dominant_ref_span_bp` = first to last B73 coordinate of the drawn ribbons on it.
- `dominant_ref_covered_bp` and `dominant_query_covered_bp` = union lengths of those ribbons.
- **Orientation** comes from the dominant-chromosome anchors: forward if >= 90 % are `+`, reverse if <= 10 %, mixed
  otherwise.
- **Inverted blocks** are drawn ribbons on the dominant chromosome whose strand is opposite to the record's majority
  strand.

**Candidate inter-chromosomal record.** The record has >= {sec_min_anchors} anchors **and** >= {sec_min_span:,} bp of
drawn ribbons on each of at least two B73 chromosomes. These thresholds match the secondary-chromosome rule of the
project's chromosome assignment. The column value is `yes`.

**Minor off-chromosome block** (`candidate_interchromosomal = minor`). The record is not `yes`, but it has at least
one drawn ribbon (>= {min_ribbon_anchors} anchors) on a chromosome other than its dominant one.
- The block is **terminal** if it lies within {terminal_bp:,} bp of a record end. A terminal block is the signature of
  a join between chromosomes.
- Otherwise it is **internal**: an inserted or transposed segment, or a duplicated or paralogous gene cluster.
- `off_chr_block_detail` lists each block.

Both tiers flag a record **for review**. It does not establish a misjoin: reference
structural differences, duplicated genes, collapsed repeats and true translocations give the same picture. Read
support at the junction decides.

**Redundancy candidate.** >= {redundant_pct:.0f} % of the record's drawn B73 interval on its dominant chromosome is
also covered by a single other record with more anchors on that chromosome. This marks a possible extra allele or
duplicated copy. It is a flag only.

**Plotted records.** A record is drawn if it has at least one drawn ribbon. All other records stay in sequences.tsv,
as `anchors below ribbon threshold` or `no anchors`.

**Order (GENESPACE-like).**
1. Dominant B73 chromosome (chr1 to chr10).
2. Within that chromosome, the median B73 position of the record's anchors on it.
3. Then longer records first.

`plot_order` in sequences.tsv is that index. Brackets on the ribbon side of each assembly track show the
dominant-chromosome groups.

**Scale.** One bp scale for every panel and figure. The B73 chromosomes are drawn proportional to length, with fixed
gaps.
"""

METHODS_ASM = """
## Assembly order (`--order {order}`)

**Grouping.** Each record takes the chromosome label of the unit's `assignment` table (`assigned_chr`). Records
without a label (or with a label that is not a reference chromosome) form the group `unassigned`, drawn last if they
carry a drawn ribbon. Groups follow reference chromosome order. `assignment_agrees` compares the label with this
run's dominant chromosome.

**Order inside a group.** FASTA (index) order, which is the order of the delivered file. It is **not** a physical
order: no scaffolding, linkage, Hi-C or optical map evidence is used, and the gaps between records are drawn at a
fixed width. The figures are chromosome-grouped contig views, not chromosome-scale assemblies.

**Orientation.** Records keep the orientation of the delivered FASTA. A group is mirrored as a whole (order inside the
group reversed and every record reverse-complemented for display) only when >= {group_flip_pct:.0f} % of the anchors
that its records have on the assigned chromosome are on `-` (`group_reversed`, `group_frac_minus`). Single records are
never flipped, so a record in opposite orientation to B73 shows a twisted ribbon, and local inversions stay visible.

**Off-chromosome emphasis.** A ribbon is outlined when its B73 chromosome differs from the record's assigned
chromosome (dominant chromosome for unassigned records).

**Caveat.** If the assignment table was itself derived from reference anchors, grouping by it is not independent of
the reference. Only the order inside a group and the orientation are assembly-native in that case.
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "run"])
    ap.add_argument("--units", required=True)
    ap.add_argument("--ref-fai", required=True, help="FASTA index of the reference the anchors were made against")
    ap.add_argument("--ref-pattern", default=r"chr\d+", help="reference sequences to show (regex, full match)")
    ap.add_argument("--out", help="new output directory (must not exist)")
    ap.add_argument("--min-ribbon-anchors", type=int, default=3)
    ap.add_argument("--max-gap-bp", type=int, default=5_000_000)
    ap.add_argument("--sec-min-anchors", type=int, default=10)
    ap.add_argument("--sec-min-span-bp", type=int, default=1_000_000)
    ap.add_argument("--terminal-bp", type=int, default=2_000_000,
                    help="an off-chromosome block this close to a record end is 'terminal' (join signature)")
    ap.add_argument("--redundant-frac", type=float, default=0.5)
    ap.add_argument("--label-min-bp", type=int, default=10_000_000)
    ap.add_argument("--order", choices=["synteny", "assembly", "both"], default="synteny")
    ap.add_argument("--group-flip-frac", type=float, default=0.8,
                    help="assembly order: mirror a chromosome group only if this fraction of its anchors is '-'")
    ap.add_argument("--compare", default="", help="comma-separated lines for compare figures (--order both)")
    ap.add_argument("--no-plots", action="store_true", help="tables only")
    a = ap.parse_args(argv)
    units = read_units(a.units)
    chroms, clen = ref_chroms(a.ref_fai, a.ref_pattern)
    if len(chroms) > len(PALETTE):
        die(f"{len(chroms)} reference chromosomes; the palette has {len(PALETTE)}")
    P = params(a)
    compare = [x for x in a.compare.split(",") if x]
    bad = sorted(set(compare) - {u["line"] for u in units})
    if bad:
        die(f"--compare: unknown line(s) {', '.join(bad)}")
    if compare and a.order != "both":
        die("--compare needs --order both")
    if a.command == "check":
        print(f"reference: {len(chroms)} chromosomes ({', '.join(chroms)})", file=sys.stderr)
        load(units, chroms, P)
        print("check OK (nothing written)", file=sys.stderr)
        return 0
    if not a.out:
        die("run needs --out")
    if os.path.exists(a.out):
        die(f"{a.out} exists; choose a new output directory")
    start = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    data = load(units, chroms, P)
    partial = a.out + ".partial"
    if os.path.exists(partial):
        die(f"{partial} exists (an earlier run stopped); move it aside first")
    os.makedirs(partial)
    for u in units:
        same = [x for x in units if x["line"] == u["line"]]
        u["track_label"] = u.get("track_label") or (u["line"] if len(same) == 1 else u["unit"].split(".")[-1])
    seq_rows = [r for u in units for r in sorted(data[u["unit"]]["rows"],
                                                 key=lambda r: (r.get("plot_order") or 10 ** 9, natural(r["record"])))]
    write_tsv(os.path.join(partial, "sequences.tsv"), seq_rows,
              COLS + (ASM_COLS if any(u.get("assignment") for u in units) else []))
    rib_rows = [dict(line=u["line"], unit=u["unit"], record=r["record"], ref_chr=r["rc"], ref_start=r["rs"],
                     ref_end=r["re"], query_start=r["qs"], query_end=r["qe"], strand=r["strand"], gene_anchors=r["n"],
                     drawn="yes" if r["drawn"] else "no") for u in units for r in data[u["unit"]]["ribbons"]]
    write_tsv(os.path.join(partial, "ribbons.tsv"), rib_rows,
              ["line", "unit", "record", "ref_chr", "ref_start", "ref_end", "query_start", "query_end", "strand",
               "gene_anchors", "drawn"])
    summ = [unit_summary(u, data[u["unit"]]["rows"], chroms, data[u["unit"]]["straddle"], data[u["unit"]]["other_ref"])
            for u in units]
    write_tsv(os.path.join(partial, "unit_summary.tsv"), summ, list(summ[0]))
    inputs = []
    for u in units:
        for k in ("anchors", "fai", "coordinate_map", "name_map", "assignment"):
            if u.get(k):
                inputs.append(dict(unit=u["unit"], kind=k, path=os.path.abspath(u[k]), bytes=os.path.getsize(u[k]),
                                   sha256=sha256(u[k])))
    inputs.append(dict(unit="*", kind="ref_fai", path=os.path.abspath(a.ref_fai), bytes=os.path.getsize(a.ref_fai),
                       sha256=sha256(a.ref_fai)))
    inputs.append(dict(unit="*", kind="units", path=os.path.abspath(a.units), bytes=os.path.getsize(a.units),
                       sha256=sha256(a.units)))
    write_tsv(os.path.join(partial, "INPUTS.tsv"), inputs, ["unit", "kind", "path", "bytes", "sha256"])
    modes = ("synteny", "assembly") if a.order == "both" else (a.order,)
    figs = [] if a.no_plots else plot_all(partial, units, data, chroms, clen, P, modes, compare)
    with open(os.path.join(partial, "METHODS.md"), "w") as fh:
        fh.write(METHODS.format(stamp=start, cv=CODE_VERSION, redundant_pct=P["redundant_frac"] * 100, **P))
        if a.order != "synteny":
            fh.write(METHODS_ASM.format(**P))
    run = dict(schema="g2f-riparian/1", code_version=CODE_VERSION, start=start,
               end=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
               command=[sys.executable] + sys.argv, cwd=os.getcwd(), script=os.path.abspath(__file__),
               script_sha256=sha256(__file__), pipeline_git=git_state(), params=P, reference_chromosomes=clen,
               units=[{k: v for k, v in u.items()} for u in units],
               outputs=[os.path.relpath(p, partial) for p in figs] +
               ["sequences.tsv", "ribbons.tsv", "unit_summary.tsv", "INPUTS.tsv", "METHODS.md"],
               python=sys.version.split()[0])
    try:
        import matplotlib
        run["matplotlib"] = matplotlib.__version__
    except ImportError:
        pass
    with open(os.path.join(partial, "RUN.json"), "w") as fh:
        json.dump(run, fh, indent=1)
    os.rename(partial, a.out)
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
