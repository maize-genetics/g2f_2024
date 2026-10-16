"""Combine read-support, reference-alignment and graph evidence into one
ranked candidate table per assembly, and draw graph neighbourhoods for the
top candidates.

Final categories (never "confirmed error" from one signal):
  suspected_misassembly        clipped-read breakpoint without spanning reads, or
                               a reference junction where read continuity breaks
  likely_biological_sv         reference junction continuously spanned by reads,
                               or a breakpoint spanned by about half the reads
  collapsed_repeat_candidate   ~2c or higher depth (graded by mixed-site enrichment)
  false_duplication_candidate  half depth with multi-mapping reads (near-identical copies)
  retained_haplotypic_copy     half depth with unique reads (an allele kept twice;
                               expected in an unpurged assembly of a heterozygous sample)
  low_support                  little or no read support
  unresolved                   signals conflict or sit in repeats longer than the reads
"""
import gzip
import json
import os

from . import graph as G

CONF_SCORE = {"high": 3, "medium": 2, "low": 1, "info": 0}
CAT_WEIGHT = {"suspected_misassembly": 5, "collapsed_repeat_candidate": 3, "false_duplication_candidate": 3,
              "low_support": 2, "unresolved": 2, "likely_biological_sv": 1, "retained_haplotypic_copy": 0}


def _tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


def _windows(path):
    by = {}
    with gzip.open(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            by.setdefault(r["contig"], []).append(r)
    return by


def _overlapping(rows, s, e):
    return [r for r in rows if int(r["start"]) < e and int(r["end"]) > s]


def junction_support(j, windows, bps_by_contig, c, sci):
    """Grade a reference junction by read continuity across it."""
    f = sci["evidence"]["junction_flank"]
    s, e = int(j["q_start"]) - f, int(j["q_end"]) + f
    ws = _overlapping(windows.get(j["contig"], []), s, e)
    if not ws:
        return "unresolved", "low", "no read-support windows (alignment stage not run?)"
    span_min = min(int(float(w["span_min"])) for w in ws)
    mq0 = max(1 - float(w["depth_hq"]) / float(w["depth_all"]) if float(w["depth_all"]) else 0 for w in ws)
    mis = [b for b in bps_by_contig.get(j["contig"], []) if s <= int(b["pos"]) <= e
           and b["call"] == "suspected_misjoin"]
    ev = f"min spanning reads {span_min} ({span_min / c:.2f}c) over {e - s} bp; max multi-mapping {mq0:.2f}"
    if mis:
        ev += f"; {len(mis)} suspected-misjoin clip cluster(s) inside"
    if span_min >= 0.3 * c and not mis:
        return "likely_biological_sv", "high", ev + "; reads continuous across the junction"
    if span_min == 0 and mis:
        return "suspected_misassembly", "high", ev
    if span_min <= 0.1 * c and mq0 >= sci["regions"]["mq0_high"]:
        return "unresolved", "medium", ev + "; repeat longer than reads at the junction"
    if span_min <= 0.1 * c:
        return "suspected_misassembly", "medium", ev
    return "unresolved", "low", ev


REGION_CATS = {"collapsed_repeat_or_duplicate": "collapsed_repeat_candidate",
               "near_identical_duplicate": "false_duplication_candidate",
               "haploid_depth_region": "retained_haplotypic_copy",
               "unsupported_sequence": "low_support", "low_support": "low_support",
               "low_support_contig_end": "low_support"}
BP_CATS = {"suspected_misjoin": "suspected_misassembly", "likely_heterozygous_sv": "likely_biological_sv",
           "weakly_spanned_breakpoint": "unresolved"}


def run(regions_dir, read_dir, ref_dirs, graph_dir, utg_gfa, sci, out_dir, label):
    with open(os.path.join(regions_dir, "summary.json")) as fh:
        c = json.load(fh)["expected_depth"]
    regions = _tsv(os.path.join(regions_dir, "regions.tsv"))
    bps = _tsv(os.path.join(regions_dir, "breakpoints.tsv"))
    windows = _windows(os.path.join(read_dir, "windows.tsv.gz"))
    bps_by = {}
    for b in bps:
        bps_by.setdefault(b["contig"], []).append(b)
    paths = G.load_paths(os.path.join(graph_dir, "contig_utg_path.tsv")) if graph_dir else {}
    trans = {}
    if graph_dir:
        for t in _tsv(os.path.join(graph_dir, "path_transitions.tsv")):
            trans.setdefault(t["contig"], []).append(t)

    cands = []
    for r in regions:
        if r["confidence"] == "low" and r["call"] != "unsupported_sequence":
            continue
        cands.append({"contig": r["contig"], "start": int(r["start"]), "end": int(r["end"]), "source": "reads",
                      "signal": r["call"], "category": REGION_CATS[r["call"]], "confidence": r["confidence"],
                      "read_evidence": r["evidence"], "ref_evidence": "", "graph_evidence": ""})
    for b in bps:
        if b["call"] not in BP_CATS or b["confidence"] not in ("high", "medium"):
            continue
        p = int(b["pos"])
        cands.append({"contig": b["contig"], "start": max(0, p - 1000), "end": p + 1000, "source": "reads",
                      "signal": b["call"], "category": BP_CATS[b["call"]], "confidence": b["confidence"],
                      "read_evidence": b["evidence"], "ref_evidence": "", "graph_evidence": "",
                      "_sa_partner": b.get("sa_partner", "NA")})
    for ref, d in ref_dirs.items():
        for j in _tsv(os.path.join(d, "junctions.tsv")):
            cat, conf, ev = junction_support(j, windows, bps_by, c, sci)
            cands.append({"contig": j["contig"], "start": int(j["q_start"]), "end": int(j["q_end"]),
                          "source": f"ref:{ref}", "signal": j["kind"], "category": cat, "confidence": conf,
                          "read_evidence": ev,
                          "ref_evidence": f"{j['kind']} vs {ref}: {j['left']} | {j['right']}", "graph_evidence": ""})
    # cross-support: a read breakpoint inside a reference junction raises both
    for x in cands:
        if x["source"] == "reads" and x["category"] == "suspected_misassembly":
            hits = [y for y in cands if y["source"].startswith("ref:") and y["contig"] == x["contig"]
                    and y["start"] - 50000 <= x["start"] <= y["end"] + 50000]
            if hits:
                x["ref_evidence"] = "; ".join(h["ref_evidence"] for h in hits[:2])
                x["confidence"] = "high"
    # graph context
    for x in cands:
        if not paths:
            break
        segs = G.segments_in(paths, x["contig"], x["start"], x["end"])
        tr = [t for t in trans.get(x["contig"], []) if x["start"] - 20000 <= int(t["pos"]) <= x["end"] + 20000]
        unl = sum(1 for t in tr if t["linked"] == "0")
        br = sum(1 for t in tr if int(t["exit_degree"]) >= 2)
        x["graph_evidence"] = (f"{len(segs)} unitig(s) {','.join(segs[:4])}{'...' if len(segs) > 4 else ''}; "
                               f"{len(tr)} path transition(s), {unl} unlinked, {br} at branching ends")
        x["_segs"] = segs
        if x["category"] == "suspected_misassembly" and unl:
            x["graph_evidence"] += "; contig path jumps between unlinked unitigs here"
    # split reads whose other part aligns inside a collapsed region: the junction a
    # missing repeat copy leaves (any distance away, even on another contig)
    collapses = [y for y in cands if y["category"] == "collapsed_repeat_candidate"]
    for x in cands:
        partner = x.get("_sa_partner", "NA")
        if x["category"] != "suspected_misassembly" or partner == "NA":
            continue
        pc, pb = partner.rsplit(":", 1)
        pb = int(pb)
        hit = [y for y in collapses if y["contig"] == pc and y["start"] - 10000 <= pb <= y["end"]]
        if hit:
            x["read_evidence"] += (f"; split reads continue in collapsed region {pc}:{hit[0]['start']:,}-"
                                   f"{hit[0]['end']:,}: consistent with a missing repeat copy")
    for x in cands:
        x["n_evidence_types"] = (1 + bool(x["ref_evidence"] and x["source"] == "reads")
                                 + ("unlinked" in x["graph_evidence"] and "0 unlinked" not in x["graph_evidence"]))
        x["score"] = (CONF_SCORE[x["confidence"]] * 10 + CAT_WEIGHT[x["category"]] * 3 + 8 * (x["n_evidence_types"] - 1)
                      + min(5, (x["end"] - x["start"]) / 1e5))
    cands.sort(key=lambda x: -x["score"])
    for i, x in enumerate(cands, 1):
        x["id"] = f"{label}-C{i:05d}"

    # neighbourhood figures for the top candidates (excluding expected haplotypic copies)
    figs = 0
    if paths and utg_gfa:
        utg = G.parse_gfa(utg_gfa, want_a=False)
        with open(os.path.join(graph_dir, "summary.json")) as fh:
            mode = json.load(fh)["rd_reference"]
        os.makedirs(os.path.join(out_dir, "neighbourhoods"), exist_ok=True)
        top = [x for x in cands if x["category"] != "retained_haplotypic_copy" and x.get("_segs")]
        for x in top[: sci["graph"]["top_candidates"]]:
            seeds = set(x["_segs"])
            nodes = G.neighborhood(utg, seeds, sci["graph"]["neighborhood_radius"], sci["graph"]["max_nodes"])
            path_nodes = set(r[3] for r in paths.get(x["contig"], []))
            base = os.path.join(out_dir, "neighbourhoods", x["id"])
            G.write_subgraph(utg, nodes, base + ".gfa")
            G.render(utg, nodes, seeds, path_nodes, mode,
                     f"{x['id']} {x['category']} ({x['confidence']})\n{x['contig']}:{x['start']:,}-{x['end']:,}",
                     base)
            x["figure"] = f"neighbourhoods/{x['id']}.svg"
            figs += 1
    cols = ["id", "contig", "start", "end", "category", "confidence", "n_evidence_types", "source", "signal", "read_evidence",
            "ref_evidence", "graph_evidence", "figure"]
    with open(os.path.join(out_dir, "candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in cands:
            fh.write("\t".join(str(x.get(k, "")) for k in cols) + "\n")
    tally = {}
    for x in cands:
        k = f"{x['category']}|{x['confidence']}"
        tally.setdefault(k, {"n": 0, "bp": 0})
        tally[k]["n"] += 1
        tally[k]["bp"] += x["end"] - x["start"]
    summary = {"candidates": len(cands), "by_category": tally, "neighbourhood_figures": figs,
               "expected_depth": c}
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary
