"""Per-sample consolidated report from validated stage outputs only.

Writes a new directory report/<stamp>/ each time (never overwrites) with
metrics.tsv (long format: assembly, section, metric, value, source),
report.md, report.html and figures/. Missing or unvalidated stages are listed
as such; nothing is inferred from file existence alone.
"""
import csv
import datetime
import html
import json
import os
import re
import shutil

from . import plots
from .runner import readtext, stage_state
from .stages import asm_dir, pair_dir

LIMITATIONS = [
    "Read-depth classes are relative to one depth baseline c (homozygous / collapsed-haplotype depth, shared by all assemblies of a sample); GC- or repeat-related coverage bias "
    "can mimic low support or collapse.",
    "HiFi reads (N50 typically 13-20 kb) cannot span maize LTR-retrotransposon arrays and long tandem repeats; "
    "breakpoints and low continuity inside such repeats are reported as unresolved, not as errors.",
    "Differences from the reference are expected: maize inbreds differ from B73 by large presence/absence "
    "variation, inversions and transposon insertions. A junction is called a suspected misassembly only "
    "when read continuity also breaks.",
    "Mixed-site density is computed from the MD tag of unique alignments (substitutions and deletion starts); "
    "insertions are not counted.",
    "Graph depth classes use the hifiasm rd tag relative to hifiasm's homozygous peak (peak_hom); in inbred "
    "lines the rd mode matches peak_hom, but rd is not the same quantity as read-alignment depth.",
    "Graph neighbourhoods are drawn from noseq unitig graphs; segment names are specific to one hifiasm run.",
    "Read-retention comparison between two hifiasm runs assumes both used the same reads and error correction "
    "(checked: fraction of B contig reads also present in A).",
]


def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except OSError:
        return None


def _validated(cfg, stages, sid):
    st = stages.get(sid)
    return st is not None and stage_state(st)[0] == "validated"


def existing_metrics(cfg, a):
    """Row of a project metrics TSV selected by the assembly's existing_metrics keys."""
    path, keys = cfg.get("existing_metrics_tsv"), a.get("existing_metrics") or {}
    if not path or not keys or not os.path.exists(path):
        return {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if all(row.get(k) == str(v) for k, v in keys.items()):
                return row
    return {}


EXISTING = [("busco_odb12.2_C", "BUSCO poales_odb12.2 C %"), ("busco_odb12.2_S", "BUSCO poales_odb12.2 S %"),
            ("busco_odb12.2_D", "BUSCO poales_odb12.2 D %"), ("busco_odb12.2_M", "BUSCO poales_odb12.2 M %"),
            ("busco_odb10_C", "BUSCO poales_odb10 C %"), ("busco_odb10_D", "BUSCO poales_odb10 D %"),
            ("qv", "Merqury QV"), ("kmer_completeness_pct", "Merqury k-mer completeness %"),
            ("homdepth_kmers_2copies_pct", "hom-depth k-mers found twice %"),
            ("homdepth_kmers_absent_pct", "hom-depth k-mers absent %"),
            ("kmer_genome_est_bp", "k-mer genome size estimate (bp)"),
            ("p_ctg_vs_kmer_est_pct", "size vs k-mer estimate %")]


def collect(cfg, stages, sname):
    s = cfg["samples"][sname]
    rows, status = [], []

    def add(asm, section, metric, value, source):
        rows.append({"assembly": asm, "section": section, "metric": metric, "value": value, "source": source})

    for aname, a in s["assemblies"].items():
        lab = a["label"]
        for sid in [x for x in stages if x.split(":")[1:3] == [sname, aname] and not x.startswith("pair_")]:
            status.append((sid, *stage_state(stages[sid])))
        if _validated(cfg, stages, f"prepare:{sname}:{aname}"):
            st = _load(os.path.join(asm_dir(cfg, sname, aname, "prepare"), "stats.json"))
            for k in ("bp", "contigs", "N50", "L50", "N90", "L90", "largest"):
                add(lab, "contiguity", k, st.get(k), "prepare/stats.json")
            log = st.get("assembler_log", {})
            for k in ("version", "hom_cov_threshold", "purge_cov_threshold", "peak_hom", "peak_het"):
                if k in log:
                    add(lab, "assembler", k, ",".join(log[k]) if isinstance(log[k], list) else log[k],
                        "assembler stderr log")
        ex = existing_metrics(cfg, a)
        for k, name in EXISTING:
            if ex.get(k) not in (None, "", "NA"):
                add(lab, "existing QC", name, ex[k], os.path.basename(cfg["existing_metrics_tsv"] or ""))
        if _validated(cfg, stages, f"map_hifi:{sname}:{aname}"):
            v = _load(os.path.join(asm_dir(cfg, sname, aname, "map_hifi"), "VALIDATED.json"))["validation"]
            add(lab, "read support", "primary reads mapped %", v["primary_mapped_pct"], "map_hifi flagstat")
            add(lab, "read support", "supplementary alignments", v["supplementary"], "map_hifi flagstat")
        if _validated(cfg, stages, f"regions:{sname}:{aname}"):
            r = _load(os.path.join(asm_dir(cfg, sname, aname, "regions"), "summary.json"))
            add(lab, "read support", "depth baseline c (homozygous depth)", r["expected_depth"], r["expected_depth_source"])
            for k, v in r["class_bp"].items():
                add(lab, "read depth class (bp)", k, v, "regions/window_classes")
            add(lab, "read support", "multi-mapping windows (bp)", r["multimapping_bp"], "regions")
            add(lab, "read support", "mixed sites per Mb tested", round(r["mixed_sites_total"] / max(1, r["tested_bp"]) * 1e6, 1), "read_support")
            add(lab, "read support", "consensus-error sites per Mb tested", round(r["error_sites_per_mb"], 2), "read_support")
            for k, v in sorted(r["breakpoints"].items()):
                add(lab, "breakpoints (n)", k, v["n"], "regions/breakpoints.tsv")
        if _validated(cfg, stages, f"graph:{sname}:{aname}"):
            g = _load(os.path.join(asm_dir(cfg, sname, aname, "graph"), "summary.json"))
            for k in ("segments", "links", "tips", "branching_ends", "simple_bubbles", "bubble_branch_bp", "rd_mode"):
                add(lab, "unitig graph", k, g["utg"][k], "graph/summary.json")
            add(lab, "unitig graph", "rd reference (homozygous depth)", g["rd_reference"], g["rd_reference_source"])
            add(lab, "unitig graph", "contig links (p_ctg GFA)", g["ctg_links"], "graph/summary.json")
            add(lab, "unitig graph", "unlinked contig path transitions", g["unlinked_transitions"], "graph")
            for k, v in g["contig_bp_by_utg_rd_class"].items():
                add(lab, "contig bp by unitig rd class", k, v, "graph")
        for rname in cfg["references"]:
            if _validated(cfg, stages, f"ref_analysis:{sname}:{aname}:{rname}"):
                r = _load(os.path.join(asm_dir(cfg, sname, aname, f"ref_{rname}_analysis"), "summary.json"))
                sec = f"vs {rname}"
                add(lab, sec, "assembly bp assigned to a chromosome", r["assigned_bp"], "ref_analysis")
                add(lab, sec, "reference chromosome bp covered %", round(100 * r["reference_covered_frac"], 2), "ref_analysis")
                add(lab, sec, "multi-chromosome contigs", r["multi_chromosome_contigs"], "ref_analysis")
                for k, v in r["junctions"].items():
                    add(lab, sec, f"junctions: {k}", v, "ref_analysis")
                add(lab, sec, "reference bp in uncovered runs", r["absent_bp"], "ref_analysis")
        if _validated(cfg, stages, f"evidence:{sname}:{aname}"):
            e = _load(os.path.join(asm_dir(cfg, sname, aname, "evidence"), "summary.json"))
            for k, v in sorted(e["by_category"].items()):
                add(lab, "candidates (n)", k, v["n"], "evidence/candidates.tsv")
                add(lab, "candidates (bp)", k, v["bp"], "evidence/candidates.tsv")
    pairs = []
    for pr in s["pairs"]:
        sid = f"pair_analysis:{sname}:{pr['a']}:{pr['b']}"
        status.append((sid, *stage_state(stages[sid])) if sid in stages else (sid, "absent", ""))
        if _validated(cfg, stages, sid):
            d = pair_dir(cfg, sname, pr["a"], pr["b"], "analysis")
            pairs.append((pr, _load(os.path.join(d, "summary.json")), _load(os.path.join(d, "read_retention_summary.json")), d))
    return rows, status, pairs


def _figures(cfg, stages, sname, figdir):
    """Comparative figures across the sample's assemblies (validated inputs only)."""
    s = cfg["samples"][sname]
    figs = []
    hs, bars, cov = [], [], []
    for aname, a in s["assemblies"].items():
        if _validated(cfg, stages, f"regions:{sname}:{aname}"):
            r = _load(os.path.join(asm_dir(cfg, sname, aname, "regions"), "summary.json"))
            h = [l.split("\t") for l in readtext(os.path.join(asm_dir(cfg, sname, aname, "read_support"), "depth_hist.tsv")).splitlines()[1:]]
            hs.append((a["label"], [int(x[0]) for x in h], [int(x[1]) for x in h], r["expected_depth"]))
            bars.append((a["label"], r["class_bp"]))
    if hs:
        plots.depth_histograms(hs, os.path.join(figdir, "depth_histograms"))
        plots.class_bars(bars, os.path.join(figdir, "depth_classes"))
        figs += [("Read-depth distribution", "depth_histograms.svg"), ("Sequence by read-depth class", "depth_classes.svg")]
    for rname, r in cfg["references"].items():
        panels, chroms = [], None
        for aname, a in s["assemblies"].items():
            d = asm_dir(cfg, sname, aname, f"ref_{rname}_analysis")
            if not _validated(cfg, stages, f"ref_analysis:{sname}:{aname}:{rname}"):
                continue
            shutil.copy(os.path.join(d, "dotplot_genome.svg"), os.path.join(figdir, f"dotplot_{rname}_{aname}.svg"))
            figs.append((f"Dot plot vs {r['label']}: {a['label']}", f"dotplot_{rname}_{aname}.svg"))
            c = {}
            for line in readtext(os.path.join(d, "reference_coverage.tsv")).splitlines()[1:]:
                p = line.split("\t")
                c.setdefault(p[0], []).append((int(p[1]), float(p[3])))
            cov.append((a["label"], c))
            chroms = sorted(c, key=plots._chrom_key)
            from .refalign import load_alignments
            assign = {}
            for line in readtext(os.path.join(d, "contig_assignment.tsv")).splitlines()[1:]:
                p = line.split("\t")
                assign[p[0]] = {"contig": p[0], "length": int(p[1]), "chrom": p[4], "ref_median_pos": int(p[7])}
            by_q = load_alignments(os.path.join(asm_dir(cfg, sname, aname, f"ref_{rname}_align"), "aln.paf.gz"), cfg["science"])
            panels.append((a["label"], by_q, assign))
        if cov:
            plots.coverage_tracks(cov, chroms, os.path.join(figdir, f"refcov_{rname}"), r["label"])
            figs.append((f"{r['label']} covered by each assembly", f"refcov_{rname}.svg"))
        if len(panels) >= 2:
            from .stages import _ref_lengths_and_gaps
            lengths, _ = _ref_lengths_and_gaps(cfg, r)
            for c in chroms:
                plots.dotplot_chrom_pair(c, lengths[c], panels, cfg["science"]["ref_align"]["plot_min_len"],
                                         r["label"], os.path.join(figdir, f"dotplot_{rname}_{c}"))
                figs.append((f"{c}: assemblies side by side vs {r['label']}", f"dotplot_{rname}_{c}.svg"))
    return figs


def _md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def _fmt(v):
    if isinstance(v, float):
        return f"{v:,.3f}".rstrip("0").rstrip(".")
    if isinstance(v, int) or (isinstance(v, str) and re.fullmatch(r"\d{5,}", v)):
        return f"{int(v):,}"
    return str(v)


def run(cfg, stages, sname, stamp):
    out = os.path.join(cfg["root"], "samples", sname, "report", stamp)
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir)
    rows, status, pairs = collect(cfg, stages, sname)
    labels = [a["label"] for a in cfg["samples"][sname]["assemblies"].values()]
    with open(os.path.join(out, "metrics.tsv"), "w") as fh:
        fh.write("sample\tassembly\tsection\tmetric\tvalue\tsource\n")
        for r in rows:
            fh.write(f"{sname}\t{r['assembly']}\t{r['section']}\t{r['metric']}\t{r['value']}\t{r['source']}\n")
    # wide table
    wide, order = {}, []
    for r in rows:
        k = (r["section"], r["metric"])
        if k not in wide:
            wide[k] = {}
            order.append(k)
        wide[k][r["assembly"]] = r["value"]
    figs = _figures(cfg, stages, sname, figdir)
    md = [f"# Assembly validation report: {sname}", "",
          f"Generated {datetime.datetime.now().astimezone().isoformat(timespec='seconds')} "
          f"from profile `{cfg['profile']}` (science hash `{cfg['science_hash']}`). "
          "Values come only from validated stages; see Stage status.", "",
          "## Metrics side by side", "",
          _md_table(["section", "metric"] + labels,
                    [[k[0], k[1]] + [_fmt(wide[k].get(l, "")) for l in labels] for k in order])]
    for pr, ps, ret, d in pairs:
        a_l = cfg["samples"][sname]["assemblies"][pr["a"]]["label"]
        b_l = cfg["samples"][sname]["assemblies"][pr["b"]]["label"]
        for f in ("sequence_classes.svg", "redundant_identity.svg"):
            if os.path.exists(os.path.join(d, f)):
                shutil.copy(os.path.join(d, f), os.path.join(figdir, f"pair_{pr['a']}_{pr['b']}_{f}"))
                figs.append((f"{a_l} vs {b_l}: {f[:-4].replace('_', ' ')}", f"pair_{pr['a']}_{pr['b']}_{f}"))
        md += ["", f"## {a_l} sequence classified against {b_l}", "",
               _md_table(["class", "bp", f"% of {a_l}"],
                         [[k, _fmt(v), f"{100 * v / ps['a_bp']:.2f}"] for k, v in ps["class_bp"].items()]),
               "", f"{b_l} sequence not covered by any {a_l} alignment: {_fmt(ps['b_only_bp'])} bp.", "",
               "Read-depth cross-tabulation (bp; A = class in "
               f"{a_l}, B = class of the matching {b_l} region):", "",
               _md_table(["class", "bp"], sorted([[k, _fmt(v)] for k, v in ps["depth_crosstab_bp"].items()],
                                                 key=lambda x: x[0]))]
        if ret:
            md += ["", "Read retention (hifiasm contig read lists):", "",
                   _md_table(["metric", "value"], [[k, json.dumps(v) if isinstance(v, dict) else _fmt(v)]
                                                   for k, v in ret.items()])]
    for aname, a in cfg["samples"][sname]["assemblies"].items():
        cand = os.path.join(asm_dir(cfg, sname, aname, "evidence"), "candidates.tsv")
        if _validated(cfg, stages, f"evidence:{sname}:{aname}"):
            shutil.copy(cand, os.path.join(out, f"candidates_{aname}.tsv"))
            top = [l.split("\t") for l in readtext(cand).splitlines()[1:]]
            top = [t for t in top if t[4] not in ("retained_haplotypic_copy",)][:25]
            md += ["", f"## Top candidates: {a['label']}", "",
                   f"Full table: `candidates_{aname}.tsv`. Expected allelic copies (retained_haplotypic_copy) "
                   "are excluded from this list.", "",
                   _md_table(["id", "location", "category", "confidence", "read evidence", "reference", "graph"],
                             [[t[0], f"{t[1]}:{int(t[2]):,}-{int(t[3]):,}", t[4], t[5], t[8], t[9], t[10]] for t in top])]
            nb = os.path.join(asm_dir(cfg, sname, aname, "evidence"), "neighbourhoods")
            if os.path.isdir(nb):
                shutil.copytree(nb, os.path.join(figdir, f"neighbourhoods_{aname}"))
    md += ["", "## Figures", ""] + [f"- {t}: `figures/{f}`" for t, f in figs]
    md += ["", "## Limitations", ""] + [f"- {x}" for x in LIMITATIONS]
    md += ["", "## Stage status", "", _md_table(["stage", "state", "detail"], status)]
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    _html(out, sname, md, figs)
    return out


def _html(out, sname, md_lines, figs):
    """Minimal self-contained HTML rendering of the Markdown tables + figures."""
    body, in_table = [], False
    for line in "\n".join(md_lines).splitlines():
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            if set("".join(cells)) <= set("-"):
                continue
            if not in_table:
                body.append("<table>")
                in_table = True
                body.append("<tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in cells) + "</tr>")
            else:
                body.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in cells) + "</tr>")
            continue
        if in_table:
            body.append("</table>")
            in_table = False
        if line.startswith("## Figures"):
            body.append("<h2>Figures</h2>")
            for t, f in figs:
                body.append(f"<figure><img src='figures/{html.escape(f)}' alt='{html.escape(t)}'>"
                            f"<figcaption>{html.escape(t)}</figcaption></figure>")
            continue
        if line.startswith("- ") and ": `figures/" in line:
            continue
        if line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("- "):
            body.append(f"<li>{html.escape(line[2:])}</li>")
        elif line.strip():
            body.append(f"<p>{html.escape(line)}</p>")
    if in_table:
        body.append("</table>")
    css = ("body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 16px;color:#0b0b0b;"
           "background:#fcfcfb}table{border-collapse:collapse;font-size:12px;margin:1em 0;display:block;overflow-x:auto}"
           "td,th{border:1px solid #e4e3df;padding:3px 6px;text-align:left}th{background:#f1f0ec}"
           "figure{margin:1.5em 0}img{max-width:100%}figcaption{font-size:12px;color:#52514e}")
    with open(os.path.join(out, "report.html"), "w") as fh:
        fh.write(f"<!doctype html><html><head><meta charset='utf-8'><title>Assembly validation {html.escape(sname)}"
                 f"</title><style>{css}</style></head><body>{''.join(body)}</body></html>\n")
