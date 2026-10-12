"""Fingerprinting analysis driver: sketch records -> matrices, assignments, report.

Input: the JSON records written by g2f.fpsketch.sketch (one per read file),
each with label (genotype name; for queries the EXPECTED genotype, may be
empty), tech (hifi | ont) and run_id. HiFi records form the reference panel;
ONT records are the queries. Nothing here moves, renames or relabels files:
the output is evidence for a human decision.
"""
import json
import os
import time

from . import fingerprint as fp

DEFAULTS = {
    "min_count_floor": 2,
    "max_copy": 3.0,
    "interval_z": 3.29,           # two-sided 99.9 % Wilson interval
    # Pairwise comparisons (and private sets) with fewer hashes than this are
    # reported as unresolved instead of scored.
    "min_informative_hashes": 200,
    "min_sensitivity": 0.0,       # queries below this are not used for calibration
    # Haploid genome size: expected depth = bases sketched / genome size,
    # which anchors the single-copy peak of each reference's spectrum.
    "genome_size": 2_300_000_000,
}


def _load_records(paths):
    recs = []
    for p in paths:
        with open(p) as fh:
            r = json.load(fh)
        r["_json"] = os.path.abspath(p)
        recs.append(r)
    return recs


def _query_set(halves, universe=None):
    keys = set()
    for d in halves.values():
        keys.update(d)
    return keys & universe if universe is not None else keys


def _restrict(halves, universe):
    return {h: set(d) & universe for h, d in halves.items()}


def run_k(records, k, params, log=print):
    """Full analysis for one k-mer size. Returns a results dict (JSON-serializable)."""
    P = dict(DEFAULTS, **(params or {}))
    z = P["interval_z"]
    hifi = [r for r in records if r["tech"] == "hifi"]
    ont = [r for r in records if r["tech"] == "ont"]
    by_g = {}
    for r in hifi:
        by_g.setdefault(r["label"], []).append(r)
    if len(by_g) < 2:
        raise ValueError("fingerprinting needs HiFi data of at least two genotypes")

    # ---- pass 1: full references --------------------------------------
    t0 = time.time()
    refs, run_bases = {}, {}
    for g in sorted(by_g):
        counts = {}
        for r in by_g[g]:
            for d in fp.load_sig(r["signature"], k).values():
                fp.add_counts(counts, d)
        run_bases[g] = sum(r["sampling"]["bases_used"] for r in by_g[g])
        refs[g] = fp.Reference(g, counts, P, [r["run_id"] for r in by_g[g]], bases=run_bases[g])
        log(f"[k={k}] reference {g}: {len(by_g[g])} run(s), band {refs[g].band}, "
            f"solid {len(refs[g].solid)}, present {refs[g].n_present}")
        del counts
    panel = fp.Panel(refs)
    log(f"[k={k}] panel built in {time.time() - t0:.0f}s; core {len(panel.core)}; private "
        + ", ".join(f"{g}={len(panel.private[g])}" for g in panel.genotypes))

    # ---- pass 2: HiFi positive/negative controls -----------------------
    controls = []
    for g in sorted(by_g):
        runs = by_g[g]
        loaded = [(r, fp.load_sig(r["signature"], k)) for r in runs]
        variants = []
        if len(runs) >= 2:
            for i, (r, halves) in enumerate(loaded):
                counts = {}
                for j, (_, other) in enumerate(loaded):
                    if j != i:
                        for d in other.values():
                            fp.add_counts(counts, d)
                variants.append(("leave_one_run_out", r["run_id"],
                                 fp.Reference(g, counts, P, [x["run_id"] for j, (x, _) in
                                                             enumerate(loaded) if j != i],
                                              kind=f"loo:{r['run_id']}",
                                              bases=sum(x["sampling"]["bases_used"]
                                                        for j, (x, _) in enumerate(loaded)
                                                        if j != i)),
                                 _query_set(halves)))
        counts_b, query_a = {}, set()
        for _, halves in loaded:
            fp.add_counts(counts_b, halves.get("B", {}))
            query_a.update(halves.get("A", {}))
        variants.append(("split_half", "halfA_vs_halfB",
                         fp.Reference(g, counts_b, P, [r["run_id"] for r in runs], kind="halfB",
                                      bases=sum(r["sampling"]["halves"]["B"]["bases"]
                                                for r in runs)),
                         query_a))
        del loaded, counts_b
        for ctype, qid, ref, query in variants:
            vp = panel.with_variant(g, ref)
            res = fp.score(vp, query & vp.universe, z)
            controls.append({"truth": g, "control_type": ctype, "query": f"{g}:{qid}",
                             "reference_band": ref.band, "result": res})
            sp = fp.support(res, g, P["min_informative_hashes"])
            log(f"[k={k}] control {g}:{qid} ({ctype}) s={res['sensitivity']:.3f} "
                f"has_min={sp['has_min']:.3f} (vs {sp['weakest']}) lacks_max="
                f"{sp['lacks_max']:.3f} (vs {sp['strongest_other']})")

    # ---- ONT queries: runs, replicate halves, pooled per expected label --
    queries = []
    pooled = {}
    for r in sorted(ont, key=lambda x: (x["label"], x["run_id"])):
        halves = _restrict(fp.load_sig(r["signature"], k), panel.universe)
        full = set().union(*halves.values())
        pooled.setdefault(r["label"], set()).update(full)
        queries.append({"query": r["run_id"], "level": "run", "expected": r["label"],
                        "record": r, "result": fp.score(panel, full, z)})
        for h in sorted(halves):
            queries.append({"query": f"{r['run_id']}:half{h}", "level": "half",
                            "expected": r["label"], "record": r,
                            "result": fp.score(panel, halves[h], z)})
    for label in sorted(pooled):
        n_runs = sum(1 for r in ont if r["label"] == label)
        queries.append({"query": f"{label or 'unlabeled'}:all_ont_runs", "level": "pooled",
                        "expected": label, "n_runs": n_runs,
                        "result": fp.score(panel, pooled[label], z)})

    # ---- calibration & classification ---------------------------------
    mn = P["min_informative_hashes"]
    run_level = [q["result"] for q in queries if q["level"] == "run"]
    calib = fp.calibrate([(c["truth"], c["result"]) for c in controls], run_level,
                         P["min_sensitivity"], mn)
    for c in controls:
        c["assignment"] = fp.classify(c["result"], c["truth"], calib, mn)
    for q in queries:
        q["assignment"] = fp.classify(q["result"], q["expected"], calib, mn)

    # consistency across runs and halves of the same expected label
    for label in sorted(pooled):
        grp = [q for q in queries if q["expected"] == label and q["level"] in ("run", "half")]
        bests = {q["assignment"]["best"] for q in grp}
        classes = {q["assignment"]["classification"] for q in grp}
        for q in queries:
            if q["expected"] == label:
                q["consistency"] = {"n_queries": len(grp), "best_matches": sorted(bests),
                                    "classes": sorted(classes),
                                    "consistent_best": len(bests) == 1}

    relatedness = pairwise_relatedness(refs)
    return {"k": k, "params": P, "panel": panel_table(panel, by_g, run_bases),
            "calibration": calib, "controls": controls, "queries": queries,
            "relatedness": relatedness, "genotypes": panel.genotypes}


def panel_table(panel, by_g, run_bases):
    rows = []
    for g in panel.genotypes:
        ref = panel.refs[g]
        rows.append({"genotype": g, "runs": ",".join(ref.runs), "n_runs": len(ref.runs),
                     "bases_sketched": run_bases[g],
                     "expected_depth": round(ref.expected_depth, 2) if ref.expected_depth else None,
                     "spectrum_peak": ref.band["peak"],
                     "solid_lo": ref.band["lo"], "solid_hi": ref.band["hi"],
                     "band_note": ref.band["note"], "n_present": ref.n_present,
                     "n_solid": len(ref.solid), "n_private": len(panel.private[g]),
                     "private_fraction": round(len(panel.private[g]) / len(ref.solid), 5)
                     if ref.solid else 0, "n_core": len(panel.core)})
    return rows


def pairwise_relatedness(refs):
    """Fraction of A's solid k-mers found anywhere in B's reads (asymmetric)."""
    names = sorted(refs)
    out = []
    for a in names:
        sa = refs[a].solid
        for b in names:
            if a == b:
                continue
            shared = sum(1 for h in sa if h in refs[b].present)
            out.append({"a": a, "b": b, "a_solid": len(sa),
                        "fraction_of_a_in_b": round(shared / len(sa), 5) if sa else None})
    return out


# --------------------------------------------------------------------------
# Writing outputs
# --------------------------------------------------------------------------

def _fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.5g}"
    return str(v)


def write_tsv(path, rows, cols):
    with open(path + ".tmp", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c)) for c in cols) + "\n")
    os.replace(path + ".tmp", path)


ASSIGN_COLS = ["query", "level", "tech", "expected", "classification", "best", "best_has_min",
               "best_lacks_max", "best_net", "second", "second_net", "margin",
               "expected_rank", "expected_has_min", "expected_lacks_max", "expected_z",
               "best_private_score", "expected_private_score", "sensitivity", "threshold",
               "unresolved_pairs", "possible_admixture", "reason", "consistent_best", "best_matches_all_runs_halves",
               "bases_sketched"]


def assignment_rows(res):
    rows = []
    for q in res["queries"]:
        a = q["assignment"]
        rows.append(dict(a, query=q["query"], level=q["level"], tech="ont",
                         expected=q["expected"] or "",
                         consistent_best=q.get("consistency", {}).get("consistent_best"),
                         best_matches_all_runs_halves=",".join(
                             q.get("consistency", {}).get("best_matches", [])),
                         bases_sketched=(q["record"]["sampling"]["bases_used"]
                                         if q.get("record") else None)))
    return rows


def control_rows(res):
    rows = []
    for c in res["controls"]:
        a = c["assignment"]
        rows.append(dict(a, query=c["query"], level=c["control_type"], tech="hifi",
                         expected=c["truth"]))
    return rows


def write_k(res, outdir):
    os.makedirs(outdir, exist_ok=True)
    gs = res["genotypes"]
    ont_q = [q for q in res["queries"]]
    for name, key in (("similarity_matrix.tsv", "score"), ("containment_matrix.tsv", "containment")):
        rows = []
        for q in ont_q + [{"query": c["query"], "level": c["control_type"], "expected": c["truth"],
                           "result": c["result"]} for c in res["controls"]]:
            row = {"query": q["query"], "level": q["level"], "expected": q["expected"],
                   "sensitivity": q["result"]["sensitivity"]}
            row.update({g: q["result"]["genotypes"][g][key] for g in gs})
            rows.append(row)
        write_tsv(os.path.join(outdir, name), rows, ["query", "level", "expected", "sensitivity"] + gs)
    long = []
    for src, items in (("ont", ont_q), ("control", res["controls"])):
        for q in items:
            for g in gs:
                v = q["result"]["genotypes"][g]
                long.append(dict(v, query=q["query"], source=src,
                                 level=q.get("level") or q.get("control_type"),
                                 expected=q.get("expected") or q.get("truth"), genotype=g,
                                 sensitivity=q["result"]["sensitivity"],
                                 core_hits=q["result"]["core_hits"], n_core=q["result"]["n_core"]))
    write_tsv(os.path.join(outdir, "scores_long.tsv"), long,
              ["query", "source", "level", "expected", "genotype", "n_private", "private_hits",
               "recovery", "sensitivity", "score", "score_lo", "score_hi", "containment",
               "core_hits", "n_core"])
    pl = []
    for src, items in (("ont", ont_q), ("control", res["controls"])):
        for q in items:
            for (a, b), v in q["result"]["pairs"].items():
                pl.append({"query": q["query"], "source": src,
                           "level": q.get("level") or q.get("control_type"),
                           "expected": q.get("expected") or q.get("truth"), "genotype": a,
                           "versus": b, "n_distinct": v["n"], "hits": v["hits"],
                           "score": v["score"], "lo": v["lo"], "hi": v["hi"]})
    write_tsv(os.path.join(outdir, "pairwise_long.tsv"), pl,
              ["query", "source", "level", "expected", "genotype", "versus", "n_distinct",
               "hits", "score", "lo", "hi"])
    write_tsv(os.path.join(outdir, "assignments.tsv"), assignment_rows(res), ASSIGN_COLS)
    write_tsv(os.path.join(outdir, "controls.tsv"), control_rows(res), ASSIGN_COLS[:-3])
    write_tsv(os.path.join(outdir, "panel.tsv"), res["panel"], list(res["panel"][0]))
    write_tsv(os.path.join(outdir, "hifi_relatedness.tsv"), res["relatedness"],
              ["a", "b", "a_solid", "fraction_of_a_in_b"])
    with open(os.path.join(outdir, "calibration.json"), "w") as fh:
        json.dump(res["calibration"], fh, indent=2)
        fh.write("\n")


def summary_md(results, primary_k, records, params, resources=None):
    """Human-readable report (primary k, with agreement at other k)."""
    res = results[primary_k]
    cal = res["calibration"]
    L = ["# Genotype fingerprinting QC summary", "",
         "Evidence for or against the genotype labels of ONT read sets, from k-mers "
         "private to each HiFi genotype. Nothing was moved or relabelled.", "",
         f"* Primary k = {primary_k}; other k: "
         f"{', '.join(str(k) for k in results if k != primary_k) or 'none'}",
         f"* HiFi panel: {len(res['genotypes'])} genotypes ({', '.join(res['genotypes'])})",
         f"* ONT read sets: {sum(1 for r in records if r['tech'] == 'ont')}", "",
         "## Calibration (from the data)", "",
         "Scores are pairwise: the fraction of genotype G's k-mers that distinguish it from "
         "genotype H which the reads contain, normalized by the reads' recovery of k-mers "
         "shared by all genotypes.", "",
         f"* Positive controls (HiFi reads vs their own genotype represented by independent "
         f"runs or the other replicate half): n = {cal['n_positive']}, min = "
         f"{_fmt(cal['positive_min'])}, median = {_fmt(cal['positive_median'])}",
         f"* ONT best-match positives (label-free): min = {_fmt(cal['ont_best_positive_min'])}, "
         f"median = {_fmt(cal['ont_best_positive_median'])}",
         f"* Negative controls (HiFi reads vs alleles of other genotypes): n = "
         f"{cal['n_negative_controls']}, max = {_fmt(cal['negative_max'])}, median = "
         f"{_fmt(cal['negative_median'])}",
         f"* ONT null (each ONT run vs genotypes other than its best match; reported, not used "
         f"for the threshold): n = {cal['n_negative_ont']}, median = "
         f"{_fmt(cal['ont_negative_median'])}, 95th percentile = {_fmt(cal['ont_negative_p95'])}, "
         f"max = {_fmt(cal['negative_max_ont'])}",
         f"* Discriminating: **{cal['discriminating']}** ({cal['reason']}); threshold = "
         f"{_fmt(cal['threshold'])}", "",
         "## ONT assignments", "",
         "has = lowest pairwise score of the genotype's distinct alleles vs any other "
         "genotype (~1 if the reads are that genotype); lacks = highest score of another "
         "genotype's distinct alleles (~0 if so). sens. = query sensitivity (core k-mers "
         "recovered).", "",
         "| query | level | expected | best | has | lacks | expected: has / lacks | sens. "
         "| class | reason |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for q in res["queries"]:
        a = q["assignment"]
        L.append(f"| {q['query']} | {q['level']} | {q['expected'] or '-'} | {a['best']} | "
                 f"{_fmt(a['best_has_min'])} | {_fmt(a['best_lacks_max'])} | "
                 f"{_fmt(a['expected_has_min'])} / {_fmt(a['expected_lacks_max'])} | "
                 f"{_fmt(a['sensitivity'])} | **{a['classification']}** | {a['reason']} |")
    L += ["", "## HiFi controls (truth known)", "",
          "| control | type | has | lacks | sens. | class |", "|---|---|---|---|---|---|"]
    for c in res["controls"]:
        a = c["assignment"]
        L.append(f"| {c['query']} | {c['control_type']} | {_fmt(a['expected_has_min'])} | "
                 f"{_fmt(a['expected_lacks_max'])} | {_fmt(a['sensitivity'])} | "
                 f"{a['classification']} |")
    L += ["", "## Panel", "", "| genotype | runs | bases | peak | solid band | solid | private |",
          "|---|---|---|---|---|---|---|"]
    for p in res["panel"]:
        L.append(f"| {p['genotype']} | {p['n_runs']} | {p['bases_sketched'] / 1e9:.1f} Gb | "
                 f"{p['spectrum_peak']} | {p['solid_lo']}-{p['solid_hi']} | {p['n_solid']} | "
                 f"{p['n_private']} |")
    if len(results) > 1:
        L += ["", "## Agreement across k", ""]
        for q in res["queries"]:
            cls = {k: next(x for x in results[k]["queries"] if x["query"] == q["query"])
                   ["assignment"]["classification"] for k in results}
            if len(set(cls.values())) > 1:
                L.append(f"* {q['query']}: " + ", ".join(f"k={k}: {v}" for k, v in cls.items()))
        if L[-1] == "":
            L.append("All queries have the same classification at every k.")
    if resources:
        L += ["", "## Resources", "", "```", json.dumps(resources, indent=2), "```"]
    L += ["", "Interpretation and method: docs/fingerprinting.md", ""]
    return "\n".join(L)


def analyze(record_paths, outdir, ksizes, primary_k, params, provenance=None, log=print):
    t0 = time.time()
    records = _load_records(record_paths)
    for r in records:
        if r["tech"] not in ("hifi", "ont"):
            raise ValueError(f"{r['_json']}: tech must be hifi or ont")
    os.makedirs(outdir, exist_ok=True)
    results = {}
    for k in ksizes:
        results[k] = run_k(records, k, params, log=log)
        write_k(results[k], os.path.join(outdir, f"k{k}"))
    # top-level primary assignments with cross-k agreement
    rows = assignment_rows(results[primary_k])
    for row in rows:
        for k in ksizes:
            if k != primary_k:
                other = next(x for x in results[k]["queries"] if x["query"] == row["query"])
                row[f"class_k{k}"] = other["assignment"]["classification"]
    cols = ASSIGN_COLS + [f"class_k{k}" for k in ksizes if k != primary_k]
    write_tsv(os.path.join(outdir, "assignments.tsv"), rows, cols)
    gs = results[primary_k]["genotypes"]
    sim = [dict({"query": q["query"], "level": q["level"], "expected": q["expected"],
                 "sensitivity": q["result"]["sensitivity"]},
                **{g: q["result"]["genotypes"][g]["score"] for g in gs})
           for q in results[primary_k]["queries"]]
    write_tsv(os.path.join(outdir, "similarity_matrix.tsv"), sim,
              ["query", "level", "expected", "sensitivity"] + gs)
    import resource
    usage = {"wall_s": round(time.time() - t0, 1),
             "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    prov = dict(provenance or {})
    prov.update({
        "ksizes": list(ksizes), "primary_k": primary_k,
        "analysis_params": dict(DEFAULTS, **(params or {})),
        "inputs": [{k: r[k] for k in ("label", "tech", "run_id", "input", "format", "signature",
                                      "sourmash_version", "params", "sampling", "runtime_s")}
                   | {"sketch_record": r["_json"]} for r in records],
        "calibration": {k: results[k]["calibration"] for k in ksizes},
        "analysis_resources": usage,
    })
    with open(os.path.join(outdir, "parameters_and_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2, default=str)
        fh.write("\n")
    sketch_res = {"sketch_cpu_s_total": round(sum(r.get("runtime_s", 0) for r in records), 1),
                  "bases_sketched_total": sum(r["sampling"]["bases_used"] for r in records)}
    with open(os.path.join(outdir, "qc_summary.md"), "w") as fh:
        fh.write(summary_md(results, primary_k, records, params, dict(usage, **sketch_res)))
    return results
