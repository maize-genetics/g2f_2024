#!/usr/bin/env python3
"""Summary tables for the supplementary fingerprinting diagnostics (read-only).

Reads an existing fingerprint report (pairwise_long.tsv, assignments.tsv) and
the outputs of g2f_fp_diagnostics.py, and writes small TSVs:

  ont_null_calibration.tsv   ONT platform-matched null of pairwise "foreign
                             allele" scores, from genotypes whose ONT data show
                             no excess (leave-genotype-out for the queries
                             tested), by sensitivity bin and k
  foreign_scores_by_k.tsv    per query x other genotype foreign scores at
                             k = 21 and k = 31, with each value's excess over
                             the ONT null (p99 of the matching sensitivity bin)
  read_cooccurrence.tsv      per-read clustering of foreign markers
  run_sensitivity.tsv        has/lacks vs sensitivity, read length and bases
  generic_foreign_excess.tsv foreign k-mers above the clean-line ONT null, and
                             the dosage of that excess (pooled queries)
  detection_limits_k21.tsv   indicative smallest detectable whole-genome
                             admixture fraction per (pooled query, panel line)

Standard library only.
"""
import argparse
import csv
import glob
import json
import math
import os
import statistics
from math import comb

BINS = [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 1.01)]


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def wt(path, rows, cols=None):
    cols = cols or list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("NA" if r.get(c) is None else
                               (f"{r[c]:.5g}" if isinstance(r[c], float) else str(r[c]))
                               for c in cols) + "\n")


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(math.ceil(q * len(xs))) - 1)] if xs else None


def binof(s):
    for lo, hi in BINS:
        if lo <= s < hi:
            return f"{lo:.2f}-{min(hi, 1):.2f}"
    return None


def zt_lambda(mean):
    if mean is None or mean <= 1.0:
        return 0.0
    lo, hi = 1e-9, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid / (1 - math.exp(-mid)) < mean:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def zt_mean(lam):
    return lam / (1 - math.exp(-lam)) if lam > 0 else 1.0


def generic_excess(D, clean, test, min_excess=100):
    """Foreign k-mers above the clean-line ONT null and their dosage (pooled queries).

    Null = median normalized recovery / dosage of the same class in the clean
    lines' pooled ONT. Dosage range: null hits at the null dosage (low) or at
    abundance 1 (high). Whole-genome admixture at fraction f gives dosage ~f."""
    rows = []
    for k in (21, 31):
        C = rd(os.path.join(D, f"k{k}", "class_recovery.tsv"))
        idx = {(r["query"], r["class"]): r for r in C}
        for g in test:
            q = f"{g}:all_ont_runs"
            core = idx[(q, "core")]
            s, lc = float(core["sensitivity"]), float(core["lambda_core"])
            for c in ["foreign_all"] + [f"foreign_in{i}" for i in range(1, 9)]:
                r = idx[(q, c)]
                n, h, m = int(r["n_class"]), int(r["hits"]), float(r["mean_abund_hits"])
                nr = statistics.median(float(idx[(f"{x}:all_ont_runs", c)]["norm_recovery"])
                                       for x in clean)
                nd = statistics.median(float(idx[(f"{x}:all_ont_runs", c)]["dose_vs_core"])
                                       for x in clean)
                hn = n * nr * s
                ex = h - hn
                lo = hi = None
                if ex >= min_excess:
                    lo = zt_lambda((h * m - hn * zt_mean(nd * lc)) / ex) / lc
                    hi = zt_lambda((h * m - hn) / ex) / lc
                rows.append({"k": k, "query": q, "class": c, "sensitivity": s, "lambda_core": lc,
                             "own_single_dose": float(idx[(q, "own_single(0.65-1.5peak)")]
                                                      ["dose_vs_core"]),
                             "n_class": n, "hits": h, "null_hits": hn, "excess_hits": ex,
                             "excess_fraction": ex / h if h else None,
                             "excess_dose_low": lo, "excess_dose_high": hi})
    return rows


def detection_limits(report, D, test, exclude):
    """Indicative smallest whole-genome admixture fraction f of line X detectable in
    each pooled ONT query: expected private score (null median + recovery(f)/s)
    exceeds the clean-line null max by 3 binomial SD; recovery(f) = 1 - exp(-f *
    lambda), lambda = core depth x own-private dosage. Assumes independent k-mers
    (anti-conservative)."""
    S = rd(os.path.join(report, "k21", "scores_long.tsv"))
    C = rd(os.path.join(D, "k21", "class_recovery.tsv"))
    lam = {r["query"]: float(r["lambda_core"]) for r in C if r["class"] == "core"}
    own = {r["query"]: float(r["dose_vs_core"]) for r in C if r["class"] == "own_private"
           and r["dose_vs_core"] != "NA"}
    null = [float(r["score"]) for r in S if r["source"] == "ont" and r["level"] == "pooled"
            and r["expected"] not in test and r["genotype"] != r["expected"]
            and r["query"].split(":")[0] not in exclude]
    nmed, nmax = statistics.median(null), max(null)
    gs = sorted({r["genotype"] for r in S})
    npriv = {r["genotype"]: int(r["n_private"]) for r in S}
    rows = []
    for q in sorted({r["query"] for r in S if r["source"] == "ont" and r["level"] == "pooled"}):
        exp = next(r["expected"] for r in S if r["query"] == q)
        s = next(float(r["sensitivity"]) for r in S if r["query"] == q)
        L = lam[q] * own[q]
        row = {"query": q, "sensitivity": s, "lambda_private": L,
               "null_median": nmed, "null_max": nmax}
        for X in gs:
            if X == exp:
                row[X] = None
                continue
            f = None
            for i in range(1, 1001):
                ff = i / 1000
                p = 1 - math.exp(-ff * L)
                pt = min(0.999, p + nmed * s)
                if nmed + p / s >= nmax + 3 * math.sqrt(pt * (1 - pt) / npriv[X]) / s:
                    f = ff
                    break
            row[X] = f
        rows.append(row)
    return rows, gs


def foreign_rows(report, k):
    """ONT queries: foreign score of every other genotype H vs the expected G."""
    sens = {}
    for r in rd(os.path.join(report, f"k{k}", "scores_long.tsv")):
        sens[r["query"]] = float(r["sensitivity"])
    out = []
    for r in rd(os.path.join(report, f"k{k}", "pairwise_long.tsv")):
        if r["source"] != "ont" or r["versus"] != r["expected"]:
            continue
        out.append({"query": r["query"], "level": r["level"], "expected": r["expected"],
                    "other": r["genotype"], "n_distinct": int(r["n_distinct"]),
                    "hits": int(r["hits"]), "score": float(r["score"]),
                    "sensitivity": sens[r["query"]]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True)
    ap.add_argument("--diagnostics", required=True)
    ap.add_argument("--test", required=True,
                    help="comma-separated genotypes under investigation (excluded from the null)")
    ap.add_argument("--exclude-queries", default="",
                    help="comma-separated run ids with insufficient data (excluded from the null)")
    a = ap.parse_args()
    test = set(a.test.split(","))
    excl = set(a.exclude_queries.split(","))
    D = a.diagnostics

    def is_excl(q):
        return q.split(":")[0] in excl

    # ---- ONT null and foreign scores at k = 21 and 31 -----------------------
    calib, fk, null_max = [], {}, {}
    for k in (21, 31):
        rows = foreign_rows(a.report, k)
        null = [r for r in rows if r["expected"] not in test and not is_excl(r["query"])
                and r["level"] in ("run", "half")]
        for lo, hi in BINS + [(0.0, 1.01)]:
            xs = [r["score"] for r in null if lo <= r["sensitivity"] < hi]
            lab = f"{lo:.2f}-{min(hi, 1):.2f}" if (lo, hi) != (0.0, 1.01) else "all"
            if not xs:
                continue
            row = {"k": k, "sensitivity_bin": lab, "n_scores": len(xs),
                   "n_queries": len({r["query"] for r in null if lo <= r["sensitivity"] < hi}),
                   "median": statistics.median(xs), "p95": pct(xs, 0.95),
                   "p99": pct(xs, 0.99), "max": max(xs)}
            calib.append(row)
            null_max[(k, lab)] = row["max"]
        pooled = [r["score"] for r in rows if r["expected"] not in test and r["level"] == "pooled"
                  and not is_excl(r["query"])]
        null_max[(k, "pooled_queries")] = max(pooled)
        calib.append({"k": k, "sensitivity_bin": "pooled_queries", "n_scores": len(pooled),
                      "n_queries": len({r["query"] for r in rows if r["level"] == "pooled"
                                        and r["expected"] not in test}),
                      "median": statistics.median(pooled), "p95": pct(pooled, 0.95),
                      "p99": pct(pooled, 0.99), "max": max(pooled)})
        for r in rows:
            fk[(r["query"], r["other"], k)] = r
    wt(os.path.join(D, "ont_null_calibration.tsv"), calib)

    out = []
    for (q, other, k), r in fk.items():
        if k != 21:
            continue
        r31 = fk.get((q, other, 31))
        # pooled queries sum several runs (more error k-mers): compare with pooled nulls
        pooled_q = r["level"] == "pooled"
        nmax21 = null_max.get((21, "pooled_queries" if pooled_q else binof(r["sensitivity"])))
        nmax31 = (null_max.get((31, "pooled_queries" if pooled_q else binof(r31["sensitivity"])))
                  if r31 else None)
        out.append({"query": q, "level": r["level"], "expected": r["expected"], "other": other,
                    "sensitivity": r["sensitivity"],
                    "hits_k21": r["hits"], "n_k21": r["n_distinct"], "score_k21": r["score"],
                    "ont_null_max_k21": nmax21,
                    "excess_k21": (r["score"] - nmax21) if nmax21 is not None else None,
                    "hits_k31": r31["hits"] if r31 else None,
                    "n_k31": r31["n_distinct"] if r31 else None,
                    "score_k31": r31["score"] if r31 else None,
                    "ont_null_max_k31": nmax31,
                    "excess_k31": (r31["score"] - nmax31) if (r31 and nmax31 is not None) else None})
    out.sort(key=lambda x: (x["expected"], x["level"], x["query"], x["other"]))
    wt(os.path.join(D, "foreign_scores_by_k.tsv"), out)

    # ---- has / lacks vs sensitivity and sequencing descriptors -----------------
    recs = {}
    with open(os.path.join(a.report, "sketch_records.txt")) as fh:
        for p in (x.strip() for x in fh if x.strip()):
            with open(p) as g:
                j = json.load(g)
            if j["tech"] == "ont":
                recs[j["run_id"]] = j["sampling"]
    rs = []
    for k in (21, 31):
        for r in rd(os.path.join(a.report, f"k{k}", "assignments.tsv")):
            run = r["query"].split(":")[0]
            smp = recs.get(run, {})
            rs.append({"k": k, "query": r["query"], "level": r["level"], "expected": r["expected"],
                       "class": r["classification"], "sensitivity": float(r["sensitivity"]),
                       "has": float(r["expected_has_min"]), "lacks": float(r["expected_lacks_max"]),
                       "bases_used": smp.get("bases_used") if r["level"] == "run" else None,
                       "read_len_n50": smp.get("read_len_n50") if r["level"] != "pooled" else None})
    wt(os.path.join(D, "run_sensitivity.tsv"), rs)

    # ---- read-level co-occurrence ---------------------------------------------
    rc = []
    for f in sorted(glob.glob(os.path.join(D, "k21", "reads_*.tsv"))):
        R = rd(f)
        lab = os.path.basename(f).split("_")[1]
        run = os.path.basename(f)[len("reads_" + lab + "_"):-4]
        for r in R:
            for c in r:
                if c != "mean_q":
                    r[c] = int(r[c])
            r["mean_q"] = float(r["mean_q"])
        gs = [c[2:] for c in R[0] if c.startswith("f_")]
        tf = sum(r["n_foreign"] for r in R)
        to = sum(r["n_own"] for r in R)
        obs2 = exp2 = obsP = expP = 0.0
        n2 = 0
        for lo, hi in ((0, 12), (12, 15), (15, 18), (18, 99)):       # quality strata
            S = [r for r in R if lo <= r["mean_q"] < hi]
            sf, so = sum(r["n_foreign"] for r in S), sum(r["n_own"] for r in S)
            if not sf + so:
                continue
            p = sf / (sf + so)
            for r in S:
                n = r["n_foreign"] + r["n_own"]
                if n < 2:
                    continue
                n2 += 1
                obs2 += r["n_foreign"] >= 2
                exp2 += sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(2, n + 1))
                obsP += r["n_foreign"] == n
                expP += p ** n
        fr = [r for r in R if r["n_foreign"]]
        fl = sorted((r["n_foreign"] for r in fr), reverse=True)
        multi = [r for r in R if r["n_foreign"] >= 2]
        row = {"expected": lab, "run": run, "reads_with_informative": len(R),
               "own_markers": to, "foreign_markers": tf, "foreign_fraction": tf / (tf + to),
               "reads_with_foreign": len(fr),
               "share_top1pct_reads": sum(fl[:max(1, len(fl) // 100)]) / tf,
               "max_foreign_per_read": fl[0] if fl else 0,
               "reads_ge2_informative": n2,
               "ge2_foreign_obs": int(obs2), "ge2_foreign_exp": exp2,
               "ge2_ratio": obs2 / exp2 if exp2 else None,
               "all_foreign_obs": int(obsP), "all_foreign_exp": expP,
               "all_foreign_ratio": obsP / expP if expP else None,
               "medianQ_foreign_reads": statistics.median(r["mean_q"] for r in fr),
               "medianQ_own_only_reads": statistics.median(r["mean_q"] for r in R
                                                           if r["n_own"] and not r["n_foreign"])}
        for g in gs:
            pg = sum(r["f_" + g] for r in R) / tf
            row[f"multi_all_{g}_obs"] = sum(1 for r in multi if r["f_" + g] == r["n_foreign"])
            row[f"multi_all_{g}_exp"] = sum(pg ** r["n_foreign"] for r in multi)
        rc.append(row)
    if rc:
        cols = list(dict.fromkeys(c for r in rc for c in r))
        wt(os.path.join(D, "read_cooccurrence.tsv"), rc, cols)

    clean = sorted({r["expected"] for r in rd(os.path.join(a.report, "k21", "assignments.tsv"))
                    if r["level"] == "pooled" and r["expected"] not in test
                    and float(r["sensitivity"]) >= 0.9})
    wt(os.path.join(D, "generic_foreign_excess.tsv"), generic_excess(D, clean, sorted(test)))
    dl, gs = detection_limits(a.report, D, test, excl)
    wt(os.path.join(D, "detection_limits_k21.tsv"), dl,
       ["query", "sensitivity", "lambda_private", "null_median", "null_max"] + gs)


if __name__ == "__main__":
    main()
