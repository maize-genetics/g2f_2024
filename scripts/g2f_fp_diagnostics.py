#!/usr/bin/env python3
"""Supplementary, read-only diagnostics for the genotype fingerprinting report.

Explains elevated "foreign allele" (lacks) scores of ONT read sets using only
the existing sketches (no reads are read). It rebuilds the reference panel
with the unchanged production code (g2f.fingerprint) and the settings of an
existing report, then, for every ONT query and its expected genotype G:

  classes   recovery of k-mer classes relative to G: core, G's informative
            k-mers stratified by G's HiFi abundance (heterozygous-like vs
            single-copy), and "foreign" k-mers (solid in another genotype,
            never seen in G's HiFi reads) stratified by how many panel
            genotypes carry them, plus per-genotype private k-mers.
  dose      zero-truncated-Poisson depth of the recovered hashes of each class
            relative to core k-mers (pooled queries): ~1 = same dosage as the
            genome itself; << 1 = a minor component (low-level mixture).
  replication  of each class between independent runs / replicate halves,
            relative to core k-mers (sequencing errors do not replicate).
  null      ONT pairwise null from the genotypes without elevated scores, as a
            function of sensitivity, to compare like with like.

Nothing in the production pipeline or its outputs is modified. Outputs go
to a new directory. Standard library only (plus g2f.fingerprint).

  python3.12 scripts/g2f_fp_diagnostics.py \
      --report ../results/fingerprint/report-742a05a563 \
      --outdir ../results/fingerprint/diagnostics-2026-10-12 --k 21
"""
import argparse
import json
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2f import fingerprint as fp          # noqa: E402
from g2f import fpanalysis                  # noqa: E402


# --------------------------------------------------------------------------

def zt_lambda(mean):
    """Poisson mean whose zero-truncated mean is `mean` (>= 1)."""
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


def write_tsv(path, rows, cols=None):
    cols = cols or list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fpanalysis._fmt(r.get(c)) for c in cols) + "\n")


def load_records(report):
    with open(os.path.join(report, "sketch_records.txt")) as fh:
        paths = [x.strip() for x in fh if x.strip()]
    recs = fpanalysis._load_records(paths)
    with open(os.path.join(report, "parameters_and_provenance.json")) as fh:
        prov = json.load(fh)
    return recs, prov


def params_of(prov):
    for key in ("params", "analysis_params", "parameters"):
        if isinstance(prov.get(key), dict):
            return dict(fpanalysis.DEFAULTS, **prov[key])
    return dict(fpanalysis.DEFAULTS)


# --------------------------------------------------------------------------

def build(records, k, P, log):
    hifi = [r for r in records if r["tech"] == "hifi"]
    by_g = {}
    for r in hifi:
        by_g.setdefault(r["label"], []).append(r)
    refs, counts_of = {}, {}
    for g in sorted(by_g):
        counts = {}
        for r in by_g[g]:
            for d in fp.load_sig(r["signature"], k).values():
                fp.add_counts(counts, d)
        bases = sum(r["sampling"]["bases_used"] for r in by_g[g])
        refs[g] = fp.Reference(g, counts, P, [r["run_id"] for r in by_g[g]], bases=bases)
        counts_of[g] = counts
        log(f"k={k} {g}: band {refs[g].band} solid {len(refs[g].solid)}")
    panel = fp.Panel(refs)
    solid_n = {}
    for g in panel.genotypes:
        for h in refs[g].solid:
            solid_n[h] = solid_n.get(h, 0) + 1
    return panel, counts_of, solid_n


def ont_queries(records, k):
    """[(query, level, label, counts dict)] for runs, halves and pooled labels."""
    out, pooled = [], {}
    for r in sorted((r for r in records if r["tech"] == "ont"),
                    key=lambda x: (x["label"], x["run_id"])):
        halves = fp.load_sig(r["signature"], k)
        full = {}
        for d in halves.values():
            fp.add_counts(full, d)
        out.append((r["run_id"], "run", r["label"], full, r))
        for h in sorted(halves):
            out.append((f"{r['run_id']}:half{h}", "half", r["label"], halves[h], r))
        fp.add_counts(pooled.setdefault(r["label"], {}), full)
    for lab in sorted(pooled):
        out.append((f"{lab}:all_ont_runs", "pooled", lab, pooled[lab], None))
    return out


def marker_classes(panel, counts_of, solid_n, g):
    """Hash classes relative to genotype g (dict name -> set)."""
    ref = panel.refs[g]
    bit, full = panel.bit[g], sum(panel.bit.values())
    peak = ref.band["peak"]
    cls = {"core": panel.core}
    own = {h for h in ref.solid if panel.mask[h] != full}
    cls["own_informative"] = own
    cls["own_private"] = panel.private[g]
    ab = counts_of[g]
    cls["own_lowdose(<0.65peak)"] = {h for h in own if ab[h] < 0.65 * peak}
    cls["own_single(0.65-1.5peak)"] = {h for h in own if 0.65 * peak <= ab[h] <= 1.5 * peak}
    cls["own_multi(>1.5peak)"] = {h for h in own if ab[h] > 1.5 * peak}
    foreign = {h for h in panel.universe if not panel.mask[h] & bit}
    cls["foreign_all"] = foreign
    for c in range(1, len(panel.genotypes)):
        cls[f"foreign_in{c}"] = {h for h in foreign if solid_n[h] == c}
    for h2 in panel.genotypes:
        if h2 != g:
            cls[f"private:{h2}"] = panel.private[h2]
            cls[f"foreign_solid_in:{h2}"] = foreign & panel.refs[h2].solid
    return cls


def class_rows(qname, level, label, q, classes, n_core):
    core_hits = [q[h] for h in classes["core"] if h in q]
    s = len(core_hits) / n_core if n_core else 0.0
    lam_core = zt_lambda(statistics.fmean(core_hits)) if core_hits else 0.0
    rows = []
    for name, hs in classes.items():
        ab = [q[h] for h in hs if h in q]
        m = statistics.fmean(ab) if ab else None
        lam = zt_lambda(m) if m else 0.0
        rows.append({"query": qname, "level": level, "expected": label, "class": name,
                     "n_class": len(hs), "hits": len(ab),
                     "recovery": len(ab) / len(hs) if hs else None,
                     "sensitivity": s,
                     "norm_recovery": (len(ab) / len(hs) / s) if (hs and s) else None,
                     "mean_abund_hits": m,
                     "frac_hits_abund_ge2": (sum(1 for a in ab if a >= 2) / len(ab)) if ab else None,
                     "lambda_zt": lam, "lambda_core": lam_core,
                     "dose_vs_core": (lam / lam_core) if lam_core else None})
    return rows


def replication(qa, qb, classes, any_present):
    """For hashes of a class found in qa: fraction also found in qb, / same for core."""
    sa, sb = set(qa), set(qb)
    core_a = classes["core"] & sa
    rc = len(core_a & sb) / len(core_a) if core_a else None
    out = {"core_rep": rc}
    for name in ("own_informative", "own_lowdose(<0.65peak)", "foreign_all", "foreign_in1",
                 "foreign_in2", "foreign_in4", "foreign_in6", "foreign_in8"):
        if name not in classes:
            continue
        a = classes[name] & sa
        r = len(a & sb) / len(a) if a else None
        out[name] = r
        out[name + "_rel"] = (r / rc) if (r is not None and rc) else None
        out[name + "_n"] = len(a)
    novel = sa - any_present
    out["novel_n"] = len(novel)
    out["novel_rep"] = len(novel & sb) / len(novel) if novel else None
    out["novel_rep_rel"] = (out["novel_rep"] / rc) if (out["novel_rep"] is not None and rc) else None
    return out


def hifi_spectra(panel, counts_of, maxa=40):
    rows = []
    for g in panel.genotypes:
        hist = fp.spectrum(counts_of[g])
        peak = panel.refs[g].band["peak"]
        tot = sum(v for a, v in hist.items() if a >= 2)
        row = {"genotype": g, "peak": peak}
        for a in range(1, maxa + 1):
            row[str(a)] = hist.get(a, 0)
        lo_half = [a for a in range(max(2, round(0.35 * peak)), round(0.65 * peak) + 1)]
        main = [a for a in range(round(0.8 * peak), round(1.3 * peak) + 1)]
        row["half_peak_mass"] = sum(hist.get(a, 0) for a in lo_half)
        row["main_peak_mass"] = sum(hist.get(a, 0) for a in main)
        row["half_over_main"] = row["half_peak_mass"] / row["main_peak_mass"]
        row["n_ge2"] = tot
        rows.append(row)
    return rows


def lowdose_composition(panel, counts_of, g):
    """For G's informative k-mers by HiFi abundance bin: fraction present in each other line."""
    ref, ab = panel.refs[g], counts_of[g]
    peak, full = ref.band["peak"], sum(panel.bit.values())
    bins = {"low(<0.65peak)": lambda a: a < 0.65 * peak,
            "single(0.65-1.5peak)": lambda a: 0.65 * peak <= a <= 1.5 * peak}
    rows = []
    for bname, f in bins.items():
        hs = [h for h in ref.solid if panel.mask[h] != full and f(ab[h])]
        row = {"genotype": g, "bin": bname, "n": len(hs),
               "mean_lines_present": statistics.fmean(
                   bin(panel.mask[h]).count("1") - 1 for h in hs) if hs else None}
        for h2 in panel.genotypes:
            if h2 != g:
                b = panel.bit[h2]
                row[h2] = sum(1 for h in hs if panel.mask[h] & b) / len(hs) if hs else None
        rows.append(row)
    return rows


def line_strata(panel, solid_n, classes, label, max_c=3):
    """{(H, c): (foreign k-mers of `label` carried by c genotypes incl. H, ... excl. H)}."""
    foreign = classes["foreign_all"]
    by_c = {c: [h for h in foreign if solid_n[h] == c] for c in range(1, max_c + 1)}
    out = {}
    for h2 in panel.genotypes:
        if h2 != label:
            sol = panel.refs[h2].solid
            for c, hs in by_c.items():
                out[(h2, c)] = ([h for h in hs if h in sol], [h for h in hs if h not in sol])
    return out


def line_specific(qname, level, label, q, panel, strata, max_c=3):
    """Excess of genotype H's alleles beyond what allele sharing explains.

    Within foreign k-mers carried by c (<= max_c) panel genotypes, compare
    those carried by H with those not carried by H. Under pure allele-frequency
    effects both recover at the same rate; contamination by (or a haplotype
    shared with) H gives an excess. The excess hits' depth (zero-truncated
    Poisson, by subtraction) relative to core gives its dosage: ~1 genomic,
    << 1 a minor DNA component."""
    core_ab = [q[h] for h in panel.core if h in q]
    s = len(core_ab) / len(panel.core)
    lam_core = zt_lambda(statistics.fmean(core_ab)) if core_ab else 0.0
    rows = []
    for h2 in panel.genotypes:
        if h2 == label:
            continue
        n_w = hit_w = ab_w = n_o = hit_o = ab_o = 0
        exp_hits = exp_ab = 0.0
        for c in range(1, max_c + 1):
            w, o = strata[(h2, c)]
            hw = [q[h] for h in w if h in q]
            ho = [q[h] for h in o if h in q]
            rate_o = len(ho) / len(o) if o else 0.0
            mean_o = statistics.fmean(ho) if ho else 0.0
            n_w += len(w); hit_w += len(hw); ab_w += sum(hw)
            n_o += len(o); hit_o += len(ho); ab_o += sum(ho)
            exp_hits += len(w) * rate_o
            exp_ab += len(w) * rate_o * mean_o
        ex_hits = hit_w - exp_hits
        ex_mean = (ab_w - exp_ab) / ex_hits if ex_hits > 0 else None
        lam = zt_lambda(ex_mean) if ex_mean else None
        rows.append({"query": qname, "level": level, "expected": label, "versus": h2,
                     "max_sharing": max_c, "sensitivity": s,
                     "n_with": n_w, "hits_with": hit_w,
                     "norm_rec_with": hit_w / n_w / s if (n_w and s) else None,
                     "n_without": n_o, "hits_without": hit_o,
                     "norm_rec_without": hit_o / n_o / s if (n_o and s) else None,
                     "expected_hits_with": exp_hits, "excess_hits": ex_hits,
                     "excess_rel": ex_hits / exp_hits if exp_hits else None,
                     "excess_mean_abund": ex_mean, "lambda_core": lam_core,
                     # dosage only when the excess is large enough to estimate
                     "excess_dose_vs_core": (lam / lam_core) if (lam is not None and lam_core
                                                                   and ex_hits >= 100) else None})
    return rows


def spikein(panel, solid_n, classes, queries, specs, fractions, calib, P, seed, log):
    """In-silico admixture of real ONT sketches (positive controls for contamination).

    host and donor are pooled ONT sketches of two genotypes. The donor's
    abundances are binomially thinned so that it makes up a fraction f of the
    mixture's genome depth (core-k-mer depth ratio), then added to the host.
    Each mixture is scored with the production scoring/classification code and
    the line-specific excess test."""
    import random
    rng = random.Random(seed)
    pooled = {lab: q for qn, lev, lab, q, _ in queries if lev == "pooled"}

    def lam(q):
        ab = [q[h] for h in panel.core if h in q]
        return zt_lambda(statistics.fmean(ab)) if ab else 0.0

    rows = []
    for spec in specs:
        host, donor = spec.split(":")
        hq, dq = pooled[host], pooled[donor]
        lh, ld = lam(hq), lam(dq)
        strata = line_strata(panel, solid_n, classes[host], host)
        for f in fractions:
            t = min(1.0, f / (1 - f) * lh / ld) if f > 0 else 0.0
            mix = dict(hq)
            for h, a in dq.items():
                k = sum(1 for _ in range(a) if rng.random() < t)
                if k:
                    mix[h] = mix.get(h, 0) + k
            res = fp.score(panel, set(mix) & panel.universe, P["interval_z"])
            asg = fp.classify(res, host, calib, P["min_informative_hashes"])
            sp = fp.support(res, host, P["min_informative_hashes"])
            pr = res["pairs"][(donor, host)]
            ls = [r for r in line_specific("mix", "spikein", host, mix, panel, strata)
                  if r["versus"] == donor][0]
            others = sorted((v["score"] for (g, b), v in res["pairs"].items()
                             if b == host and g != donor))
            rows.append({"host": host, "donor": donor, "fraction": f, "thinning": t,
                         "sensitivity": res["sensitivity"], "has_min": sp["has_min"],
                         "lacks_max": sp["lacks_max"], "lacks_max_line": sp["strongest_other"],
                         "donor_pair_score": pr["score"], "donor_pair_lo": pr["lo"],
                         "other_lines_median": statistics.median(others) if others else None,
                         "donor_private_score": res["genotypes"][donor]["score"],
                         "class": asg["classification"],
                         "possible_admixture": asg["possible_admixture"],
                         "donor_excess_rel": ls["excess_rel"],
                         "donor_excess_dose": ls["excess_dose_vs_core"]})
            log(f"spike-in {host}+{donor} f={f}: lacks {sp['lacks_max']:.3f} "
                f"({sp['strongest_other']}), class {asg['classification']}")
    return rows


# --------------------------------------------------------------------------
# Read-level co-occurrence (streams reads; bounded)
# --------------------------------------------------------------------------

def _fastq_with_qual(lines):
    it = iter(lines)
    for header in it:
        if not header.startswith("@"):
            continue
        seq = next(it).rstrip("\n")
        next(it)
        qual = next(it).rstrip("\n")
        yield header[1:].split(None, 1)[0], seq, qual


def reads_scan(path, label, panel, classes, k, scaled, max_bases, min_len, out_tsv, log):
    """Per-read counts of informative sampled k-mers relative to `label`."""
    from sourmash import MinHash
    from g2f import archive
    gs = [g for g in panel.genotypes if g != label]
    own = classes["own_informative"]
    low = classes["own_lowdose(<0.65peak)"]
    foreign = classes["foreign_all"]
    core = panel.core
    solid = {g: panel.refs[g].solid for g in gs}
    lines = archive.open_text(path) if archive.is_tar_name(os.path.basename(path).lower()) \
        else open(path)
    mh = MinHash(n=0, ksize=k, scaled=scaled, seed=42)
    used = n_reads = 0
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(["read", "length", "mean_q", "n_hashes", "n_core", "n_own",
                            "n_own_low", "n_foreign"] + [f"f_{g}" for g in gs]) + "\n")
        for name, seq, qual in _fastq_with_qual(lines):
            if len(seq) < min_len:
                continue
            n_reads += 1
            used += len(seq)
            mh.clear()
            mh.add_sequence(seq, force=True)
            hs = mh.hashes
            nf = [h for h in hs if h in foreign]
            no = sum(1 for h in hs if h in own)
            if nf or no:
                # mean Phred from ASCII (error-probability mean, reported as Q)
                pe = sum(10 ** (-(ord(c) - 33) / 10) for c in qual[::10]) / max(1, len(qual[::10]))
                mq = -10 * math.log10(pe) if pe > 0 else 60
                fh.write("\t".join(map(str, [
                    n_reads, len(seq), round(mq, 1), len(hs),
                    sum(1 for h in hs if h in core), no,
                    sum(1 for h in hs if h in low), len(nf)]
                    + [sum(1 for h in nf if h in solid[g]) for g in gs])) + "\n")
            if max_bases and used >= max_bases:
                break
    log(f"{path}: {n_reads} reads, {used} bases scanned")
    return {"reads": n_reads, "bases": used}


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="existing fingerprint report directory")
    ap.add_argument("--outdir", required=True, help="new output directory")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--reads", nargs="*", metavar="PATH:LABEL",
                    help="read files to scan per read (bounded by --reads-max-bases)")
    ap.add_argument("--reads-max-bases", type=int, default=1_500_000_000)
    ap.add_argument("--spikein", nargs="*", metavar="HOST:DONOR",
                    help="in-silico admixture of pooled ONT sketches (then exit)")
    ap.add_argument("--spikein-fractions", default="0,0.01,0.02,0.05,0.1,0.2,0.5")
    ap.add_argument("--only-reads", action="store_true",
                    help="skip the sketch-level tables (already written)")
    a = ap.parse_args()
    out = os.path.join(a.outdir, f"k{a.k}")
    os.makedirs(out, exist_ok=True)
    log = lambda m: (sys.stderr.write(m + "\n"), sys.stderr.flush())  # noqa: E731

    records, prov = load_records(a.report)
    P = params_of(prov)
    log(f"params: {P}")
    panel, counts_of, solid_n = build(records, a.k, P, log)
    any_present = set()
    for g in panel.genotypes:
        any_present |= panel.refs[g].present
    n_core = len(panel.core)
    queries = ont_queries(records, a.k)
    classes = {g: marker_classes(panel, counts_of, solid_n, g) for g in panel.genotypes}

    if a.spikein:
        with open(os.path.join(a.report, f"k{a.k}", "calibration.json")) as fh:
            calib = json.load(fh)
        fr = [float(x) for x in a.spikein_fractions.split(",")]
        rows = spikein(panel, solid_n, classes, queries, a.spikein, fr, calib, P, 1, log)
        write_tsv(os.path.join(out, "spikein_admixture.tsv"), rows)
        return

    if a.only_reads:
        for spec in a.reads or []:
            path, label = spec.rsplit(":", 1)
            tag = os.path.basename(path).split(".")[0]
            reads_scan(path, label, panel, classes[label], a.k, 2000, a.reads_max_bases, 1000,
                       os.path.join(out, f"reads_{label}_{tag}.tsv"), log)
        return

    # class recovery / dose for every ONT query vs its expected genotype
    rows = []
    for qn, lev, lab, q, _ in queries:
        if lab in classes:
            rows += class_rows(qn, lev, lab, q, classes[lab], n_core)
    write_tsv(os.path.join(out, "class_recovery.tsv"), rows)

    # genotype-specific excess beyond allele sharing (runs and pooled)
    ls = []
    for g in panel.genotypes:
        strata = line_strata(panel, solid_n, classes[g], g)
        for qn, lev, lab, q, _ in queries:
            if lab == g and lev in ("run", "pooled") and len(q) > 1000:
                ls += line_specific(qn, lev, lab, q, panel, strata)
        log(f"line-specific excess: {g} done")
    write_tsv(os.path.join(out, "line_specific_excess.tsv"), ls)

    # optional read-level co-occurrence scans (bounded prefix of each file)
    for spec in a.reads or []:
        path, label = spec.rsplit(":", 1)
        tag = os.path.basename(path).split(".")[0]
        reads_scan(path, label, panel, classes[label], a.k, 2000, a.reads_max_bases, 1000,
                   os.path.join(out, f"reads_{label}_{tag}.tsv"), log)

    # replication between runs of the same label, and between halves of each run
    rep = []
    by_lab = {}
    for qn, lev, lab, q, r in queries:
        by_lab.setdefault(lab, []).append((qn, lev, q))
    for lab, items in by_lab.items():
        runs = [(qn, q) for qn, lev, q in items if lev == "run"]
        halves = [(qn, q) for qn, lev, q in items if lev == "half"]
        pairs = [(x, y, "between_runs") for i, x in enumerate(runs) for y in runs[i + 1:]]
        pairs += [(x, y, "between_runs") for i, y in enumerate(runs) for x in runs[i + 1:]]
        hb = {}
        for qn, q in halves:
            hb.setdefault(qn.rsplit(":", 1)[0], []).append((qn, q))
        for run, hs in hb.items():
            if len(hs) == 2:
                pairs.append((hs[0], hs[1], "between_halves"))
        for (na, qa), (nb, qb), kind in pairs:
            if len(qa) < 1000 or len(qb) < 1000:
                continue
            row = {"expected": lab, "kind": kind, "a": na, "b": nb}
            row.update(replication(qa, qb, classes[lab], any_present))
            rep.append(row)
    write_tsv(os.path.join(out, "replication.tsv"), rep)

    write_tsv(os.path.join(out, "hifi_spectra.tsv"), hifi_spectra(panel, counts_of))
    comp = []
    for g in panel.genotypes:
        comp += lowdose_composition(panel, counts_of, g)
    write_tsv(os.path.join(out, "hifi_informative_by_dose.tsv"), comp,
              ["genotype", "bin", "n", "mean_lines_present"] + panel.genotypes)

    # sizes of the foreign sets relative to each G, by sharing count (denominators)
    comp_rows = []
    for g in panel.genotypes:
        for h2 in panel.genotypes:
            if h2 == g:
                continue
            d = classes[g][f"foreign_solid_in:{h2}"]
            row = {"genotype": g, "versus": h2, "n_distinct": len(d),
                   "mean_lines_solid": statistics.fmean(solid_n[h] for h in d) if d else None}
            for c in range(1, len(panel.genotypes)):
                row[f"frac_in{c}"] = sum(1 for h in d if solid_n[h] == c) / len(d) if d else None
            comp_rows.append(row)
    write_tsv(os.path.join(out, "distinct_set_composition.tsv"), comp_rows)
    with open(os.path.join(out, "panel_info.json"), "w") as fh:
        json.dump({"k": a.k, "n_core": n_core, "params": P,
                   "bands": {g: panel.refs[g].band for g in panel.genotypes},
                   "report": os.path.abspath(a.report)}, fh, indent=2)
    log("done")


if __name__ == "__main__":
    main()
