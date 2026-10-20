"""Classify read-support windows and clip clusters into candidate regions.

Every call is a candidate with its evidence and a confidence level, never a
confirmed error. Read depth is expressed relative to one depth baseline c per
sample: the depth of sequence represented once in the assembly that receives the
reads of both haplotypes (homozygous regions, or heterozygous regions collapsed
to one allele). This is hifiasm's peak_hom / the homozygous k-mer peak, the same
for every assembly of the sample. It is NOT the modal depth of each assembly:
in an unpurged assembly of a heterozygous sample both alleles are kept, each
receives about c/2, and the depth mode is then the one-allele peak.

  none / low   (< 0.15c / < 0.30c)  weak support: possible erroneous sequence,
                                     or a coverage bias
  half         (0.30-0.70c)         one haplotype's reads: a heterozygous region
                                     present twice (both alleles kept), or one
                                     copy of a near-identical duplicate
                                     (multi-mapping reads split between copies)
  normal       (0.70-1.40c)
  dup          (1.40-2.50c)         collapsed repeat / collapsed duplicate
  multi        (> 2.50c)            high-copy collapse or repeat

Supporting evidence used to grade confidence:
  - mixed-site density (NucFreq-style): elevated in collapses of diverged copies
    and in heterozygous regions represented once
  - multi-mapping fraction (MAPQ < min_mapq): near-identical copies in the assembly
  - read continuity (span_min) and clipped-read clusters for breakpoints
"""
import gzip
import json

import numpy as np

CLASSES = ["none", "low", "half", "normal", "dup", "multi"]


class CalibrationError(ValueError):
    """The depth baseline cannot be resolved safely; the stage must not classify."""


def expected_depth(setting, hist_mode, kmer_peak=None, tolerance=0.15):
    """Resolve the depth baseline c; returns (value, source, warning).

    A number is used as given (warning if it differs from the k-mer peak).
    "auto" takes the read-depth histogram mode only after checking it against
    the k-mer peak; it raises CalibrationError when no k-mer peak is given or the
    mode differs from it by more than `tolerance` (e.g. a mode at half the k-mer
    peak in an assembly that keeps both alleles). Classifying against an
    unverified baseline would mislabel whole depth classes.
    """
    warn = ""
    if setting not in (None, "auto"):
        c = float(setting)
        if c <= 0:
            raise CalibrationError(f"expected_depth must be positive, got {setting!r}")
        if kmer_peak and abs(c - kmer_peak) / kmer_peak > tolerance:
            warn = f"configured expected_depth {c:g} differs from k-mer peak {kmer_peak} by more than {tolerance:.0%}"
        return c, "configured", warn
    c = float(hist_mode)
    if not kmer_peak:
        raise CalibrationError(
            f"expected_depth is 'auto' but no kmer_peak is configured to verify the read-depth mode ({c:g}); "
            "set kmer_peak (hifiasm peak_hom or the homozygous k-mer peak) or a numeric expected_depth")
    if abs(c - kmer_peak) / kmer_peak > tolerance:
        hint = (" The mode is close to half the k-mer peak: the assembly probably keeps both alleles of "
                "heterozygous regions, so the mode is the one-allele depth." if abs(c / kmer_peak - 0.5) <= tolerance / 2
                else "")
        raise CalibrationError(
            f"read-depth mode {c:g} differs from k-mer peak {kmer_peak} by more than {tolerance:.0%}.{hint} "
            f"Set expected_depth to the shared homozygous depth (normally {kmer_peak}) for every assembly of "
            "the sample, or correct kmer_peak.")
    return c, "depth-histogram mode (checked against k-mer peak)", warn


def depth_class(ratio, b):
    if ratio < b["none"]:
        return "none"
    if ratio < b["low"]:
        return "low"
    if ratio < b["half_max"]:
        return "half"
    if ratio < b["normal_max"]:
        return "normal"
    if ratio < b["dup_max"]:
        return "dup"
    return "multi"


def read_windows(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    cols = {h: i for i, h in enumerate(header)}
    w = {"contig": np.array([r[0] for r in rows], dtype=object)}
    for h in header[1:]:
        w[h] = np.array([float(r[cols[h]]) for r in rows])
    return w


def classify_windows(w, c, sci):
    b = sci["regions"]["depth_classes"]
    ratio = w["depth_all"] / c
    cls = np.array([depth_class(x, b) for x in ratio], dtype=object)
    with np.errstate(divide="ignore", invalid="ignore"):
        mq0 = np.where(w["depth_all"] > 0, 1 - w["depth_hq"] / w["depth_all"], 0.0)
        mixed_kb = np.where(w["tested_bp"] > 0, w["mixed_sites"] / w["tested_bp"] * 1000, 0.0)
    return ratio, cls, mq0, mixed_kb


def class_bp(w, cls):
    length = w["end"] - w["start"]
    return {k: int(length[cls == k].sum()) for k in CLASSES}


def merge_regions(w, cls, keep, gap):
    """Merge consecutive windows (same contig) whose class is in `keep`,
    allowing up to `gap` other windows inside a region. Returns index lists."""
    regions, cur, last, gaps = [], [], None, 0
    for i in range(len(cls)):
        same_contig = cur and w["contig"][i] == w["contig"][cur[-1]]
        if cls[i] in keep:
            if cur and same_contig and gaps <= gap:
                cur.append(i)
            else:
                if cur:
                    regions.append(cur)
                cur = [i]
            gaps = 0
        elif cur:
            gaps += 1
            if gaps > gap or not same_contig:
                regions.append(cur)
                cur, gaps = [], 0
    if cur:
        regions.append(cur)
    return regions


def _region_record(w, idx, kind, ratio, mq0, mixed_kb, med_mixed, c, sci, clen):
    i0, i1 = idx[0], idx[-1]
    s, e = int(w["start"][i0]), int(w["end"][i1])
    sl = slice(i0, i1 + 1)
    length = w["end"][sl] - w["start"][sl]
    dr = float(np.average(ratio[sl], weights=length))
    mq = float(np.average(mq0[sl], weights=length))
    mk = float(np.average(mixed_kb[sl], weights=length))
    span_min = int(w["span_min"][sl].min())
    enrich = mk / med_mixed if med_mixed > 0 else (float("inf") if mk > 0 else 0.0)
    r = sci["regions"]
    near_end = s < sci["read_support"]["end_margin"] or e > clen - sci["read_support"]["end_margin"]
    evidence = [f"depth {dr:.2f}c", f"multi-mapping {mq:.2f}", f"mixed sites {mk:.2f}/kb ({enrich:.1f}x median)",
                f"min spanning reads {span_min}"]
    if kind == "collapse":
        if dr >= 1.75 and enrich >= r["het_enriched"] and (e - s) >= 20000:
            conf = "high"
        elif dr >= 1.75 or enrich >= r["het_enriched"]:
            conf = "medium"
        else:
            conf = "low"
        call = "collapsed_repeat_or_duplicate"
    elif kind == "half":
        if mq >= r["mq0_high"]:
            call, conf = "near_identical_duplicate", ("high" if (e - s) >= 20000 else "medium")
            evidence.append("reads split between near-identical copies")
        else:
            call, conf = "haploid_depth_region", ("medium" if (e - s) >= 20000 else "low")
            evidence.append("reads of one haplotype only: a second allele copy or a heterozygous SV")
    else:  # low support
        if near_end:
            call, conf = "low_support_contig_end", "low"
        elif span_min == 0 and dr < r["depth_classes"]["none"]:
            call, conf = "unsupported_sequence", "high"
        else:
            call, conf = "low_support", ("medium" if span_min <= 1 else "low")
        if mq >= r["mq0_high"]:
            evidence.append("mostly multi-mapping reads (repeat): support may be underestimated")
            conf = "low"
    return {"contig": w["contig"][i0], "start": s, "end": e, "length": e - s, "call": call,
            "confidence": conf, "depth_ratio": round(dr, 3), "multimap_frac": round(mq, 3),
            "mixed_per_kb": round(mk, 3), "span_min": span_min, "near_contig_end": int(near_end),
            "evidence": "; ".join(evidence)}


def classify_breakpoints(clusters, c, sci):
    """Grade clip clusters. clusters: list of dicts with CLUSTER_COLS keys."""
    bp = sci["regions"]["breakpoint"]
    out = []
    lo, hi = bp["het_sv_span_frac"]
    for k in clusters:
        n, span = int(k["n_clip"]), int(k["span"])
        clip_frac, span_frac = n / c, span / c
        ev = [f"{n} clipped ({k['n_split']} split, {k['n_left']} left / {k['n_right']} right)",
              f"{span} spanning unique reads ({span_frac:.2f}c)", f"depth {int(k['depth_all'])}"]
        if k.get("sa_partner", "NA") != "NA":
            ev.append(f"{k['sa_partner_reads']} split reads continue at {k['sa_partner']}")
        if int(k["near_end"]):
            call, conf = "contig_end_clipping", "info"
        elif n >= bp["min_clip_reads"] and clip_frac >= bp["high_clip_frac"] and span_frac <= bp["high_max_span_frac"]:
            call, conf = "suspected_misjoin", "high" if int(k["n_split"]) >= bp["min_clip_reads"] else "medium"
        elif n >= bp["min_clip_reads"] and lo * 0.5 <= clip_frac and lo <= span_frac <= hi:
            call, conf = "likely_heterozygous_sv", "medium"
        elif n >= bp["min_clip_reads"] and span_frac < 0.5:
            call, conf = "weakly_spanned_breakpoint", "medium"
        else:
            call, conf = "minor_clipping", "low"
        out.append({**k, "call": call, "confidence": conf, "clip_frac": round(clip_frac, 3),
                    "span_frac": round(span_frac, 3), "evidence": "; ".join(ev)})
    return out


def read_clusters(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def run(windows_path, clusters_path, summary_path, contig_lengths, sci, out_dir, kmer_peak=None,
        expected_setting=None):
    with open(summary_path) as fh:
        rs = json.load(fh)
    setting = expected_setting if expected_setting is not None else sci["regions"]["expected_depth"]
    c, source, warn = expected_depth(setting, rs["depth_mode_all"], kmer_peak)
    if c <= 0:
        raise ValueError("expected depth is 0; no reads aligned?")
    w = read_windows(windows_path)
    ratio, cls, mq0, mixed_kb = classify_windows(w, c, sci)
    normal = cls == "normal"
    med_mixed = float(np.median(mixed_kb[normal])) if normal.any() else 0.0
    gap = sci["regions"]["merge_gap_windows"]
    regions = []
    for kind, keep in (("collapse", {"dup", "multi"}), ("half", {"half"}), ("low", {"none", "low"})):
        for idx in merge_regions(w, cls, keep, gap):
            regions.append(_region_record(w, idx, kind, ratio, mq0, mixed_kb, med_mixed, c, sci,
                                          contig_lengths[w["contig"][idx[0]]]))
    regions.sort(key=lambda r: (r["contig"], r["start"]))
    bps = classify_breakpoints(read_clusters(clusters_path), c, sci)

    cols = ["contig", "start", "end", "length", "call", "confidence", "depth_ratio", "multimap_frac",
            "mixed_per_kb", "span_min", "near_contig_end", "evidence"]
    with open(f"{out_dir}/regions.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in regions:
            fh.write("\t".join(str(r[k]) for k in cols) + "\n")
    bcols = ["contig", "pos", "call", "confidence", "n_clip", "n_split", "span", "depth_all", "clip_frac",
             "span_frac", "near_end", "sa_partner", "sa_partner_reads", "evidence"]
    with open(f"{out_dir}/breakpoints.tsv", "w") as fh:
        fh.write("\t".join(bcols) + "\n")
        for k in bps:
            fh.write("\t".join(str(k[x]) for x in bcols) + "\n")
    with gzip.open(f"{out_dir}/window_classes.tsv.gz", "wt") as fh:
        fh.write("contig\tstart\tend\tdepth_ratio\tclass\tmultimap_frac\tmixed_per_kb\n")
        for i in range(len(cls)):
            fh.write(f"{w['contig'][i]}\t{int(w['start'][i])}\t{int(w['end'][i])}\t{ratio[i]:.3f}\t{cls[i]}\t"
                     f"{mq0[i]:.3f}\t{mixed_kb[i]:.3f}\n")

    def tally(items, key="call"):
        t = {}
        for x in items:
            k = f"{x[key]}|{x['confidence']}"
            t.setdefault(k, {"n": 0, "bp": 0})
            t[k]["n"] += 1
            t[k]["bp"] += int(x.get("length", 0))
        return t

    summary = {"expected_depth": c, "expected_depth_source": source, "warning": warn,
               "kmer_peak": kmer_peak, "depth_mode_all": rs["depth_mode_all"],
               "median_mixed_per_kb_normal": med_mixed,
               "class_bp": class_bp(w, cls),
               "multimapping_bp": int((w["end"] - w["start"])[mq0 >= sci["regions"]["mq0_high"]].sum()),
               "mixed_sites_total": rs["mixed_sites"], "error_sites_total": rs["error_sites"],
               "tested_bp": rs["tested_bp"],
               "error_sites_per_mb": rs["error_sites"] / max(1, rs["tested_bp"]) * 1e6,
               "regions": tally(regions), "breakpoints": tally(bps)}
    with open(f"{out_dir}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary
