"""Whole-genome alignment of an assembly to a reference (minimap2 PAF):
contig assignment, synteny segments, structural junctions and reference
coverage.

A junction is a boundary inside a contig between two synteny segments that
disagree in chromosome, orientation or order. It is a candidate only: maize
inbreds carry real inversions, translocated segments and presence/absence
variation relative to B73. Junctions are graded later against HiFi read
continuity (evidence stage), and are compared between assemblies.
"""
import json
import re

from .paf import merge_intervals, read_paf


def load_alignments(paf_path, sci):
    p = sci["ref_align"]
    by_q = {}
    for a in read_paf(paf_path, min_len=p["plot_min_len"], min_mapq=0):
        by_q.setdefault(a["q"], []).append(a)
    return by_q


def assign_contigs(by_q, contig_lengths, sci, chrom_re):
    """Per contig: aligned bp per chromosome (structure-grade alignments)."""
    p = sci["ref_align"]
    rx = re.compile(chrom_re)
    out = {}
    for q, qlen in contig_lengths.items():
        per = {}
        for a in by_q.get(q, []):
            if a["mapq"] < p["min_mapq"] or a["qe"] - a["qs"] < p["min_block_len"]:
                continue
            t = a["t"] if rx.match(a["t"]) else "other"
            per.setdefault(t, []).append((a["qs"], a["qe"]))
        bp = {t: sum(e - s for s, e in merge_intervals(iv)) for t, iv in per.items()}
        total = sum(bp.values())
        ranked = sorted(bp.items(), key=lambda x: -x[1])
        best = ranked[0][0] if ranked else ""
        second = [t for t, b in ranked[1:] if t != "other" and b >= p["multi_chr_min_bp"]
                  and b >= p["multi_chr_min_frac"] * total]
        ref_pos = 0
        if best and best != "other":
            pos = [(a["ts"] + a["te"]) / 2 for a in by_q.get(q, []) if a["t"] == best]
            ref_pos = sorted(pos)[len(pos) // 2] if pos else 0
        out[q] = {"contig": q, "length": qlen, "aligned_bp": total, "aligned_frac": total / qlen if qlen else 0,
                  "chrom": best, "chrom_bp": ranked[0][1] if ranked else 0, "other_chroms": ",".join(second),
                  "ref_median_pos": int(ref_pos)}
    return out


def synteny_segments(alns, sci):
    """Greedy collinear merge of one contig's alignments (sorted by query)."""
    p = sci["ref_align"]
    alns = sorted((a for a in alns if a["mapq"] >= p["min_mapq"] and a["qe"] - a["qs"] >= p["min_block_len"]),
                  key=lambda a: a["qs"])

    def compatible(seg, a):
        if a["t"] != seg["t"] or a["strand"] != seg["strand"]:
            return False
        if seg["strand"] == "+":
            gap = a["ts"] - seg["te"]
        else:
            gap = seg["ts"] - a["te"]
        return -p["synteny_merge_gap"] <= gap <= p["synteny_merge_gap"]

    segs = []
    for a in alns:
        if segs and compatible(segs[-1], a):
            s = segs[-1]
            s["qe"] = max(s["qe"], a["qe"])
            s["ts"], s["te"] = min(s["ts"], a["ts"]), max(s["te"], a["te"])
            s["n"] += 1
            s["aligned"] += a["qe"] - a["qs"]
        else:
            segs.append({"t": a["t"], "strand": a["strand"], "qs": a["qs"], "qe": a["qe"], "ts": a["ts"],
                         "te": a["te"], "n": 1, "aligned": a["qe"] - a["qs"]})
    big = [s for s in segs if s["qe"] - s["qs"] >= p["min_segment_len"]]
    merged = []  # re-merge large segments separated only by small ones
    for s in big:
        if merged and compatible(merged[-1], s):
            m = merged[-1]
            m["qe"], m["ts"], m["te"] = max(m["qe"], s["qe"]), min(m["ts"], s["ts"]), max(m["te"], s["te"])
            m["n"] += s["n"]
            m["aligned"] += s["aligned"]
        else:
            merged.append(dict(s))
    return merged


def junctions(contig, segs, sci):
    p = sci["ref_align"]
    out = []
    for a, b in zip(segs, segs[1:]):
        if a["t"] != b["t"]:
            kind = "inter_chromosomal"
        elif a["strand"] != b["strand"]:
            kind = "inversion"
        else:
            jump = (b["ts"] - a["te"]) if a["strand"] == "+" else (a["ts"] - b["te"])
            if abs(jump) <= p["ref_jump"]:
                continue
            kind = "intra_chromosomal_jump"
        s, e = sorted((a["qe"], b["qs"]))
        out.append({"contig": contig, "q_start": s, "q_end": e, "kind": kind,
                    "left": f"{a['t']}:{a['ts']}-{a['te']}({a['strand']})",
                    "right": f"{b['t']}:{b['ts']}-{b['te']}({b['strand']})",
                    "left_len": a["qe"] - a["qs"], "right_len": b["qe"] - b["qs"]})
    return out


def reference_coverage(by_q, ref_lengths, sci, chrom_re):
    """Covered fraction of each reference bin by structure-grade alignments."""
    p = sci["ref_align"]
    rx = re.compile(chrom_re)
    per_t = {}
    for alns in by_q.values():
        for a in alns:
            if a["mapq"] >= p["min_mapq"] and rx.match(a["t"]):
                per_t.setdefault(a["t"], []).append((a["ts"], a["te"]))
    binsz = p["coverage_bin"]
    rows = []
    for t, tlen in ref_lengths.items():
        if not rx.match(t):
            continue
        iv = merge_intervals(per_t.get(t, []))
        nb = (tlen + binsz - 1) // binsz
        cov = [0] * nb
        for s, e in iv:
            for b in range(s // binsz, (e - 1) // binsz + 1):
                lo, hi = max(s, b * binsz), min(e, (b + 1) * binsz)
                cov[b] += hi - lo
        for b in range(nb):
            blen = min(binsz, tlen - b * binsz)
            rows.append((t, b * binsz, b * binsz + blen, cov[b] / blen))
    return rows


def absent_regions(cov_rows, n_gap_bins, min_bins, max_cov=0.05):
    """Runs of >= min_bins consecutive bins with coverage <= max_cov (skipping
    bins that are mostly assembly gaps in the reference)."""
    out, run = [], []
    for r in cov_rows + [("", 0, 0, 1.0)]:
        gap = n_gap_bins.get((r[0], r[1]), 0.0) > 0.5
        if r[3] <= max_cov and not gap and (not run or run[-1][0] == r[0]):
            run.append(r)
            continue
        if len(run) >= min_bins:
            out.append({"chrom": run[0][0], "start": run[0][1], "end": run[-1][2],
                        "length": run[-1][2] - run[0][1]})
        run = [r] if (r[3] <= max_cov and not gap and r[0]) else []
    return out


def run(paf_path, contig_lengths, ref_lengths, ref_ngaps, sci, chrom_re, out_dir):
    by_q = load_alignments(paf_path, sci)
    assign = assign_contigs(by_q, contig_lengths, sci, chrom_re)
    all_j, seg_rows = [], []
    for q in contig_lengths:
        segs = synteny_segments(by_q.get(q, []), sci)
        for s in segs:
            seg_rows.append([q, s["qs"], s["qe"], s["t"], s["ts"], s["te"], s["strand"], s["n"], s["aligned"]])
        all_j.extend(junctions(q, segs, sci))
        assign[q]["n_segments"] = len(segs)
    cov = reference_coverage(by_q, ref_lengths, sci, chrom_re)
    absent = absent_regions(cov, ref_ngaps, sci["ref_align"]["min_absent_bins"])

    acols = ["contig", "length", "aligned_bp", "aligned_frac", "chrom", "chrom_bp", "other_chroms",
             "ref_median_pos", "n_segments"]
    with open(f"{out_dir}/contig_assignment.tsv", "w") as fh:
        fh.write("\t".join(acols) + "\n")
        for q in sorted(assign, key=lambda x: -assign[x]["length"]):
            a = assign[q]
            fh.write("\t".join(str(round(a[k], 4) if isinstance(a[k], float) else a[k]) for k in acols) + "\n")
    with open(f"{out_dir}/synteny_segments.tsv", "w") as fh:
        fh.write("contig\tq_start\tq_end\tref\tref_start\tref_end\tstrand\tn_alignments\taligned_bp\n")
        for r in seg_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    jcols = ["contig", "q_start", "q_end", "kind", "left", "right", "left_len", "right_len"]
    with open(f"{out_dir}/junctions.tsv", "w") as fh:
        fh.write("\t".join(jcols) + "\n")
        for j in all_j:
            fh.write("\t".join(str(j[k]) for k in jcols) + "\n")
    with open(f"{out_dir}/reference_coverage.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tcovered_frac\tref_gap_frac\n")
        for r in cov:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{ref_ngaps.get((r[0], r[1]), 0.0):.4f}\n")
    with open(f"{out_dir}/reference_absent.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tlength\n")
        for r in absent:
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['length']}\n")
    total = sum(contig_lengths.values())
    chrom_ref = sum(l for t, l in ref_lengths.items() if re.match(chrom_re, t))
    covered = sum((r[2] - r[1]) * r[3] for r in cov)
    multi = [a for a in assign.values() if a["other_chroms"]]
    summary = {
        "assembly_bp": total,
        "assigned_bp": sum(a["length"] for a in assign.values() if a["chrom"] not in ("", "other")),
        "aligned_bp": sum(a["aligned_bp"] for a in assign.values()),
        "reference_chrom_bp": chrom_ref, "reference_covered_bp": int(covered),
        "reference_covered_frac": covered / chrom_ref if chrom_ref else 0,
        "multi_chromosome_contigs": len(multi),
        "multi_chromosome_contig_bp": sum(a["length"] for a in multi),
        "junctions": {k: sum(1 for j in all_j if j["kind"] == k)
                      for k in ("inter_chromosomal", "inversion", "intra_chromosomal_jump")},
        "absent_regions": len(absent), "absent_bp": sum(r["length"] for r in absent),
    }
    with open(f"{out_dir}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary, by_q, assign
