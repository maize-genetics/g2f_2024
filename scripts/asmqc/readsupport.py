"""HiFi read support along an assembly, from a coordinate-sorted, indexed
alignment file (minimap2 map-hifi, MD tags, no secondary alignments).

One pass over the alignments, in parallel chunks, produces per window:
  depth_all    mean depth of primary + supplementary alignments
  depth_hq     the same for MAPQ >= min_mapq ("unique")
  span_min     minimum number of MAPQ >= min_mapq alignments that extend at
               least span_flank bp on both sides of a position (read continuity)
  mixed_sites  positions where a fraction het_frac of unique reads disagree
               with the assembly base (heterozygous-like, NucFreq-style)
  error_sites  positions where >= error_frac of unique reads disagree (likely
               consensus errors)
and clusters of clipped / split alignment ends (candidate breakpoints).

Mismatch positions come from the MD tag alone (reference coordinates), so no
per-base pileup is needed. Deletions count once, at their first base;
insertions are not counted.
"""
import gzip
import json
import multiprocessing as mp
import re

import numpy as np

MD_RE = re.compile(r"(\d+)|\^([A-Za-z]+)|([A-Za-z])")
WINDOW_COLS = ["contig", "start", "end", "depth_all", "depth_hq", "depth_all_median", "span_min",
               "span_mean", "zero_bp", "tested_bp", "mixed_sites", "error_sites"]
CLUSTER_COLS = ["contig", "pos", "n_clip", "n_clip_hq", "n_left", "n_right", "n_split", "depth_all",
                "depth_hq", "span", "near_end", "contig_len", "sa_partner", "sa_partner_reads"]
SA_BIN = 10000  # supplementary-alignment partner locations are pooled in bins of this size
HIST_MAX = 1000


def md_events(md, ref_start):
    """Reference positions of mismatches and deletion starts encoded in an MD tag."""
    pos, out = ref_start, []
    for num, dele, mis in MD_RE.findall(md):
        if num:
            pos += int(num)
        elif dele:
            out.append(pos)
            pos += len(dele)
        else:
            out.append(pos)
            pos += 1
    return out


def chunks_for(contigs, window, chunk):
    """[(contig, length, start, end)] with chunk boundaries on window multiples."""
    if chunk % window:
        raise ValueError("chunk must be a multiple of window")
    out = []
    for name, length in contigs:
        for s in range(0, length, chunk):
            out.append((name, length, s, min(s + chunk, length)))
    return out


def _window_stats(arr_all, arr_hq, span, mm, start, end, window, p):
    rows = []
    lo, hi = p["het_frac"]
    tested = arr_hq >= p["het_min_depth"]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(tested, mm / np.maximum(arr_hq, 1), 0.0)
    alt_ok = mm >= p["het_min_alt"]
    mixed = tested & alt_ok & (frac >= lo) & (frac < hi)
    err = tested & alt_ok & (frac >= p["error_frac"])
    for ws in range(start, end, window):
        a, b = ws - start, min(ws + window, end) - start
        da = arr_all[a:b]
        rows.append([ws, ws + (b - a), round(float(da.mean()), 2), round(float(arr_hq[a:b].mean()), 2),
                     float(np.median(da)), int(span[a:b].min()), round(float(span[a:b].mean()), 2),
                     int((da == 0).sum()), int(tested[a:b].sum()), int(mixed[a:b].sum()), int(err[a:b].sum())])
    return rows


def process_chunk(args):
    path, ref, contig, clen, start, end, p = args
    import pysam
    L = end - start
    d_all = np.zeros(L + 1, np.int32)
    d_hq = np.zeros(L + 1, np.int32)
    d_span = np.zeros(L + 1, np.int32)
    mm = np.zeros(L, np.int32)
    clips = []  # (pos, side 0=left/1=right, split, hq, SA partner "contig:bin" or "")
    f = p["span_flank"]
    with pysam.AlignmentFile(path, reference_filename=ref) as af:
        for r in af.fetch(contig, start, end):
            if r.is_unmapped or r.is_secondary or r.is_qcfail or r.is_duplicate:
                continue
            rs, re_ = r.reference_start, r.reference_end
            a, b = max(rs, start), min(re_, end)
            if a < b:
                d_all[a - start] += 1
                d_all[b - start] -= 1
            hq = r.mapping_quality >= p["min_mapq"]
            if hq:
                if a < b:
                    d_hq[a - start] += 1
                    d_hq[b - start] -= 1
                sa, sb = max(rs + f, start), min(re_ - f, end)
                if sa < sb:
                    d_span[sa - start] += 1
                    d_span[sb - start] -= 1
                if r.has_tag("MD"):
                    ev = [x - start for x in md_events(r.get_tag("MD"), rs) if start <= x < end]
                    if ev:
                        np.add.at(mm, ev, 1)
            ct = r.cigartuples
            split = r.has_tag("SA")
            partner = ""
            if split:  # first other alignment of the read: contig and position bin
                sa = r.get_tag("SA").split(";")[0].split(",")
                partner = f"{sa[0]}:{int(sa[1]) // SA_BIN * SA_BIN}"
            if ct[0][0] in (4, 5) and ct[0][1] >= p["clip_min"] and start <= rs < end:
                clips.append((rs, 0, split, hq, partner))
            if ct[-1][0] in (4, 5) and ct[-1][1] >= p["clip_min"] and start <= re_ < end:
                clips.append((re_, 1, split, hq, partner))
    arr_all = np.cumsum(d_all)[:L]
    arr_hq = np.cumsum(d_hq)[:L]
    span = np.cumsum(d_span)[:L]
    rows = [[contig] + r for r in _window_stats(arr_all, arr_hq, span, mm, start, end, p["window"], p)]
    clusters = []
    if clips:
        clips.sort()
        groups, cur = [], [clips[0]]
        for c in clips[1:]:
            if c[0] - cur[-1][0] <= p["cluster_dist"]:
                cur.append(c)
            else:
                groups.append(cur)
                cur = [c]
        groups.append(cur)
        for g in groups:
            if len(g) < p["min_cluster_reads"]:
                continue
            pos = int(np.median([c[0] for c in g]))
            i = min(max(pos - start, 0), L - 1)
            near = pos < p["end_margin"] or pos > clen - p["end_margin"]
            partners = {}
            for c in g:
                if c[4]:
                    partners[c[4]] = partners.get(c[4], 0) + 1
            top = max(partners, key=partners.get) if partners else ""
            clusters.append([contig, pos, len(g), sum(c[3] for c in g), sum(1 for c in g if c[1] == 0),
                             sum(1 for c in g if c[1] == 1), sum(1 for c in g if c[2]), int(arr_all[i]),
                             int(arr_hq[i]), int(span[i]), int(near), clen, top or "NA", partners.get(top, 0)])
    h_all = np.bincount(np.minimum(arr_all, HIST_MAX), minlength=HIST_MAX + 1)
    h_hq = np.bincount(np.minimum(arr_hq, HIST_MAX), minlength=HIST_MAX + 1)
    return rows, clusters, h_all, h_hq


def run(aln_path, ref_fasta, contigs, params, workers, out_prefix):
    """Write <prefix>windows.tsv.gz, clip_clusters.tsv.gz, depth_hist.tsv, summary.json."""
    jobs = [(aln_path, ref_fasta, c, l, s, e, params)
            for c, l, s, e in chunks_for(contigs, params["window"], params["chunk"])]
    h_all = np.zeros(HIST_MAX + 1, np.int64)
    h_hq = np.zeros(HIST_MAX + 1, np.int64)
    n_win = n_clu = 0
    totals = {"mixed_sites": 0, "error_sites": 0, "tested_bp": 0, "zero_bp": 0}
    with gzip.open(out_prefix + "windows.tsv.gz", "wt") as fw, \
            gzip.open(out_prefix + "clip_clusters.tsv.gz", "wt") as fc:
        fw.write("\t".join(WINDOW_COLS) + "\n")
        fc.write("\t".join(CLUSTER_COLS) + "\n")
        ctx = mp.get_context("fork")
        with ctx.Pool(max(1, workers)) as pool:
            for rows, clusters, ha, hh in pool.imap(process_chunk, jobs, chunksize=1):
                for r in rows:
                    fw.write("\t".join(map(str, r)) + "\n")
                    totals["zero_bp"] += r[8]
                    totals["tested_bp"] += r[9]
                    totals["mixed_sites"] += r[10]
                    totals["error_sites"] += r[11]
                n_win += len(rows)
                for c in clusters:
                    fc.write("\t".join(map(str, c)) + "\n")
                n_clu += len(clusters)
                h_all += ha
                h_hq += hh
    with open(out_prefix + "depth_hist.tsv", "w") as fh:
        fh.write("depth\tbases_all\tbases_hq\n")
        for d in range(HIST_MAX + 1):
            fh.write(f"{d}\t{h_all[d]}\t{h_hq[d]}\n")
    summary = {"windows": n_win, "clip_clusters": n_clu, "assembly_bp": int(h_all.sum()),
               "depth_mode_all": depth_mode(h_all), "depth_mode_hq": depth_mode(h_hq),
               "mean_depth_all": float((np.arange(HIST_MAX + 1) * h_all).sum() / max(1, h_all.sum())),
               **totals}
    with open(out_prefix + "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def depth_mode(hist, min_depth=5):
    """Mode of a per-base depth histogram (3-bin smoothing, depth >= min_depth)."""
    h = np.asarray(hist, float)
    if h[min_depth:].sum() == 0:
        return 0
    sm = np.convolve(h, np.ones(3) / 3, mode="same")
    sm[:min_depth] = 0
    sm[-2:] = 0  # the HIST_MAX overflow bin and its smoothed neighbour
    return int(np.argmax(sm))
