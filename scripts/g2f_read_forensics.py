#!/usr/bin/env python3
"""Read-level genotype forensics: is a long-read library one diploid individual,
and do two libraries (e.g. HiFi and ONT) come from the same individual?

Works from reads aligned to an assembly of the same line (a collapsed /
purged primary, so that both alleles of a heterozygous region map to one
copy). Read-only; every output goes to a new directory.

Subcommands (each writes TSV, optionally gzipped):

  windows   sampled windows (BED) on long contigs: a systematic genome sample.
  spans     merged reference spans (BED) of primary alignments in a BAM/CRAM.
  snv       HiFi SNV sites over a BED (bcftools mpileup, parallel chunks):
            chrom, pos, ref, alt, ref_n, alt_n, dp. Only sites whose first
            alternative allele has >= --min-alt reads are kept.
  lookup    allele counts of one library at given positions (bcftools).
  af        allele-fraction summary of SNV sites at single-copy depth:
            histogram, binomial-mixture fits (1-4 components, BIC),
            per-window heterozygosity map.
  haps      local haplotype count: clusters of heterozygous sites spanned
            by single reads; distinct allele strings with >= --min-support
            reads. A single diploid individual has <= 2 per cluster except
            at collapsed paralogs (measured in inbred controls).
  ontobs    walk long-read alignments (MD tags required): substitutions
            with base quality >= --min-bq at positions that are not SNV
            sites of the reference library, and observations at the SNV
            sites; per-read and per-window denominators.
  ontsum    join ontobs with lookup counts; per read / per window rates of
            library-only alleles (absent from the other library at
            adequate depth) and allele use at heterozygous sites.

Interpretation needs controls: run the same steps on lines whose libraries
are known to agree. Nothing here is a contamination percentage.
"""
import argparse
import bisect
import collections
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

BASES = "ACGT"


def zopen(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def read_fai(path, min_len=0):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[1]) >= min_len:
                out.append((f[0], int(f[1])))
    return out


def read_bed(path):
    out = []
    with zopen(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track")):
                f = line.split("\t")
                out.append((f[0], int(f[1]), int(f[2])))
    return out


def merge_intervals(iv, gap=0):
    """Merge (chrom, start, end) intervals; input in any order."""
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2] + gap:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


# --------------------------------------------------------------------------
# windows / spans

def make_windows(contigs, win, step, offset=0):
    """Windows of length `win` every `step` bp on each (name, length) contig."""
    out = []
    for name, ln in contigs:
        s = offset
        while s + win <= ln:
            out.append((name, s, s + win))
            s += step
    return out


def cmd_windows(a):
    w = make_windows(read_fai(a.fai, a.min_contig), a.win, a.step, a.offset)
    with open(a.out, "w") as fh:
        for c, s, e in w:
            fh.write(f"{c}\t{s}\t{e}\n")
    log(f"{len(w)} windows, {sum(e - s for _, s, e in w)} bp")


def cmd_spans(a):
    import pysam
    iv = []
    with pysam.AlignmentFile(a.bam, reference_filename=a.ref) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
                continue
            iv.append((r.reference_name, max(0, r.reference_start - a.pad), r.reference_end + a.pad))
    m = merge_intervals(iv, a.gap)
    with open(a.out, "w") as fh:
        for c, s, e in m:
            fh.write(f"{c}\t{s}\t{e}\n")
    log(f"{len(iv)} alignments -> {len(m)} spans, {sum(e - s for _, s, e in m)} bp")


# --------------------------------------------------------------------------
# bcftools wrappers

def parse_ad(alt_field, ad_field):
    """'C,T,<*>' + '10,5,1,0' -> {'REF':10, 'C':5, 'T':1}; '<*>' dropped."""
    alts = alt_field.split(",")
    ad = [int(x) if x != "." else 0 for x in ad_field.split(",")]
    out = {"REF": ad[0] if ad else 0}
    for al, n in zip(alts, ad[1:]):
        if al != "<*>":
            out[al] = n
    return out


def _chunks(items, n, max_piece=20_000_000):
    """Split (chrom, start, end) regions into about n groups of similar total
    length; regions longer than the group size (or max_piece) are cut first,
    so one long contig cannot hold up a parallel run."""
    n = max(1, min(n, len(items) * 1000))
    total = sum(e - s for _, s, e in items)
    size = max(1, min(math.ceil(total / n), max_piece))
    pieces = []
    for c, s, e in items:
        while e - s > size:
            pieces.append((c, s, s + size))
            s += size
        pieces.append((c, s, e))
    out, cur, acc = [], [], 0
    for it in pieces:
        cur.append(it)
        acc += it[2] - it[1]
        if acc >= size:
            out.append(cur)
            cur, acc = [], 0
    if cur:
        out.append(cur)
    return out


def _run_pipeline(cmds, out_path):
    procs = []
    prev = None
    with open(out_path, "w") as out:
        for i, c in enumerate(cmds):
            last = i == len(cmds) - 1
            p = subprocess.Popen(c, stdin=prev, stdout=out if last else subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            if prev is not None:
                prev.close()
            prev = p.stdout
            procs.append(p)
        errs = [p.stderr.read().decode() for p in procs]
        rcs = [p.wait() for p in procs]
        for p in procs:
            p.stderr.close()
            if p.stdout:
                p.stdout.close()
    if any(rcs):
        raise RuntimeError(f"pipeline failed {rcs}: {cmds}\n{''.join(errs)[-2000:]}")


def mpileup_cmd(a, regions_file, targets_file=None):
    c = [a.bcftools, "mpileup", "-f", a.ref, "-R", regions_file, "-a", "FORMAT/AD,FORMAT/DP",
         "-q", str(a.min_mapq), "-Q", str(a.min_bq), "-I", "-d", str(a.max_depth), "-Ou"]
    if targets_file:
        c += ["-T", targets_file]
    return c + [a.aln]


def cmd_snv(a):
    regions = merge_intervals(read_bed(a.bed))
    chunks = _chunks(regions, a.threads * 4)
    tmp = tempfile.mkdtemp(prefix="snv.", dir=a.tmpdir)
    query = [a.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n"]
    filt = [a.bcftools, "view", "-i", f"FMT/AD[0:1]>={a.min_alt}", "-Ou"]

    def job(i):
        bed = os.path.join(tmp, f"c{i:04d}.bed")
        with open(bed, "w") as fh:
            for c, s, e in chunks[i]:
                fh.write(f"{c}\t{s}\t{e}\n")
        out = os.path.join(tmp, f"c{i:04d}.tsv")
        _run_pipeline([mpileup_cmd(a, bed), filt, query], out)
        return out

    with ThreadPoolExecutor(a.threads) as ex:
        outs = list(ex.map(job, range(len(chunks))))
    n = 0
    with zopen(a.out, "wt") as fo:
        fo.write("chrom\tpos\tref\talt\tref_n\talt_n\tdp\n")
        for o in outs:
            with open(o) as fh:
                for line in fh:
                    c, p, ref, alt, ad, dp = line.rstrip("\n").split("\t")
                    if len(ref) != 1:
                        continue
                    d = parse_ad(alt, ad)
                    al = alt.split(",")[0]
                    if al not in BASES:
                        continue
                    fo.write(f"{c}\t{p}\t{ref}\t{al}\t{d['REF']}\t{d[al]}\t{dp}\n")
                    n += 1
            os.remove(o)
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    log(f"{n} SNV sites over {sum(e - s for _, s, e in regions)} bp in {len(chunks)} chunks")


def cmd_lookup(a):
    """Counts per base at positions (chrom, pos 1-based) within regions."""
    regions = merge_intervals(read_bed(a.bed))
    pos = collections.defaultdict(list)
    with zopen(a.positions) as fh:
        for line in fh:
            c, p = line.split("\t")[:2]
            if p.strip().isdigit():
                pos[c].append(int(p))
    chunks = _chunks(regions, a.threads * 4)
    tmp = tempfile.mkdtemp(prefix="lookup.", dir=a.tmpdir)
    query = [a.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n"]

    def job(i):
        bed = os.path.join(tmp, f"c{i:04d}.bed")
        tgt = os.path.join(tmp, f"c{i:04d}.pos")
        k = 0
        with open(bed, "w") as fb, open(tgt, "w") as ft:
            for c, s, e in chunks[i]:
                fb.write(f"{c}\t{s}\t{e}\n")
                ps = sorted(set(pos.get(c, ())))
                for p in ps[bisect.bisect_right(ps, s):bisect.bisect_right(ps, e)]:
                    ft.write(f"{c}\t{p}\n")
                    k += 1
        out = os.path.join(tmp, f"c{i:04d}.tsv")
        if k:
            _run_pipeline([mpileup_cmd(a, bed, tgt), query], out)
        else:
            open(out, "w").close()
        return out

    with ThreadPoolExecutor(a.threads) as ex:
        outs = list(ex.map(job, range(len(chunks))))
    n = 0
    with zopen(a.out, "wt") as fo:
        fo.write("chrom\tpos\tref\tdp\tA\tC\tG\tT\n")
        for o in outs:
            with open(o) as fh:
                for line in fh:
                    c, p, ref, alt, ad, dp = line.rstrip("\n").split("\t")
                    if len(ref) != 1:
                        continue
                    d = parse_ad(alt, ad)
                    cnt = {b: d.get(b, 0) for b in BASES}
                    if ref in cnt:
                        cnt[ref] = d["REF"]
                    fo.write(f"{c}\t{p}\t{ref}\t{dp}\t" + "\t".join(str(cnt[b]) for b in BASES) + "\n")
                    n += 1
            os.remove(o)
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    log(f"{n} of {sum(len(set(v)) for v in pos.values())} distinct positions looked up")


# --------------------------------------------------------------------------
# allele fractions

def load_sites(path, dp_lo, dp_hi, min_alt):
    """SNV sites at single-copy depth: list of (chrom, pos, alt_n, dp_used)."""
    out = []
    with zopen(path) as fh:
        next(fh)
        for line in fh:
            c, p, ref, alt, rn, an, dp = line.rstrip("\n").split("\t")
            rn, an = int(rn), int(an)
            n = rn + an
            if an >= min_alt and dp_lo <= n <= dp_hi:
                out.append((c, int(p), an, n))
    return out


def binom_logpmf(k, n, p):
    p = min(max(p, 1e-9), 1 - 1e-9)
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log(1 - p))


def fit_binom_mixture(data, k, iters=200, init=None, tol=1e-7):
    """EM for a k-component binomial mixture on (alt, n) pairs.
    Returns (loglik, [(weight, p), ...]) with components sorted by p."""
    if not data:
        return float("-inf"), []
    ps = list(init) if init else [0.15 + 0.7 * (i + 0.5) / k for i in range(k)]
    ws = [1.0 / k] * k
    # group identical (alt, n) pairs: much faster on large sets
    grp = collections.Counter(data)
    keys = list(grp)
    cnt = [grp[x] for x in keys]
    tot = sum(cnt)
    ll_old = float("-inf")
    ll = ll_old
    for _ in range(iters):
        sw = [0.0] * k
        sk = [0.0] * k
        sn = [0.0] * k
        ll = 0.0
        for (a, n), c in zip(keys, cnt):
            lp = [math.log(max(ws[j], 1e-300)) + binom_logpmf(a, n, ps[j]) for j in range(k)]
            m = max(lp)
            s = sum(math.exp(x - m) for x in lp)
            ll += c * (m + math.log(s))
            for j in range(k):
                r = c * math.exp(lp[j] - m) / s
                sw[j] += r
                sk[j] += r * a
                sn[j] += r * n
        ws = [x / tot for x in sw]
        ps = [sk[j] / sn[j] if sn[j] > 0 else ps[j] for j in range(k)]
        if ll - ll_old < tol * abs(ll):
            break
        ll_old = ll
    comp = sorted(zip(ws, ps), key=lambda x: x[1])
    return ll, comp


def af_hist(values, edges):
    h = [0] * (len(edges) - 1)
    for v in values:
        i = bisect.bisect_right(edges, v) - 1
        if 0 <= i < len(h):
            h[i] += 1
        elif v == edges[-1]:
            h[-1] += 1
    return h


def cmd_af(a):
    c = a.expected_depth
    lo, hi = int(math.floor(a.dp_lo * c)), int(math.ceil(a.dp_hi * c))
    sites = load_sites(a.sites, lo, hi, a.min_alt)
    afs = [x[2] / x[3] for x in sites]
    edges = [i / 50 for i in range(51)]
    h = af_hist(afs, edges)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "af_hist.tsv"), "w") as fh:
        fh.write("af_lo\taf_hi\tsites\n")
        for i, n in enumerate(h):
            fh.write(f"{edges[i]:.2f}\t{edges[i + 1]:.2f}\t{n}\n")
    # mixture fits on "heterozygous-range" sites (AF 0.1-0.9): the error
    # tail (< 0.1) and hom-alt (> 0.9) are summarised separately.
    mid = [(x[2], x[3]) for x in sites if a.fit_lo <= x[2] / x[3] <= a.fit_hi]
    if len(mid) > a.max_fit_sites:          # deterministic thinning
        step = len(mid) / a.max_fit_sites
        mid = [mid[int(i * step)] for i in range(a.max_fit_sites)]
    fits = []
    for k in range(1, 5):
        best = None
        for init in ([None] if k == 1 else [None, [0.25, 0.5, 0.75, 0.33][:k] if k <= 3 else None]):
            ll, comp = fit_binom_mixture(mid, k, init=init)
            if best is None or ll > best[0]:
                best = (ll, comp)
        npar = 2 * k - 1
        bic = -2 * best[0] + npar * math.log(max(1, len(mid)))
        fits.append({"k": k, "loglik": best[0], "bic": bic,
                     "components": [{"weight": w, "p": p} for w, p in best[1]]})
    # binomial dispersion check: observed variance of AF vs binomial at p
    ll1 = fits[0]
    p1 = ll1["components"][0]["p"] if ll1["components"] else float("nan")
    obs_var = _var([x[0] / x[1] for x in mid]) if mid else float("nan")
    exp_var = (sum(p1 * (1 - p1) / x[1] for x in mid) / len(mid)) if mid else float("nan")
    # per-window map
    win = collections.defaultdict(list)
    for ch, p, an, n in sites:
        win[(ch, (p - 1) // a.window)].append(an / n)
    with open(os.path.join(a.outdir, "af_windows.tsv"), "w") as fh:
        fh.write("chrom\tstart\tend\tsites\thet_sites\tmedian_af_het\tlow_af_sites\n")
        for (ch, w) in sorted(win):
            v = win[(ch, w)]
            het = sorted(x for x in v if a.fit_lo <= x <= a.fit_hi)
            med = het[len(het) // 2] if het else float("nan")
            low = sum(1 for x in v if x < a.fit_lo)
            fh.write(f"{ch}\t{w * a.window}\t{(w + 1) * a.window}\t{len(v)}\t{len(het)}\t{med:.3f}\t{low}\n")
    summary = {
        "sites_file": os.path.abspath(a.sites), "expected_depth": c, "depth_range": [lo, hi],
        "min_alt": a.min_alt, "sites_at_single_copy_depth": len(sites),
        "het_range": [a.fit_lo, a.fit_hi], "het_range_sites_fitted": len(mid),
        "low_af_sites": sum(1 for x in afs if x < a.fit_lo),
        "high_af_sites": sum(1 for x in afs if x > a.fit_hi),
        "het_af_observed_variance": obs_var, "het_af_binomial_variance": exp_var,
        "dispersion_ratio": obs_var / exp_var if exp_var and exp_var == exp_var else float("nan"),
        "mixture_fits": fits,
        "best_k_by_bic": min(fits, key=lambda f: f["bic"])["k"],
    }
    with open(os.path.join(a.outdir, "af_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    log(json.dumps({k: v for k, v in summary.items() if k != "mixture_fits"}))


def _var(v):
    if len(v) < 2:
        return float("nan")
    m = sum(v) / len(v)
    return sum((x - m) ** 2 for x in v) / (len(v) - 1)


# --------------------------------------------------------------------------
# read allele extraction (CIGAR walk; no MD needed)

CIG_M, CIG_I, CIG_D, CIG_N, CIG_S, CIG_H, CIG_P, CIG_EQ, CIG_X = range(9)


def read_alleles(ref_start, cigartuples, seq, quals, sites):
    """Bases of one read at sorted 0-based reference positions `sites`.
    Returns {pos: (base, qual)}; positions in deletions map to ('-', 0)."""
    out = {}
    if not sites:
        return out
    rpos, qpos = ref_start, 0
    i = bisect.bisect_left(sites, rpos)
    n = len(sites)
    for op, ln in cigartuples:
        if i >= n:
            break
        if op in (CIG_M, CIG_EQ, CIG_X):
            end = rpos + ln
            while i < n and sites[i] < end:
                if sites[i] >= rpos:
                    q = qpos + (sites[i] - rpos)
                    out[sites[i]] = (seq[q], quals[q] if quals is not None else 30)
                i += 1
            rpos = end
            qpos += ln
        elif op in (CIG_D, CIG_N):
            end = rpos + ln
            while i < n and sites[i] < end:
                if sites[i] >= rpos:
                    out[sites[i]] = ("-", 0)
                i += 1
            rpos = end
        elif op in (CIG_I, CIG_S):
            qpos += ln
    return out


def cluster_sites(positions, min_sites, max_span):
    """Greedy non-overlapping clusters of >= min_sites sorted positions
    spanning <= max_span bp."""
    out = []
    i, n = 0, len(positions)
    while i < n:
        j = i
        while j + 1 < n and positions[j + 1] - positions[i] <= max_span:
            j += 1
        if j - i + 1 >= min_sites:
            out.append(positions[i:j + 1])
            i = j + 1
        else:
            i += 1
    return out


def count_haplotypes(strings, min_support, min_frac):
    """strings: allele strings (one per read, complete over the cluster).
    Returns (n_haplotypes, supports sorted desc) for strings supported by
    >= min_support reads and >= min_frac of reads."""
    c = collections.Counter(strings)
    tot = sum(c.values())
    sup = sorted((v for v in c.values() if v >= min_support and v >= min_frac * tot), reverse=True)
    return len(sup), sup


def _haps_contig(job):
    """Haplotype counts for the clusters of one contig (worker process)."""
    import pysam
    a, ch, pos, ref_alt = job
    rows = []
    with pysam.AlignmentFile(a.aln, reference_filename=a.ref) as bam:
        for cl in cluster_sites(pos, a.min_sites, a.max_span):
            strings = []
            tot = 0
            for r in bam.fetch(ch, cl[0], cl[-1] + 1):
                if (r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq
                        or r.reference_start > cl[0] or r.reference_end <= cl[-1]):
                    continue
                tot += 1
                al = read_alleles(r.reference_start, r.cigartuples, r.query_sequence,
                                  r.query_qualities, cl)
                s = []
                for p in cl:
                    b, q = al.get(p, ("-", 0))
                    ra = ref_alt[p]
                    if q < a.min_bq or b not in ra:
                        s = None
                        break
                    s.append("0" if b == ra[0] else "1")
                if s is not None:
                    strings.append("".join(s))
            nh, sup = count_haplotypes(strings, a.min_support, a.min_frac)
            rows.append((ch, cl[0], cl[-1] + 1, len(cl), len(strings), tot, nh, sup))
    return rows


def cmd_haps(a):
    import multiprocessing
    c = a.expected_depth
    lo, hi = int(a.dp_lo * c), int(math.ceil(a.dp_hi * c))
    het = collections.defaultdict(set)
    for ch, p, an, n in load_sites(a.sites, lo, hi, a.min_alt):
        if a.het_lo <= an / n <= a.het_hi:
            het[ch].add(p - 1)
    ref_alt = collections.defaultdict(dict)
    with zopen(a.sites) as fh:
        next(fh)
        for line in fh:
            ch, p, ref, alt = line.split("\t")[:4]
            if int(p) - 1 in het.get(ch, ()):
                ref_alt[ch][int(p) - 1] = (ref, alt)
    jobs = [(a, ch, sorted(het[ch]), ref_alt[ch]) for ch in sorted(het, key=lambda x: -len(het[x]))]
    res = collections.Counter()
    nclu = 0
    with multiprocessing.Pool(max(1, a.threads)) as pool:
        results = pool.map(_haps_contig, jobs, chunksize=1)
    rows = sorted(r for rr in results for r in rr)
    with zopen(a.out, "wt") as fo:
        fo.write("chrom\tstart\tend\tsites\treads_complete\treads_total\tn_haps\tsupports\n")
        for ch, s0, e0, ns, nc, tot, nh, sup in rows:
            res[min(nh, 4)] += 1
            nclu += 1
            fo.write(f"{ch}\t{s0}\t{e0}\t{ns}\t{nc}\t{tot}\t{nh}\t{','.join(map(str, sup))}\n")
    log(f"{nclu} clusters; haplotype counts {dict(sorted(res.items()))}")


# --------------------------------------------------------------------------
# long-read observations against another library's SNV sites

def cmd_ontobs(a):
    import pysam
    snv = collections.defaultdict(set)
    with zopen(a.sites) as fh:
        next(fh)
        for line in fh:
            ch, p = line.split("\t")[:2]
            snv[ch].add(int(p) - 1)
    keep = None
    if a.read_ids:
        keep = set(x.split("\t")[0].strip() for x in zopen(a.read_ids) if x.strip())
    W = a.window
    nr = 0
    with pysam.AlignmentFile(a.aln, reference_filename=a.ref) as bam, \
            zopen(os.path.join(a.outdir, "reads.tsv.gz"), "wt") as fr, \
            zopen(os.path.join(a.outdir, "obs.tsv.gz"), "wt") as fo, \
            zopen(os.path.join(a.outdir, "windows.tsv.gz"), "wt") as fw:
        fr.write("read\tchrom\tstart\tend\tmapq\tread_len\tquery_aligned\tn_supp\tmean_bq\t"
                 "aligned_q\tsubs_q\tsubs_q_at_snv\n")
        fo.write("read\tchrom\tpos\tref\tbase\tbq\tkind\n")
        fw.write("read\tchrom\twstart\taligned_q\tsubs_q\n")
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
                continue
            if keep is not None and r.query_name not in keep:
                continue
            if r.query_alignment_length < a.min_aligned:
                continue
            ch = r.reference_name
            s = snv.get(ch, set())
            q = r.query_qualities
            seq = r.query_sequence
            alq = subs = subs_snv = 0
            wa = collections.Counter()
            ws = collections.Counter()
            for qp, rp, rb in r.get_aligned_pairs(matches_only=True, with_seq=True):
                bq = q[qp] if q is not None else 30
                base = seq[qp]
                rb = rb.upper()
                if rp in s:
                    fo.write(f"{r.query_name}\t{ch}\t{rp + 1}\t{rb}\t{base}\t{bq}\tsnv\n")
                    if base != rb and bq >= a.min_bq:
                        subs_snv += 1
                    continue
                if bq < a.min_bq:
                    continue
                alq += 1
                wa[rp // W] += 1
                if base != rb and rb in BASES and base in BASES:
                    subs += 1
                    ws[rp // W] += 1
                    fo.write(f"{r.query_name}\t{ch}\t{rp + 1}\t{rb}\t{base}\t{bq}\tsub\n")
            supp = len(r.get_tag("SA").rstrip(";").split(";")) if r.has_tag("SA") else 0
            mbq = (sum(q) / len(q)) if q is not None and len(q) else float("nan")
            fr.write(f"{r.query_name}\t{ch}\t{r.reference_start}\t{r.reference_end}\t{r.mapping_quality}\t"
                     f"{r.infer_read_length()}\t{r.query_alignment_length}\t{supp}\t{mbq:.2f}\t"
                     f"{alq}\t{subs}\t{subs_snv}\n")
            for w in sorted(wa):
                fw.write(f"{r.query_name}\t{ch}\t{w * W}\t{wa[w]}\t{ws.get(w, 0)}\n")
            nr += 1
            if a.max_reads and nr >= a.max_reads:
                break
    log(f"{nr} reads written to {a.outdir}")


def classify_hifi_site(ref_n, alt_n, dp_lo, dp_hi, het_lo, het_hi):
    n = ref_n + alt_n
    if n < dp_lo or n > dp_hi:
        return "offdepth"
    f = alt_n / n
    if het_lo <= f <= het_hi:
        return "het"
    if f > het_hi:
        return "homalt"
    return "lowaf"


def cmd_ontsum(a):
    c = a.expected_depth
    dlo, dhi = int(a.dp_lo * c), int(math.ceil(a.dp_hi * c))
    sites = {}
    with zopen(a.sites) as fh:
        next(fh)
        for line in fh:
            ch, p, ref, alt, rn, an, dp = line.rstrip("\n").split("\t")
            sites[(ch, int(p))] = (ref, alt, classify_hifi_site(int(rn), int(an), dlo, dhi, a.het_lo, a.het_hi))
    look = {}
    with zopen(a.lookup) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            look[(f[0], int(f[1]))] = dict(zip(BASES, map(int, f[4:8])))
    per = collections.defaultdict(collections.Counter)
    perw = collections.defaultdict(collections.Counter)
    W = a.window
    with zopen(os.path.join(a.obsdir, "obs.tsv.gz")) as fh:
        next(fh)
        for line in fh:
            rd, ch, p, rb, b, bq, kind = line.rstrip("\n").split("\t")
            p = int(p)
            bq = int(bq)
            k = per[rd]
            if kind == "snv":
                ref, alt, cls = sites[(ch, p)]
                if cls != "het" or bq < a.min_bq:
                    continue
                k["het_obs"] += 1
                k["het_ref" if b == ref else "het_alt" if b == alt else "het_other"] += 1
                continue
            cnt = look.get((ch, p))
            dp = sum(cnt.values()) if cnt else 0
            if dp < dlo or dp > dhi:
                k["sub_offdepth"] += 1
                continue
            if cnt.get(b, 0) >= a.min_support:
                k["sub_known"] += 1
            elif cnt.get(b, 0) == 0:
                k["sub_only"] += 1
                perw[(rd, ch, (p - 1) // W * W)]["only"] += 1
            else:
                k["sub_weak"] += 1
    reads = {}
    with zopen(os.path.join(a.obsdir, "reads.tsv.gz")) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            reads[d["read"]] = d
    flag = set()
    if a.flagged:
        flag = set(x.split("\t")[0].strip() for x in zopen(a.flagged) if x.strip())
    cols = ["sub_only", "sub_known", "sub_weak", "sub_offdepth", "het_obs", "het_ref", "het_alt", "het_other"]
    with zopen(os.path.join(a.outdir, "reads_summary.tsv.gz"), "wt") as fo:
        fo.write("\t".join(["read", "chrom", "start", "end", "mapq", "read_len", "query_aligned", "n_supp",
                            "mean_bq", "aligned_q"] + cols + ["only_per_10kb", "flagged"]) + "\n")
        for rd, d in reads.items():
            k = per.get(rd, collections.Counter())
            alq = int(d["aligned_q"])
            rate = 1e4 * k["sub_only"] / alq if alq else float("nan")
            fo.write("\t".join([rd, d["chrom"], d["start"], d["end"], d["mapq"], d["read_len"],
                                d["query_aligned"], d["n_supp"], d["mean_bq"], d["aligned_q"]]
                               + [str(k[x]) for x in cols] + [f"{rate:.3f}", "1" if rd in flag else "0"]) + "\n")
    with zopen(os.path.join(a.obsdir, "windows.tsv.gz")) as fh, \
            zopen(os.path.join(a.outdir, "windows_summary.tsv.gz"), "wt") as fo:
        next(fh)
        fo.write("read\tchrom\twstart\taligned_q\tsubs_q\tsub_only\tflagged\n")
        for line in fh:
            rd, ch, w, alq, sq = line.rstrip("\n").split("\t")
            o = perw.get((rd, ch, int(w)), {}).get("only", 0)
            fo.write(f"{rd}\t{ch}\t{w}\t{alq}\t{sq}\t{o}\t{1 if rd in flag else 0}\n")
    log(f"{len(reads)} reads summarised")


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    def common_aln(p):
        p.add_argument("--aln", required=True, help="BAM/CRAM, coordinate-sorted and indexed")
        p.add_argument("--ref", required=True, help="reference FASTA (indexed)")
        p.add_argument("--min-mapq", type=int, default=20)
        p.add_argument("--min-bq", type=int, default=20)

    def bcf(p):
        p.add_argument("--bcftools", default="bcftools")
        p.add_argument("--threads", type=int, default=4)
        p.add_argument("--max-depth", type=int, default=500)
        p.add_argument("--tmpdir", default=None)

    def depth(p):
        p.add_argument("--expected-depth", type=float, required=True,
                       help="single-copy (haploid-sum) read depth of the library on this reference")
        p.add_argument("--dp-lo", type=float, default=0.7)
        p.add_argument("--dp-hi", type=float, default=1.4)

    p = sp.add_parser("windows")
    p.add_argument("--fai", required=True)
    p.add_argument("--min-contig", type=int, default=1_000_000)
    p.add_argument("--win", type=int, default=200_000)
    p.add_argument("--step", type=int, default=1_000_000)
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_windows)

    p = sp.add_parser("spans")
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", default=None)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--pad", type=int, default=0)
    p.add_argument("--gap", type=int, default=1000)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_spans)

    p = sp.add_parser("snv")
    common_aln(p)
    bcf(p)
    p.add_argument("--bed", required=True)
    p.add_argument("--min-alt", type=int, default=2)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_snv)

    p = sp.add_parser("lookup")
    common_aln(p)
    bcf(p)
    p.add_argument("--bed", required=True, help="regions containing the positions (index jumps)")
    p.add_argument("--positions", required=True, help="TSV chrom, pos (1-based)")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_lookup)

    p = sp.add_parser("af")
    depth(p)
    p.add_argument("--sites", required=True)
    p.add_argument("--min-alt", type=int, default=3)
    p.add_argument("--fit-lo", type=float, default=0.1)
    p.add_argument("--fit-hi", type=float, default=0.9)
    p.add_argument("--max-fit-sites", type=int, default=400_000)
    p.add_argument("--window", type=int, default=1_000_000)
    p.add_argument("--outdir", required=True)
    p.set_defaults(fn=cmd_af)

    p = sp.add_parser("haps")
    common_aln(p)
    depth(p)
    p.add_argument("--sites", required=True)
    p.add_argument("--min-alt", type=int, default=3)
    p.add_argument("--het-lo", type=float, default=0.2)
    p.add_argument("--het-hi", type=float, default=0.8)
    p.add_argument("--min-sites", type=int, default=4)
    p.add_argument("--max-span", type=int, default=8000)
    p.add_argument("--min-support", type=int, default=3)
    p.add_argument("--min-frac", type=float, default=0.1)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_haps)

    p = sp.add_parser("ontobs")
    common_aln(p)
    p.add_argument("--sites", required=True, help="SNV sites of the other library (snv output)")
    p.add_argument("--min-aligned", type=int, default=10_000)
    p.add_argument("--window", type=int, default=10_000)
    p.add_argument("--max-reads", type=int, default=0)
    p.add_argument("--read-ids", default=None)
    p.add_argument("--outdir", required=True)
    p.set_defaults(fn=cmd_ontobs)

    p = sp.add_parser("ontsum")
    depth(p)
    p.add_argument("--sites", required=True)
    p.add_argument("--lookup", required=True)
    p.add_argument("--obsdir", required=True)
    p.add_argument("--het-lo", type=float, default=0.2)
    p.add_argument("--het-hi", type=float, default=0.8)
    p.add_argument("--min-bq", type=int, default=20)
    p.add_argument("--min-support", type=int, default=2)
    p.add_argument("--window", type=int, default=10_000)
    p.add_argument("--flagged", default=None, help="read IDs to mark (e.g. k-mer flagged reads)")
    p.add_argument("--outdir", required=True)
    p.set_defaults(fn=cmd_ontsum)

    a = ap.parse_args(argv)
    for d in ("outdir",):
        if getattr(a, d, None):
            os.makedirs(getattr(a, d), exist_ok=True)
    a.fn(a)


if __name__ == "__main__":
    main()
