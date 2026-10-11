"""Screen a meryl k-mer histogram for signs of residual heterozygosity.

This is a deliberately simple screening heuristic used to route samples to
human review; it does not estimate heterozygosity. GenomeScope-style model
fitting can replace it later.
"""


def read_histogram(path):
    hist = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                hist.append((int(parts[0]), int(parts[1])))
    return sorted(hist)


def _smooth(counts, xs, half):
    """Centered moving average over 2*half+1 bins (half=0: no smoothing)."""
    if half <= 0:
        return dict(counts)
    return {x: sum(counts.get(x + d, 0) for d in range(-half, half + 1)) / (2 * half + 1)
            for x in xs}


def screen(hist, half_window=(0.35, 0.65), max_ratio=0.25, min_local_max_ratio=0.05,
           smooth_frac=0.1, k=21):
    """Screen a k-mer spectrum for a heterozygous peak near half the main peak.

    Interpretation (idealised: Poisson coverage, diploid, k-mer length k,
    per-base heterozygosity h << 1/k): heterozygous k-mers number ~2*k*h*G and
    peak at half coverage with sqrt(2)-narrower spread, so the half-peak to
    main-peak height ratio is ~2*sqrt(2)*k*h. For k = 21 a ratio of 0.25
    corresponds to h ~ 0.4%. Real HiFi spectra are over-dispersed and maize
    repeats add mass at higher multiplicities, so this is only a screen.

    Local maxima are sought on a smoothed spectrum (window ~10% of the main
    peak coverage) and must reach min_local_max_ratio of the main peak, so
    sampling noise in sparse histograms is not reported as a peak."""
    if len(hist) < 10:
        return {"ok": False, "error": "histogram too short to analyse"}
    raw = dict(hist)
    xs = [x for x, _ in hist]
    start = xs[0]
    for x in xs[1:]:
        if raw.get(x, 0) > raw.get(x - 1, 0):
            break
        start = x
    peak_x = max((x for x in xs if x > start), key=lambda x: raw[x], default=None)
    if peak_x is None:
        return {"ok": False, "error": "no coverage peak found after the error region"}
    half = max(1, int(round(smooth_frac * peak_x / 2))) if smooth_frac else 0
    counts = _smooth(raw, xs, half)
    peak_h = counts[peak_x] or 1
    lo, hi = int(half_window[0] * peak_x), int(round(half_window[1] * peak_x))
    window = [(x, counts.get(x, 0)) for x in range(max(lo, start + 1), hi + 1)]
    win_max_x, win_max = max(window, key=lambda t: t[1]) if window else (None, 0)
    local_max = [x for x, c in window
                 if c > counts.get(x - 1, 0) and c > counts.get(x + 1, 0)
                 and x - 1 > start + half and c >= min_local_max_ratio * peak_h]
    ratio = win_max / peak_h
    flagged = bool(local_max) or ratio > max_ratio
    return {
        "ok": True,
        "error_region_end": start,
        "main_peak_coverage": peak_x,
        "main_peak_height": raw[peak_x],
        "smoothing_half_window": half,
        "half_peak_window": [lo, hi],
        "half_window_max_coverage": win_max_x,
        "half_window_max_height_smoothed": round(win_max, 1),
        "half_to_main_height_ratio": round(ratio, 4),
        "k": k,
        "idealised_heterozygosity_equiv": round(ratio / (2 * 2 ** 0.5 * k), 5),
        "half_window_local_maxima": local_max,
        "flagged": flagged,
        "flag_reason": ("local maximum near half the main peak" if local_max else
                        f"half-peak height ratio {ratio:.3f} > {max_ratio}" if ratio > max_ratio
                        else ""),
        "method": "heuristic screen v2 (smoothed; see scripts/g2f/kmer.py)",
    }
