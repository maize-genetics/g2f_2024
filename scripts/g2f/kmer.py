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


def screen(hist, half_window=(0.35, 0.65), max_ratio=0.25):
    if len(hist) < 10:
        return {"ok": False, "error": "histogram too short to analyse"}
    counts = dict(hist)
    xs = [x for x, _ in hist]
    # Skip the error peak: walk forward from coverage 1 until counts stop falling.
    start = xs[0]
    for x in xs[1:]:
        if counts.get(x, 0) > counts.get(x - 1, 0):
            break
        start = x
    peak_x = max((x for x in xs if x > start), key=lambda x: counts[x], default=None)
    if peak_x is None:
        return {"ok": False, "error": "no coverage peak found after the error region"}
    lo, hi = int(half_window[0] * peak_x), int(round(half_window[1] * peak_x))
    window = [(x, counts.get(x, 0)) for x in range(max(lo, start + 1), hi + 1)]
    win_max_x, win_max = max(window, key=lambda t: t[1]) if window else (None, 0)
    local_max = [x for x, c in window
                 if c > counts.get(x - 1, 0) and c > counts.get(x + 1, 0)
                 and x - 1 > start]
    ratio = win_max / counts[peak_x] if counts[peak_x] else 0
    flagged = bool(local_max) or ratio > max_ratio
    return {
        "ok": True,
        "error_region_end": start,
        "main_peak_coverage": peak_x,
        "main_peak_height": counts[peak_x],
        "half_peak_window": [lo, hi],
        "half_window_max_coverage": win_max_x,
        "half_window_max_height": win_max,
        "half_to_main_height_ratio": round(ratio, 4),
        "half_window_local_maxima": local_max,
        "flagged": flagged,
        "flag_reason": ("local maximum near half the main peak" if local_max else
                        f"half-peak height ratio {ratio:.3f} > {max_ratio}" if ratio > max_ratio
                        else ""),
        "method": "heuristic screen v1 (see scripts/g2f/kmer.py)",
    }
