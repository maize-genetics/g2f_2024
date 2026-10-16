"""Publication figures (SVG + PDF) for assembly validation."""
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

# Categorical slots (fixed order): assembly A / + strand = blue, B / - strand = orange.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
TEXT, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# Diverging around single-copy depth: below (blue) - neutral gray - above (orange/red).
CLASS_COLORS = {"none": "#1c4f8f", "low": "#4f86c6", "half": "#9cc0e8", "normal": "#c9c8c3",
                "dup": "#f2a27c", "multi": "#c8452a", "na": "#ffffff"}
CLASS_ORDER = ["none", "low", "half", "normal", "dup", "multi"]

plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": TEXT,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.grid": False, "legend.frameon": False})


def save(fig, base):
    fig.tight_layout()
    for ext in ("svg", "pdf"):
        fig.savefig(f"{base}.{ext}")
    plt.close(fig)


def _chrom_key(name):
    m = re.search(r"(\d+)$", name)
    return (int(m.group(1)) if m else 10 ** 6, name)


def _offsets(names_lengths):
    off, acc = {}, 0
    for n, l in names_lengths:
        off[n] = acc
        acc += l
    return off, acc


def contig_order(assign, chroms):
    """Contigs grouped by assigned chromosome, ordered by median reference position."""
    rank = {c: i for i, c in enumerate(chroms)}
    keep = [a for a in assign.values() if a["chrom"] in rank]
    keep.sort(key=lambda a: (rank[a["chrom"]], a["ref_median_pos"]))
    return [(a["contig"], a["length"]) for a in keep]


def _segments(by_q, xoff, yoff, min_len):
    segs, cols = [], []
    for q, alns in by_q.items():
        if q not in yoff:
            continue
        for a in alns:
            if a["t"] not in xoff or a["qe"] - a["qs"] < min_len:
                continue
            x0, x1 = xoff[a["t"]] + a["ts"], xoff[a["t"]] + a["te"]
            y0, y1 = yoff[q] + a["qs"], yoff[q] + a["qe"]
            if a["strand"] == "-":
                y0, y1 = y1, y0
            segs.append([(x0 / 1e6, y0 / 1e6), (x1 / 1e6, y1 / 1e6)])
            cols.append(SERIES[0] if a["strand"] == "+" else SERIES[1])
    return segs, cols


def _dotplot_ax(ax, by_q, assign, ref_lengths, chroms, min_len, title, ref_label, asm_label):
    xoff, xtot = _offsets([(c, ref_lengths[c]) for c in chroms])
    order = contig_order(assign, chroms)
    yoff, ytot = _offsets(order)
    segs, cols = _segments(by_q, xoff, yoff, min_len)
    lc = LineCollection(segs, colors=cols, linewidths=0.8, rasterized=len(segs) > 20000)
    ax.add_collection(lc)
    for c in chroms[1:]:
        ax.axvline(xoff[c] / 1e6, color=GRID, lw=0.5, zorder=0)
    # chromosome boundaries on y: first contig assigned to each chromosome
    seen = set()
    for q, _ in order:
        c = assign[q]["chrom"]
        if c not in seen and seen:
            ax.axhline(yoff[q] / 1e6, color=GRID, lw=0.5, zorder=0)
        seen.add(c)
    ax.set_xlim(0, xtot / 1e6)
    ax.set_ylim(0, max(ytot, 1) / 1e6)
    ax.set_xticks([(xoff[c] + ref_lengths[c] / 2) / 1e6 for c in chroms])
    ax.set_xticklabels([re.sub(r"^chr", "", c) for c in chroms])
    ax.set_xlabel(f"{ref_label} chromosome")
    ax.set_ylabel(f"{asm_label} contigs (Mb, ordered by assignment)")
    ax.set_title(title, fontsize=9, color=TEXT, loc="left")
    return len(segs)


def dotplot_genome(by_q, assign, ref_lengths, chrom_re, min_len, title, ref_label, asm_label, base):
    chroms = sorted([c for c in ref_lengths if re.match(chrom_re, c)], key=_chrom_key)
    fig, ax = plt.subplots(figsize=(7.5, 7))
    n = _dotplot_ax(ax, by_q, assign, ref_lengths, chroms, min_len, title, ref_label, asm_label)
    ax.plot([], [], color=SERIES[0], label="forward")
    ax.plot([], [], color=SERIES[1], label="reverse")
    ax.legend(loc="upper left", fontsize=7)
    save(fig, base)
    return n


def dotplot_chrom_pair(chrom, ref_len, panels, min_len, ref_label, base):
    """Side-by-side dot plots of one reference chromosome, one panel per assembly.
    panels: list of (label, by_q, assign)."""
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4.4), squeeze=False)
    for ax, (label, by_q, assign) in zip(axes[0], panels):
        contigs = [(a["contig"], a["length"]) for a in sorted(
            (a for a in assign.values() if a["chrom"] == chrom), key=lambda a: a["ref_median_pos"])]
        yoff, ytot = _offsets(contigs)
        segs, cols = _segments(by_q, {chrom: 0}, yoff, min_len)
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.9, rasterized=len(segs) > 20000))
        for q, _ in contigs[1:]:
            ax.axhline(yoff[q] / 1e6, color=GRID, lw=0.4, zorder=0)
        ax.set_xlim(0, ref_len / 1e6)
        ax.set_ylim(0, max(ytot, 1) / 1e6)
        ax.set_xlabel(f"{ref_label} {chrom} (Mb)")
        ax.set_ylabel(f"{label}: {len(contigs)} contigs (Mb)")
        ax.set_title(label, fontsize=9, loc="left", color=TEXT)
    save(fig, base)


def depth_histograms(series, base, max_ratio=3.0):
    """series: list of (label, depths list, bases list, expected depth)."""
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    for i, (label, d, n, c) in enumerate(series):
        tot = sum(n) or 1
        x = [di / c for di in d if di / c <= max_ratio]
        y = [ni / tot * 100 for di, ni in zip(d, n) if di / c <= max_ratio]
        ax.plot(x, y, color=SERIES[i], lw=2, label=f"{label} (c = {c:.0f}x)")
    for v, lab in ((0.5, "half"), (1.0, "1c"), (2.0, "2c")):
        ax.axvline(v, color=GRID, lw=1, zorder=0)
        ax.text(v, ax.get_ylim()[1] * 0.98, lab, color=MUTED, fontsize=7, ha="center", va="top")
    ax.set_xlabel("HiFi read depth / expected single-copy depth")
    ax.set_ylabel("% of assembly bases")
    ax.legend(fontsize=7)
    save(fig, base)


def class_bars(rows, base, xlabel="Assembly sequence (Mb)"):
    """Horizontal stacked bars. rows: list of (label, {class: bp})."""
    fig, ax = plt.subplots(figsize=(6.5, 0.7 + 0.55 * len(rows)))
    for y, (label, bp) in enumerate(rows):
        left = 0
        for k in CLASS_ORDER:
            v = bp.get(k, 0) / 1e6
            if v <= 0:
                continue
            ax.barh(y, v, left=left, color=CLASS_COLORS[k], edgecolor="white", linewidth=1, height=0.6)
            if v > 60:
                ax.text(left + v / 2, y, f"{v:.0f}", ha="center", va="center", fontsize=6, color=TEXT)
            left += v
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[k]) for k in CLASS_ORDER]
    ax.legend(handles, ["<0.15c", "0.15-0.3c", "0.3-0.7c (half)", "0.7-1.4c", "1.4-2.5c (dup)", ">2.5c"],
              ncol=6, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.35 if len(rows) < 3 else -0.15))
    save(fig, base)


def pair_class_bar(summary, a_label, b_label, base):
    keys = [("shared", "shared with " + b_label, "#c9c8c3"),
            ("redundant_identical", "extra copy, >=99.9% identical", "#c8452a"),
            ("redundant_divergent", "extra copy, divergent (allelic)", "#f2a27c"),
            ("in_B_other_outputs", f"only in {b_label} haplotype outputs", "#9cc0e8"),
            ("A_unique", f"not in any {b_label} output", "#1c4f8f")]
    fig, ax = plt.subplots(figsize=(6.5, 1.9))
    left = 0
    for k, lab, col in keys:
        v = summary["class_bp"].get(k, 0) / 1e6
        ax.barh(0, v, left=left, color=col, edgecolor="white", linewidth=1, height=0.5, label=f"{lab} ({v:.0f} Mb)")
        left += v
    ax.set_yticks([0])
    ax.set_yticklabels([a_label])
    ax.set_xlabel(f"{a_label} sequence classified against {b_label} (Mb)")
    ax.legend(fontsize=6, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.45))
    save(fig, base)


def identity_hist(rows, base, label):
    """rows: list of (identity, bp)."""
    fig, ax = plt.subplots(figsize=(5, 3))
    x = [r[0] * 100 for r in rows if r[0] >= 0.95]
    y = [r[1] / 1e6 for r in rows if r[0] >= 0.95]
    ax.bar(x, y, width=0.09, color=SERIES[0])
    ax.set_xlabel("identity of extra copy to its counterpart (%)")
    ax.set_ylabel("Mb")
    ax.set_title(label, fontsize=9, loc="left")
    save(fig, base)


def coverage_tracks(series, chroms, base, ref_label):
    """Reference coverage per bin. series: list of (label, {chrom: [(start, frac)]})."""
    n = len(chroms)
    fig, axes = plt.subplots(n, 1, figsize=(7, 0.75 * n + 0.6), sharex=True, squeeze=False)
    for ax, c in zip(axes[:, 0], chroms):
        for i, (label, cov) in enumerate(series):
            pts = cov.get(c, [])
            ax.plot([p[0] / 1e6 for p in pts], [p[1] for p in pts], color=SERIES[i], lw=1,
                    label=label if c == chroms[0] else None)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(re.sub(r"^chr", "", c), rotation=0, ha="right", va="center")
        ax.set_yticks([])
    axes[0, 0].legend(fontsize=7, ncol=len(series), loc="lower right")
    axes[-1, 0].set_xlabel(f"{ref_label} position (Mb); y = fraction of 1 Mb bin covered by alignments")
    save(fig, base)
