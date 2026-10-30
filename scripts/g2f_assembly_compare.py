#!/usr/bin/env python3
"""Contiguity and completeness table for a set of assemblies, from existing QC outputs (nothing is recomputed except
length statistics from the FASTA index).

Assemblies table (TSV, one row per FASTA; a haplotype pair is two rows):
  name, line, assembly_class, role, fai                                         required
  assembly_class   free text, used to group fair comparisons, e.g. inbred_primary, phased_haplotype,
                   unpurged_primary, purged_primary
  role             e.g. release_target, companion, alternative
  busco            BUSCO short_summary*.txt                                     optional
  merqury_qv, merqury_qv_row         Merqury *.qv file and the first-column value of the row to use
  merqury_cn, merqury_cn_row         Merqury *.completeness.stats file and row
  merqury_cn_pair_row                row of the same file for a haplotype pair together (e.g. 'both')
  riparian_dir, riparian_unit        a g2f_riparian run and unit for gene-anchor statistics
  notes                              copied through (e.g. documented redundancy, with its source)

Per assembly: total bp, contigs, N50 / L50, NG50 / LG50 (--genome-size), largest contig, contigs and bp fraction
>= 1 / 10 / 50 Mb, BUSCO C/S/D/F/M, Merqury QV and k-mer completeness, and (with a riparian run) anchor-bearing
records and their bp, bp inside drawn anchor blocks, unanchored bp, records per reference chromosome (median,
max) and the reference span covered by drawn anchor blocks.

Outputs (new directory): assembly_compare.tsv, METHODS.md, RUN.json, figures/ngx_by_line.{pdf,svg,png}.

Usage:
  g2f_assembly_compare.py run --assemblies A.tsv --genome-size 2.3e9 --out OUTDIR
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

CODE_VERSION = 1
REQ = ("name", "line", "assembly_class", "role", "fai")
CLASS_COL = {"inbred_primary": "#2a78d6", "phased_haplotype": "#1baf7a", "unpurged_primary": "#eb6834",
             "purged_primary": "#6250d6"}
OTHER_COL = "#8a8984"


def die(msg):
    sys.exit(f"g2f_assembly_compare: error: {msg}")


def read_tsv(path):
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def lengths(fai):
    with open(fai) as fh:
        return sorted((int(ln.split("\t")[1]) for ln in fh if ln.strip()), reverse=True)


def nx(lens, total, frac=0.5):
    acc = 0
    for i, x in enumerate(lens, 1):
        acc += x
        if acc >= total * frac:
            return x, i
    return 0, ""


def contiguity(lens, gsize):
    tot = sum(lens)
    n50, l50 = nx(lens, tot)
    ng50, lg50 = nx(lens, gsize) if gsize else ("", "")
    out = dict(total_bp=tot, contigs=len(lens), N50=n50, L50=l50, NG50=ng50, LG50=lg50, largest=lens[0] if lens else 0,
               size_vs_genome=round(tot / gsize, 3) if gsize else "")
    for t, lab in ((1e6, "1Mb"), (10e6, "10Mb"), (50e6, "50Mb")):
        sel = [x for x in lens if x >= t]
        out[f"n_ge_{lab}"] = len(sel)
        out[f"frac_bp_ge_{lab}"] = round(sum(sel) / tot, 4) if tot else 0
    return out


def busco(path):
    txt = open(path).read()
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", txt)
    if not m:
        die(f"{path}: no BUSCO summary line")
    lin = re.search(r"lineage dataset is: (\S+)", txt)
    return dict(busco_lineage=lin.group(1) if lin else "", busco_C=float(m.group(1)), busco_S=float(m.group(2)),
                busco_D=float(m.group(3)), busco_F=float(m.group(4)), busco_M=float(m.group(5)),
                busco_n=int(m.group(6)))


def merqury_row(path, key, col):
    rows = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    hit = [r for r in rows if r[0] == key] if key else (rows if len(rows) == 1 else [])
    if len(hit) != 1:
        die(f"{path}: need exactly one row {key!r} (rows: {[r[0] for r in rows]})")
    return float(hit[0][col])


def anchored(ripdir, unit):
    rows = [r for r in read_tsv(os.path.join(ripdir, "sequences.tsv")) if r["unit"] == unit]
    if not rows:
        die(f"{ripdir}: unit {unit} not found")
    ribs = [r for r in read_tsv(os.path.join(ripdir, "ribbons.tsv")) if r["unit"] == unit and r["drawn"] == "yes"]
    pl = [r for r in rows if r["plotted"] == "yes"]
    per = defaultdict(int)
    for r in pl:
        per[r["dominant_chr"]] += 1
    qcov = defaultdict(list)
    rcov = defaultdict(list)
    for r in ribs:
        qcov[r["record"]].append((int(r["query_start"]), int(r["query_end"])))
        rcov[r["ref_chr"]].append((int(r["ref_start"]), int(r["ref_end"])))

    def union(iv):
        tot, cs, ce = 0, None, None
        for s, e in sorted(iv):
            if ce is None or s > ce + 1:
                if ce is not None:
                    tot += ce - cs + 1
                cs, ce = s, e
            else:
                ce = max(ce, e)
        return tot + (ce - cs + 1 if ce is not None else 0)
    counts = sorted(per.values())
    return dict(records_anchor_bearing=len(pl), bp_anchor_bearing=sum(int(r["length_bp"]) for r in pl),
                bp_in_anchor_blocks=sum(union(v) for v in qcov.values()),
                bp_unanchored=sum(int(r["length_bp"]) for r in rows if r["plotted"] != "yes"),
                records_per_chr_median=counts[len(counts) // 2] if counts else "",
                records_per_chr_max=counts[-1] if counts else "",
                ref_bp_covered_by_blocks=sum(union(v) for v in rcov.values()),
                gene_anchors=sum(int(r["gene_anchors"]) for r in rows))


def plot(out, rows, lens_by, gsize):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    lines = list(dict.fromkeys(r["line"] for r in rows))
    ncol = 3
    nrow = (len(lines) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.3 * nrow + 0.9), sharex=True, sharey=True, squeeze=False)
    for ax, line in zip(axes.flat, lines):
        for r in [x for x in rows if x["line"] == line]:
            lens = lens_by[r["name"]]
            xs, ys, acc = [0.0], [lens[0] / 1e6], 0
            for x in lens:
                acc += x
                xs.append(min(acc / gsize, 1.2) * 100)
                ys.append(x / 1e6)
                if acc / gsize > 1.2:
                    break
            ls = "--" if r["name"].endswith("hap2") else ":" if r["role"] != "release_target" else "-"
            ax.step(xs, ys, where="pre", color=CLASS_COL.get(r["assembly_class"], OTHER_COL), lw=2 if
                    r["role"] == "release_target" else 1.4, ls=ls, label=r["name"])
        ax.axvline(50, color="#d9d8d3", lw=1, zorder=0)
        ax.set_title(line, loc="left", fontsize=10)
        ax.set_yscale("log")
        ax.grid(axis="y", color="#eeeeea", lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(fontsize=6.5, frameon=False, loc="lower left")
    for ax in axes.flat[len(lines):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel(f"cumulative % of genome size ({gsize / 1e9:.2f} Gb)")
    for ax in axes[:, 0]:
        ax.set_ylabel("contig length (Mb, log)")
    h = [Line2D([], [], color=c, lw=2, label=k.replace("_", " ")) for k, c in CLASS_COL.items()]
    h += [Line2D([], [], color=OTHER_COL, lw=1.4, ls=":", label="not a release target (dotted)"),
          Line2D([], [], color=OTHER_COL, lw=1.4, ls="--", label="hap2 (dashed)")]
    fig.legend(handles=h, loc="upper center", ncol=6, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("NGx curves (vertical line = NG50); each haplotype is its own curve", x=0.01, ha="left", y=0.985,
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.join(out, "figures"))
    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(out, "figures", f"ngx_by_line.{ext}")
        fig.savefig(p, dpi=150 if ext == "png" else None)
        paths.append(p)
    plt.close(fig)
    return paths


COLS = ["line", "name", "assembly_class", "role", "total_bp", "size_vs_genome", "contigs", "N50", "L50", "NG50",
        "LG50", "largest", "n_ge_1Mb", "frac_bp_ge_1Mb", "n_ge_10Mb", "frac_bp_ge_10Mb", "n_ge_50Mb",
        "frac_bp_ge_50Mb", "records_anchor_bearing", "bp_anchor_bearing", "bp_in_anchor_blocks", "bp_unanchored",
        "records_per_chr_median", "records_per_chr_max", "ref_bp_covered_by_blocks", "gene_anchors", "busco_lineage",
        "busco_C", "busco_S", "busco_D", "busco_F", "busco_M", "busco_n", "merqury_qv", "merqury_completeness",
        "merqury_pair_completeness", "notes"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "run"])
    ap.add_argument("--assemblies", required=True)
    ap.add_argument("--genome-size", type=float, default=0, help="for NG50 / LG50 (bp)")
    ap.add_argument("--out")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    specs = read_tsv(a.assemblies)
    if len({s["name"] for s in specs}) != len(specs):
        die("assembly names must be unique")
    rows, lens_by, inputs = [], {}, []
    for s in specs:
        miss = [k for k in REQ if not s.get(k)]
        if miss:
            die(f"{s.get('name', '?')} lacks {', '.join(miss)}")
        for k in ("fai", "busco", "merqury_qv", "merqury_cn"):
            if s.get(k):
                if not os.path.isfile(s[k]):
                    die(f"{s['name']}: {k} not found: {s[k]}")
                inputs.append(dict(name=s["name"], kind=k, path=os.path.abspath(s[k]), sha256=sha256(s[k])))
        lens = lengths(s["fai"])
        lens_by[s["name"]] = lens
        r = {k: s.get(k, "") for k in ("line", "name", "assembly_class", "role", "notes")}
        r.update(contiguity(lens, a.genome_size))
        if s.get("busco"):
            r.update(busco(s["busco"]))
        if s.get("merqury_qv"):
            r["merqury_qv"] = round(merqury_row(s["merqury_qv"], s.get("merqury_qv_row"), 3), 2)
        if s.get("merqury_cn"):
            r["merqury_completeness"] = round(merqury_row(s["merqury_cn"], s.get("merqury_cn_row"), 4), 2)
            if s.get("merqury_cn_pair_row"):
                r["merqury_pair_completeness"] = round(merqury_row(s["merqury_cn"], s["merqury_cn_pair_row"], 4), 2)
        if s.get("riparian_dir"):
            r.update(anchored(s["riparian_dir"], s["riparian_unit"]))
        rows.append(r)
    if a.command == "check":
        print(f"check OK: {len(rows)} assemblies (nothing written)", file=sys.stderr)
        return 0
    if not a.out:
        die("run needs --out")
    if os.path.exists(a.out) or os.path.exists(a.out + ".partial"):
        die(f"{a.out} (or .partial) exists")
    start = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    part = a.out + ".partial"
    os.makedirs(part)
    with open(os.path.join(part, "assembly_compare.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLS})
    figs = [] if a.no_plots or not a.genome_size else plot(part, rows, lens_by, a.genome_size)
    with open(os.path.join(part, "METHODS.md"), "w") as fh:
        fh.write(f"""# Assembly comparison: methods ({start})

Generated by `g2f_assembly_compare.py` (code version {CODE_VERSION}). Inputs with sha256: RUN.json.

- **Contiguity** from the FASTA index only. N50/L50 use the assembly's own size; NG50/LG50 use the genome size
  {a.genome_size:,.0f} bp. A haplotype FASTA is one assembly row; a pair is never summed.
- **Size classes**: contigs >= 1 / 10 / 50 Mb and the fraction of assembly bp in them.
- **BUSCO, Merqury**: read from the existing files named per row (no recomputation). Merqury completeness of one
  haplotype is expected to be below that of a pair or a primary, because heterozygous k-mers are split between
  haplotypes; `merqury_pair_completeness` is the pair together.
- **Anchors** (when a riparian run is named): records with a drawn AnchorWave block, their bp, the union of drawn
  blocks on the record (`bp_in_anchor_blocks`), bp of records without drawn blocks, records per dominant reference
  chromosome, and the reference bp covered by drawn blocks.
- Fewer contigs is not by itself better: a purged or phased assembly can be more contiguous while carrying less
  sequence, and an unpurged primary carries extra allelic copies. Compare within `assembly_class`.
""")
    with open(os.path.join(part, "RUN.json"), "w") as fh:
        json.dump(dict(schema="g2f-assembly-compare/1", code_version=CODE_VERSION, start=start,
                       command=[sys.executable] + sys.argv, cwd=os.getcwd(), genome_size=a.genome_size,
                       assemblies=os.path.abspath(a.assemblies), assemblies_sha256=sha256(a.assemblies),
                       inputs=inputs, script_sha256=sha256(__file__),
                       outputs=["assembly_compare.tsv", "METHODS.md"] + [os.path.relpath(p, part) for p in figs]),
                  fh, indent=1)
    os.rename(part, a.out)
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
