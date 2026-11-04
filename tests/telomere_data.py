"""Synthetic contigs/scaffolds with telomeric arrays for test_telomere and test_agp_revise (invented sequences only).

Layout (old = "v0.2", with two deliberate scaffold errors):
  UNIT_chr01 = ctgA(+) gap ctgB(+)   ctgA: C-rich telomere at its start + a 399-bp interstitial array; ctgB: G-rich end
  UNIT_chr02 = ctgT(+) gap ctgD(+)   ctgT: tip contig whose G-rich telomere is at its END -> faces the join (error 1)
  UNIT_chr03 = ctgE(+) gap ctgX(+)   ctgE: telomeres at both ends; ctgX appended beyond ctgE's end telomere (error 2)
  ctgU                               unplaced, short C-rich telomere at its start
Fix plan: reverse ctgT; unplace ctgX.
"""
import gzip
import hashlib
import os
import random
import shutil
import subprocess

import util

G, C = "TTTAGGG", "CCCTAAA"
RC = str.maketrans("ACGTN", "TGCAN")
BGZIP = os.environ.get("G2F_BGZIP") or shutil.which("bgzip") or "bgzip"
HAVE_TOOLS = util.HAVE_SAMTOOLS and shutil.which(BGZIP) is not None
STUB_TIDK = os.path.join(util.PIPELINE_DIR, "tests", "stubs", "tidk")

SCIENCE = """version: 1
motif: TTTAGGG
tidk:
  primary: {window: 3000, min_length: 200, max_offset: 1000, exact: false}
  length_window: 10000
  sensitivity:
    - {name: l500, window: 10000, min_length: 500, max_offset: 1000, exact: false}
  search_window: 2000
arrays: {max_gap: 50, min_hits: 5}
classes: {kb_bp: 1000, short_min_bp: 200, terminal_window: 5000, at_end: 1000}
tips: {segment_bp: 6000, outer_bp: 2000, inner_bp: [2000, 6000], near_end_bp: 3000, min_mapq: 0, preset: asm20}
closure: {min_mapq: 0, outer_bp: 1000, min_clip_bp: 500}
"""


def rc(s):
    return s.translate(RC)[::-1]


def contigs(seed=11):
    rng = random.Random(seed)
    r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    tel_c = lambda n: C * n
    tel_g = lambda n: G * n
    ctg = {
        "ctgA": tel_c(214) + r(18000) + G * 57 + r(20000),           # 1,498-bp start telomere; 399-bp ITS
        "ctgB": r(28000) + tel_g(172),                                # 1,204-bp end telomere
        "ctgT": r(8000) + tel_g(200),                                 # tip contig, 1,400-bp telomere at its END
        "ctgD": r(29000) + tel_g(143),                                # 1,001-bp end telomere
        "ctgE": tel_c(186) + r(24000) + tel_g(157),                   # telomeres at both ends
        "ctgX": r(6000),                                              # appended beyond ctgE
        "ctgU": tel_c(43) + r(5000),                                  # 301-bp start telomere, unplaced
    }
    return ctg


OLD_LAYOUT = [("UNIT_chr01", [("ctgA", "+"), ("ctgB", "+")]),
              ("UNIT_chr02", [("ctgT", "+"), ("ctgD", "+")]),
              ("UNIT_chr03", [("ctgE", "+"), ("ctgX", "+")]),
              ("ctgU", [("ctgU", "+")])]
NEW_LAYOUT = [("UNIT_chr01", [("ctgA", "+"), ("ctgB", "+")]),
              ("UNIT_chr02", [("ctgT", "-"), ("ctgD", "+")]),
              ("UNIT_chr03", [("ctgE", "+")]),
              ("ctgU", [("ctgU", "+")]),
              ("ctgX", [("ctgX", "+")])]
PLAN = [dict(action="REVERSE_ORIENTATION", object="UNIT_chr02", component="ctgT", old_orientation="+", new_orientation="-",
             new_status="PLACED", evidence_class="TELOMERE_AND_SYNTENY", evidence="synthetic"),
        dict(action="UNPLACE_TERMINAL_EXTENSION", object="UNIT_chr03", component="ctgX", old_orientation="", new_orientation="",
             new_status="UNPLACED_BEYOND_TELOMERE", evidence_class="TELOMERE", evidence="synthetic")]


def layout_agp_and_seqs(ctg, layout, gap=100):
    lines, seqs = ["##agp-version 2.1"], {}
    for obj, comps in layout:
        pos, part, s = 1, 1, []
        for i, (c, o) in enumerate(comps):
            if i:
                lines.append(f"{obj}\t{pos}\t{pos + gap - 1}\t{part}\tU\t{gap}\tscaffold\tyes\tmap")
                s.append("N" * gap)
                pos += gap
                part += 1
            L = len(ctg[c])
            lines.append(f"{obj}\t{pos}\t{pos + L - 1}\t{part}\tW\t{c}\t1\t{L}\t{o}")
            s.append(ctg[c] if o == "+" else rc(ctg[c]))
            pos += L
            part += 1
        seqs[obj] = "".join(s)
    return "\n".join(lines) + "\n", seqs


def write_fasta(path, seqs, gz=False, width=60):
    text = "".join(f">{k}\n" + "\n".join(v[i:i + width] for i in range(0, len(v), width)) + "\n" for k, v in seqs.items())
    if gz:
        raw = path[:-3]
        with open(raw, "w") as fh:
            fh.write(text)
        subprocess.run([BGZIP, "-f", raw], check=True)
        subprocess.run([util.SAMTOOLS, "faidx", path], check=True)
    else:
        with open(path, "w") as fh:
            fh.write(text)
        if util.HAVE_SAMTOOLS:
            subprocess.run([util.SAMTOOLS, "faidx", path], check=True)
    return path


def build(d, gz=True):
    """Write contigs, old/new AGPs and scaffold FASTAs, plan and science into d. Returns dict of paths."""
    ctg = contigs()
    ext = ".fa.gz" if gz else ".fa"
    p = {"dir": d, "ctg": ctg}
    p["contigs"] = write_fasta(os.path.join(d, "contigs" + ext), ctg, gz)
    for tag, lay in (("old", OLD_LAYOUT), ("new_expected", NEW_LAYOUT)):
        agp, seqs = layout_agp_and_seqs(ctg, lay)
        with open(os.path.join(d, f"{tag}.agp"), "w") as fh:
            fh.write(agp)
        p[f"{tag}_agp"] = os.path.join(d, f"{tag}.agp")
        p[f"{tag}_scaffolds"] = write_fasta(os.path.join(d, f"{tag}_scaffolds" + ext), seqs, gz)
        p[f"{tag}_seqs"] = seqs
    keys = list(PLAN[0])
    with open(os.path.join(d, "plan.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n" + "".join("\t".join(r[k] for k in keys) + "\n" for r in PLAN))
    p["plan"] = os.path.join(d, "plan.tsv")
    with open(os.path.join(d, "science.yaml"), "w") as fh:
        fh.write(SCIENCE)
    p["science"] = os.path.join(d, "science.yaml")
    return p


def md5u(s):
    return hashlib.md5(s.upper().encode()).hexdigest()
