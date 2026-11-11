"""Tests for g2f_inversion_riparian.py + g2f_inversion_riparian_verify.py (synthetic tables only; needs matplotlib)."""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import util

FIG = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_inversion_riparian.py")
VER = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_inversion_riparian_verify.py")
try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

# One locus E1 on chr3: reversed (inner) block 400-600 kb, breakpoint intervals 380-400 and 600-620 kb.
# SAMPLE_A        B73 arrangement, contig forward
# SAMPLE_B_hap1   inverted (inner segment on '-'), two junctions
# SAMPLE_B_hap2   B73 arrangement, contig forward
# SAMPLE_C        B73 arrangement, contig runs opposite to the reference (drawn "rc")
SEG = [  # unit, contig, zone, strand, b0, b1, c0, c1
    ("SAMPLE_A", "ctgA", "left", "+", 100000, 400000, 200000, 500000),
    ("SAMPLE_A", "ctgA", "mid", "+", 400000, 600000, 500000, 700000),
    ("SAMPLE_A", "ctgA", "right", "+", 600000, 900000, 700000, 1000000),
    ("SAMPLE_B_hap1", "ctgB1", "left", "+", 100000, 400000, 50000, 350000),
    ("SAMPLE_B_hap1", "ctgB1", "inner", "-", 400000, 600000, 350000, 550000),
    ("SAMPLE_B_hap1", "ctgB1", "right", "+", 600000, 800000, 550000, 750000),
    ("SAMPLE_B_hap2", "ctgB2", "left", "+", 100000, 400000, 100000, 400000),
    ("SAMPLE_B_hap2", "ctgB2", "mid", "+", 400000, 600000, 400000, 600000),
    ("SAMPLE_B_hap2", "ctgB2", "right", "+", 600000, 900000, 600000, 900000),
    ("SAMPLE_C", "ctgC", "left", "-", 100000, 400000, 600000, 900000),
    ("SAMPLE_C", "ctgC", "mid", "-", 400000, 600000, 400000, 600000),
    ("SAMPLE_C", "ctgC", "right", "-", 600000, 900000, 100000, 400000),
]
CLEN = {"SAMPLE_A": ("ctgA", 1200000), "SAMPLE_B_hap1": ("ctgB1", 800000), "SAMPLE_B_hap2": ("ctgB2", 1500000),
        "SAMPLE_C": ("ctgC", 1100000)}
ASM = [("SAMPLE_A", "SAMPLE_A", ""), ("SAMPLE_B_hap1", "SAMPLE_B", "hap1"), ("SAMPLE_B_hap2", "SAMPLE_B", "hap2"),
       ("SAMPLE_C", "SAMPLE_C", "")]
STATE = {"SAMPLE_A": ("B73", ""), "SAMPLE_B_hap1": ("INVERTED", "flagged"), "SAMPLE_B_hap2": ("B73", ""),
         "SAMPLE_C": ("B73", "")}


def write(path, header, rows):
    with open(path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(d):
    cfg, inp = f"{d}/config", f"{d}/inputs"
    os.makedirs(cfg)
    os.makedirs(inp)
    with open(f"{inp}/ref.fa.fai", "w") as fh:
        for c, n in [(f"chr{i}", 1500000 + 100000 * i) for i in range(1, 11)] + [("scaf_9", 5000)]:
            fh.write(f"{c}\t{n}\t0\t60\t61\n")
    for u, (c, n) in CLEN.items():
        with open(f"{inp}/{u}.fai", "w") as fh:
            fh.write(f"{c}\t{n}\t0\t60\t61\nother\t1000\t0\t60\t61\n")
    write(f"{inp}/chain.tsv", ["event", "unit", "contig", "zone", "strand", "b73_start", "b73_end", "contig_start",
                               "contig_end", "aligned_bp"],
          [("E1",) + s + (s[5] - s[4],) for s in SEG])
    write(f"{inp}/intervals.tsv", ["event", "b73_inner_start", "b73_inner_end", "b73_outer_start", "b73_outer_end"],
          [("E1", 400000, 600000, 380000, 620000)])
    write(f"{inp}/states.tsv", ["event", "unit", "state", "flag"], [("E1", u) + STATE[u] for u in STATE])
    write(f"{inp}/junctions.tsv", ["event", "unit", "b73_breakpoint", "contig_junction", "single_reads_crossing"],
          [("E1", "SAMPLE_B_hap1", "left", "349000-351000", 3), ("E1", "SAMPLE_B_hap1", "right", "549000-551000", 0)])
    write(f"{cfg}/sources.tsv", ["key", "path", "note"],
          [("chain_segments", f"{inp}/chain.tsv", ""), ("event_intervals", f"{inp}/intervals.tsv", ""),
           ("arrangement", f"{inp}/states.tsv", ""), ("junctions_chain", f"{inp}/junctions.tsv", ""),
           ("b73_fai", f"{inp}/ref.fa.fai", "")])
    write(f"{cfg}/loci.tsv", ["order", "event", "chrom", "window_start", "window_end", "tick_step", "heading", "note"],
          [(1, "E1", "chr3", 100000, 900000, 200000, "Chromosome 3", "")])
    write(f"{cfg}/assemblies.tsv", ["order", "unit", "line", "partition", "source", "contig_fai"],
          [(i + 1, u, ln, p, "synthetic", f"{inp}/{u}.fai") for i, (u, ln, p) in enumerate(ASM)])
    write(f"{cfg}/legend.tsv", ["row", "kind", "text"],
          [(1, "aligned", "aligned in reference orientation"), (1, "reversed", "inversion"), (1, "contig", "not aligned"),
           (2, "rc", "reverse-complemented contig")])


def run(script, *args):
    return subprocess.run([sys.executable, script] + [str(x) for x in args], capture_output=True, text=True)


@unittest.skipUnless(HAVE_MPL, "matplotlib needed")
class InversionRiparian(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        build(cls.d)
        cls.runs = {}
        for name, extra in (("genome", []), ("focus", ["--ideogram", "focus", "--headings", "hide", "--flag-tag", "hide",
                                                       "--tag-anchor", "contig", "--junction-marks", "hide"])):
            out = f"{cls.d}/{name}/fig"
            os.makedirs(f"{cls.d}/{name}")                         # the script does not create the --out directory
            r = run(FIG, "--root", cls.d, "--tables", f"{cls.d}/{name}/tables", "--out", out, *extra)
            assert r.returncode == 0, r.stdout + r.stderr
            cls.runs[name] = (out, f"{cls.d}/{name}/tables")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def verify(self, name, tables=None):
        out, tab = self.runs[name]
        return run(VER, "--root", self.d, "--tables", tables or tab, "--svg", out + ".svg")

    def test_outputs_and_verifier_pass(self):
        for name, (out, tab) in self.runs.items():
            for ext in ("svg", "pdf", "png"):
                self.assertGreater(os.path.getsize(f"{out}.{ext}"), 0, f"{name} {ext}")
            r = self.verify(name)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            for t in ("T1", "T2", "T3", "T4", "T5", "T6", "T7"):
                self.assertRegex(r.stdout, rf"{t}\t\d+ checked\t0 failed\tPASS")
        self.assertTrue(os.path.exists(f"{self.runs['focus'][1]}/guides.tsv"))
        self.assertFalse(os.path.exists(f"{self.runs['genome'][1]}/guides.tsv"))

    def test_tracks(self):
        trk = {r["unit"]: r for r in read(f"{self.runs['genome'][1]}/tracks.tsv")}
        self.assertEqual(set(trk), set(CLEN))
        for u, (c, n) in CLEN.items():
            self.assertEqual((trk[u]["contig"], int(trk[u]["contig_len"])), (c, n))
        self.assertEqual(trk["SAMPLE_C"]["strand_drawn"], "rc")
        self.assertEqual({trk[u]["strand_drawn"] for u in ("SAMPLE_A", "SAMPLE_B_hap1", "SAMPLE_B_hap2")}, {"fwd"})
        self.assertEqual(trk["SAMPLE_B_hap1"]["contig_ends_in_window"], "850000")   # 800 kb contig, offset +50 kb
        self.assertEqual(len(trk["SAMPLE_B_hap1"]["junctions"].split(";")), 2)

    def test_ribbons(self):
        rib = read(f"{self.runs['genome'][1]}/ribbons.tsv")
        pairs = {(r["upper"], r["lower"]) for r in rib}
        self.assertEqual(pairs, {("B73", "SAMPLE_A"), ("SAMPLE_A", "SAMPLE_B_hap1"), ("SAMPLE_B_hap1", "SAMPLE_B_hap2"),
                                 ("SAMPLE_B_hap2", "SAMPLE_C")})
        for r in rib:
            mid = (int(r["b73_start"]) + int(r["b73_end"])) / 2
            in_inner = 400000 <= mid <= 600000 and "SAMPLE_B_hap1" in (r["upper"], r["lower"])
            self.assertEqual(r["inner"], "1" if in_inner else "0", r)
            self.assertEqual(r["crossing"], "1" if in_inner else "0", r)   # rc of SAMPLE_C restores B73 orientation

    def test_tags(self):
        gen = {r["unit"]: float(r["tag_x_mm"]) for r in read(f"{self.runs['genome'][1]}/tags.tsv")}
        foc = {r["unit"]: float(r["tag_x_mm"]) for r in read(f"{self.runs['focus'][1]}/tags.tsv")}
        self.assertEqual(set(gen), {"SAMPLE_C", "SAMPLE_B_hap1"})      # rc tag, and the flag tag (shown by default)
        self.assertEqual(set(foc), {"SAMPLE_C"})                       # --flag-tag hide: the flag alone draws no tag
        with open(self.runs["genome"][0] + ".svg") as fh:
            self.assertIn("†", fh.read())
        with open(self.runs["focus"][0] + ".svg") as fh:
            self.assertNotIn("†", fh.read())

    def test_deterministic_tables(self):
        tab = f"{self.d}/rerun"
        r = run(FIG, "--root", self.d, "--tables", tab, "--out", f"{self.d}/rerun_fig")
        self.assertEqual(r.returncode, 0, r.stderr)
        for t in ("ribbons.tsv", "tracks.tsv", "tags.tsv"):
            with open(f"{tab}/{t}") as a, open(f"{self.runs['genome'][1]}/{t}") as b:
                self.assertEqual(a.read(), b.read(), t)

    def test_verifier_detects_tampering(self):
        src = self.runs["genome"][1]
        bad = f"{self.d}/tampered"
        shutil.copytree(src, bad)
        rows = read(f"{bad}/tracks.tsv")
        rows[0]["contig_len"] = "999"
        rib = read(f"{bad}/ribbons.tsv")
        for r in rib:
            if r["inner"] == "1":
                r["crossing"] = "0"
        for name, rr in (("tracks.tsv", rows), ("ribbons.tsv", rib)):
            write(f"{bad}/{name}", list(rr[0]), [list(x.values()) for x in rr])
        r = self.verify("genome", tables=bad)
        self.assertEqual(r.returncode, 1, r.stdout)
        self.assertRegex(r.stdout, r"T1\t\d+ checked\t1 failed\tFAIL")
        self.assertRegex(r.stdout, r"T4\t\d+ checked\t[1-9]\d* failed\tFAIL")

    def test_suggest_order(self):
        r = run(FIG, "--root", self.d, "--suggest-order")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("SAMPLE_B", r.stdout)
        self.assertIn("# pattern over ['E1']", r.stdout)


if __name__ == "__main__":
    unittest.main()
