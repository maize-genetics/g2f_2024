"""Tests for g2f_telomere (synthetic sequences only; tidk is the stub in tests/stubs unless G2F_TIDK names a real tidk >= 0.3)."""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

import util
import telomere_data as TD
import g2f_telomere as T


def wtext(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rows_text(path):
    with open(path) as fh:
        return fh.readlines()

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_telomere.py")
REAL_TIDK = os.environ.get("G2F_TIDK")
MINIMAP2 = os.environ.get("G2F_MINIMAP2") or shutil.which("minimap2")


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(*args, env=None):
    e = dict(os.environ, **(env or {}))
    return subprocess.run([sys.executable, SCRIPT] + [str(x) for x in args], capture_output=True, text=True, env=e)


class Defaults(unittest.TestCase):
    def test_yaml_equals_builtin_defaults(self):
        with open(os.path.join(util.PIPELINE_DIR, "config", "telomere_science.yaml")) as fh:
            self.assertEqual(yaml.safe_load(fh), T.DEFAULTS)


class ArrayCaller(unittest.TestCase):
    def test_chains_tolerate_a_diverged_unit_and_split_at_long_gaps(self):
        s = "A" * 30 + "TTTAGGG" * 40 + "TTAAGGG" + "TTTAGGG" * 40 + "C" * 200 + "TTTAGGG" * 10
        ch = T.chains(s, "TTTAGGG", 50, 5)
        self.assertEqual(ch[0], (30, 30 + 7 * 81, 80))
        self.assertEqual(ch[1][2], 10)
        self.assertEqual(T.chains(s, "TTTAGGG", 50, 11), ch[:1])

    def test_scan_arrays_coordinates_strands_gaps(self):
        d = tempfile.mkdtemp()
        seq = "CCCTAAA" * 300 + "ACGT" * 2500 + "N" * 100 + "ACGT" * 2500 + "TTTAGGG" * 200
        fa = os.path.join(d, "x.fa")
        wtext(fa, ">s1 desc\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n>s2\nACGT\n")
        r = T.scan_arrays(fa, d, "x", "TTTAGGG", 50, 5)
        self.assertEqual((r["records"], r["bp"]), (2, len(seq) + 4))
        rows = tsv(os.path.join(d, "x.arrays.tsv"))
        self.assertEqual([(x["strand"], int(x["start0"]), int(x["end"]), int(x["hits"])) for x in rows],
                         [("C-rich", 0, 2100, 300), ("G-rich", len(seq) - 1400, len(seq), 200)])
        self.assertEqual(int(rows[1]["dist_end"]), 0)
        g = tsv(os.path.join(d, "x.gaps.tsv"))
        self.assertEqual([(int(x["start0"]), int(x["length"])) for x in g], [(12100, 100)])

    def test_palindromic_motif_refused(self):
        with self.assertRaises(SystemExit):
            T.scan_arrays(__file__, tempfile.mkdtemp(), "x", "ACGT", 50, 5)


class Pipeline(unittest.TestCase):
    """scan contigs + old scaffolds, integrate, and check classes, lift self-check and candidates."""

    @classmethod
    def setUpClass(cls):
        if not util.HAVE_SAMTOOLS:
            raise unittest.SkipTest("samtools needed")
        cls.d = tempfile.mkdtemp()
        cls.p = TD.build(cls.d, gz=False)
        cls.tidk = REAL_TIDK or TD.STUB_TIDK
        sci = cls.p["science"]
        for lab, fa in (("ctg", cls.p["contigs"]), ("old", cls.p["old_scaffolds"]), ("new", cls.p["new_expected_scaffolds"])):
            r = run("scan", "--fasta", fa, "--label", lab, "--out", os.path.join(cls.d, f"scan_{lab}"), "--science", sci,
                    "--tidk", cls.tidk)
            assert r.returncode == 0, r.stderr
        r = run("integrate", "--contigs-scan", os.path.join(cls.d, "scan_ctg"), "--scaffolds-scan", os.path.join(cls.d, "scan_old"),
                "--agp", cls.p["old_agp"], "--label", "UNIT", "--out", os.path.join(cls.d, "int_old"), "--science", sci)
        assert r.returncode == 0, r.stderr + r.stdout
        cls.integ = os.path.join(cls.d, "int_old")

    def test_scan_record(self):
        v = jload(os.path.join(self.d, "scan_ctg", "VALIDATED.json"))
        self.assertEqual(v["step"], "scan")
        self.assertEqual(v["arrays"]["records"], 7)
        self.assertEqual(v["fasta_sha256"], T.sha256_file(self.p["contigs"]))
        self.assertTrue(all(os.path.exists(os.path.join(self.d, "scan_ctg", f)) for f in v["outputs_sha256"]))
        self.assertTrue(os.path.exists(os.path.join(self.d, "scan_ctg", "search", "search_telomeric_repeat_windows.tsv.gz")))

    def test_scan_refuses_overwrite(self):
        r = run("scan", "--fasta", self.p["contigs"], "--label", "ctg", "--out", os.path.join(self.d, "scan_ctg"),
                "--science", self.p["science"], "--tidk", self.tidk)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("refusing to overwrite", r.stderr)

    def test_old_tidk_refused(self):
        d = tempfile.mkdtemp()
        old = os.path.join(d, "tidk")
        wtext(old, "#!/bin/sh\necho 'tidk 0.2.7'\n")
        os.chmod(old, 0o755)
        r = run("scan", "--fasta", self.p["contigs"], "--label", "x", "--out", os.path.join(d, "o"), "--tidk", old)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn(">= 0.3.0", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(d, "o")))

    def test_contig_end_classes(self):
        e = {(r["contig"], r["side"]): r for r in tsv(os.path.join(self.integ, "contig_ends.tsv"))}
        self.assertEqual(e[("ctgA", "L")]["class"], "TEL_KB")
        self.assertEqual(e[("ctgA", "R")]["class"], "ABSENT")
        self.assertEqual(e[("ctgT", "R")]["class"], "TEL_KB")
        self.assertEqual(e[("ctgU", "L")]["class"], "TEL_SHORT")
        self.assertEqual(e[("ctgU", "L")]["placed"], "False")

    def test_chromosome_ends_and_errors_found(self):
        ce = {(r["chromosome"], r["end"]): r for r in tsv(os.path.join(self.integ, "chromosome_ends.tsv"))}
        self.assertEqual(len(ce), 6)
        self.assertEqual(ce[("chr01", "start")]["class"], "TEL_KB")
        self.assertEqual(ce[("chr02", "start")]["class"], "ABSENT")
        self.assertEqual(ce[("chr02", "start")]["inward_contig_end_class"], "TEL_KB")    # buried tip telomere
        self.assertEqual(ce[("chr03", "end")]["class"], "ABSENT")
        cand = {(r["contig"], r["position_in_contig"], r["faces"]) for r in tsv(os.path.join(self.integ, "internal_candidates.tsv"))}
        self.assertIn(("ctgT", "terminal", "join"), cand)
        self.assertIn(("ctgE", "terminal", "join"), cand)
        self.assertIn(("ctgA", "interstitial", "scaffold_end"), cand)   # nearest contig end = chromosome start

    def test_lift_self_check_and_summary(self):
        s = jload(os.path.join(self.integ, "summary.json"))
        self.assertEqual(s["self_checks"]["contig_arrays_lifted"], s["self_checks"]["found_identically_in_scaffold_scan"])
        self.assertEqual(s["self_checks"]["scaffold_vs_contig_end_class_mismatches"], 0)
        self.assertEqual((s["telomere_kb"], s["non_telomeric"], s["interstitial_arrays"]), (4, 2, 1))
        self.assertEqual(s["unplaced_telomeric_contigs"], ["ctgU"])
        self.assertIn("ctgE", s["contigs_telomeric_both_ends"])

    def test_integrate_refuses_mismatched_agp(self):
        r = run("integrate", "--contigs-scan", os.path.join(self.d, "scan_ctg"), "--scaffolds-scan", os.path.join(self.d, "scan_old"),
                "--agp", self.p["new_expected_agp"], "--label", "UNIT", "--out", os.path.join(self.d, "int_bad"),
                "--science", self.p["science"])
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(os.path.exists(os.path.join(self.d, "int_bad")))

    def test_compare_before_after(self):
        r = run("compare", "--before", os.path.join(self.d, "scan_old"), "--after", os.path.join(self.d, "scan_new"),
                "--label", "UNIT", "--out", os.path.join(self.d, "cmp"), "--science", self.p["science"])
        self.assertEqual(r.returncode, 0, r.stderr)
        s = jload(os.path.join(self.d, "cmp", "summary.json"))
        self.assertEqual((s["telomeric_ends_before"], s["telomeric_ends_after"], s["improved"], s["regressed"]), (4, 6, 2, 0))
        self.assertEqual((s["internal_kb_arrays_before"], s["internal_kb_arrays_after"], s["new_internal_arrays_after"]), (2, 0, 0))

    def test_closure(self):
        ctg = self.p["ctg"]
        clip_tel = "CCCTAAA" * 100
        clip_rand = TD.contigs(seed=99)["ctgX"][:700]
        sam = ["@HD\tVN:1.6\tSO:coordinate"] + [f"@SQ\tSN:{k}\tLN:{len(v)}\tM5:{TD.md5u(v)}" for k, v in ctg.items()]
        a = ctg["ctgA"][:3000]
        sam.append(f"r1\t0\tctgA\t1\t60\t700S3000M\t*\t0\t0\t{clip_tel}{a}\t*")
        b = ctg["ctgB"][-3000:]
        sam.append(f"r2\t0\tctgB\t{len(ctg['ctgB']) - 2999}\t60\t3000M700S\t*\t0\t0\t{b}{clip_rand}\t*")
        e = ctg["ctgE"][-2000:]
        sam.append(f"r3\t0\tctgE\t{len(ctg['ctgE']) - 1999}\t60\t2000M\t*\t0\t0\t{e}\t*")
        sam_p = os.path.join(self.d, "r.sam")
        wtext(sam_p, "\n".join(sam) + "\n")
        bam = os.path.join(self.d, "r.bam")
        subprocess.run([util.SAMTOOLS, "sort", "-o", bam, sam_p], check=True, capture_output=True)
        subprocess.run([util.SAMTOOLS, "index", bam], check=True)
        r = run("closure", "--integrated", self.integ, "--alignments", bam, "--contigs", self.p["contigs"],
                "--out", os.path.join(self.d, "clo"), "--science", self.p["science"])
        self.assertEqual(r.returncode, 0, r.stderr)
        c = {(x["contig"], x["side"]): x["status"] for x in tsv(os.path.join(self.d, "clo", "closure.tsv"))}
        self.assertEqual(c[("ctgA", "L")], "CLIPPED_TELOMERIC")
        self.assertEqual(c[("ctgB", "R")], "CLIPPED_NON_TELOMERIC")
        self.assertEqual(c[("ctgE", "R")], "CLOSED")
        self.assertEqual(c[("ctgD", "R")], "NO_READS")

    def test_plot(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not installed")
        r = run("plot", "--integrated", self.integ, "--out", os.path.join(self.d, "plots"))
        self.assertEqual(r.returncode, 0, r.stderr)
        for f in ("map_UNIT.pdf", "map_UNIT.svg", "end_heatmap.pdf", "VALIDATED.json"):
            self.assertTrue(os.path.exists(os.path.join(self.d, "plots", f)), f)

    @unittest.skipUnless(MINIMAP2, "minimap2 not available (set G2F_MINIMAP2)")
    def test_tips_against_corrected_reference(self):
        # reference = the corrected layout: the buried chr02 tip must be INTERIOR_or_INVERTED, ctgT's telomeric end REF_START
        ref = os.path.join(self.d, "ref.fa")
        seqs = {k: v for k, v in self.p["new_expected_seqs"].items() if "_chr" in k}
        wtext(ref, "".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        r = run("tips", "--integrated", self.integ, "--contigs", self.p["contigs"], "--scaffolds", self.p["old_scaffolds"],
                "--reference", ref, "--out", os.path.join(self.d, "tips"), "--science", self.p["science"],
                "--minimap2", MINIMAP2, "--threads", "2")
        self.assertEqual(r.returncode, 0, r.stderr)
        t = {x["segment"]: x for x in tsv(os.path.join(self.d, "tips", "tip_calls.tsv"))}
        self.assertEqual(t["ctg|ctgT|R"]["outer_call"], "REF_START")
        self.assertEqual(t["scf|UNIT_chr02|L"]["outer_call"], "INTERIOR_or_INVERTED")


class TipCall(unittest.TestCase):
    def test_overlap_window_and_calls(self):
        # PAF rows: q, ql, qs, qe, strand, t, tl, ts, te, match, block, mapq
        row = lambda qs, qe, st, t, ts, te: ["q", "1000000", str(qs), str(qe), st, t, "200000000", str(ts), str(te), "1", "1", "60"]
        start = [row(0, 900000, "+", "chr2", 1000, 901000)]
        self.assertEqual(T.tip_call(start, {"chr2": 200000000}, 0, 300000, 3000000, 5)["call"], "REF_START")
        # one long alignment beginning in the outer window still counts for the inner window (overlap, not start)
        self.assertEqual(T.tip_call(start, {"chr2": 200000000}, 300000, 1000000, 3000000, 5)["call"], "REF_START")
        inverted = [row(0, 700000, "-", "chr2", 10000, 710000)]
        self.assertEqual(T.tip_call(inverted, {"chr2": 200000000}, 0, 300000, 3000000, 5)["call"], "INTERIOR_or_INVERTED")
        end = [row(0, 500000, "-", "chr2", 199400000, 199900000)]
        self.assertEqual(T.tip_call(end, {"chr2": 200000000}, 0, 300000, 3000000, 5)["call"], "REF_END")
        self.assertEqual(T.tip_call([], {}, 0, 300000, 3000000, 5)["call"], "NO_ALIGNMENT")


if __name__ == "__main__":
    unittest.main()
