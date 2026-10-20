"""Tests for the assembly-validation workflow (scripts/g2f_asmqc.py).

Synthetic data only. The end-to-end test runs the real minimap2 and samtools
on a 500 kb planted genome (~30 s) and is skipped if they are missing. Needs
numpy, pysam, matplotlib, networkx (BioHPC: the default /usr/bin/python3).
"""
import contextlib
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml

PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PIPELINE_DIR, "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import g2f_asmqc  # noqa: E402
from asmqc import graph as G  # noqa: E402
from asmqc import paf, pairalign, readsupport, regions, runner  # noqa: E402
import asmqc_data  # noqa: E402


def find(name, *alts):
    p = shutil.which(name)
    if p:
        return p
    for a in alts:
        if os.path.exists(a):
            return a
    return None


MINIMAP2 = find("minimap2", "/programs/minimap2-2.30/minimap2")
SAMTOOLS = find("samtools", "/programs/samtools-1.20/bin/samtools")
SEQKIT = find("seqkit", "/programs/seqkit-2.13.0/seqkit")
HAVE_TOOLS = all((MINIMAP2, SAMTOOLS, SEQKIT))


def tmpdir():
    base = os.environ.get("TMPDIR") or os.path.join(PIPELINE_DIR, "..", "tmp")
    os.makedirs(base, exist_ok=True)
    return tempfile.mkdtemp(prefix="asmqc_test_", dir=base)


class Units(unittest.TestCase):
    def test_md_events(self):
        # 10 matches, mismatch A, 5 matches, deletion of 2, 3 matches, mismatch C
        self.assertEqual(readsupport.md_events("10A5^GT3C0", 100), [110, 116, 121])

    def test_depth_mode_ignores_low_and_overflow(self):
        h = [1000, 900] + [0] * 28 + [5, 50, 100, 50, 5] + [0] * 965 + [10 ** 6]
        self.assertEqual(readsupport.depth_mode(h), 32)

    def test_chunks_aligned(self):
        ch = readsupport.chunks_for([("a", 25), ("b", 5)], 5, 10)
        self.assertEqual(ch, [("a", 25, 0, 10), ("a", 25, 10, 20), ("a", 25, 20, 25), ("b", 5, 0, 5)])

    def test_interval_ops(self):
        self.assertEqual(paf.merge_intervals([(5, 8), (1, 3), (2, 6)]), [[1, 8]])
        self.assertEqual(paf.subtract([[0, 10]], [[2, 3], [5, 7]]), [[0, 2], [3, 5], [7, 10]])
        self.assertEqual(paf.subtract([[0, 10]], []), [[0, 10]])

    def test_depth_classes(self):
        b = {"none": 0.15, "low": 0.3, "half_max": 0.7, "normal_max": 1.4, "dup_max": 2.5}
        self.assertEqual([regions.depth_class(x, b) for x in (0.1, 0.2, 0.5, 1.0, 2.0, 3.0)],
                         ["none", "low", "half", "normal", "dup", "multi"])

    def test_auto_depth_fails_on_one_allele_mode(self):
        # regression: an unpurged assembly of a heterozygous sample has its depth mode at the one-allele
        # peak (17) while the homozygous peak is 34; auto must refuse rather than classify against 17
        with self.assertRaises(regions.CalibrationError) as cm:
            regions.expected_depth("auto", 17, kmer_peak=34)
        self.assertIn("half the k-mer peak", str(cm.exception))
        with self.assertRaises(regions.CalibrationError):
            regions.expected_depth("auto", 45, kmer_peak=34)  # any >15% mismatch, not only half
        with self.assertRaises(regions.CalibrationError):
            regions.expected_depth("auto", 34, kmer_peak=None)  # unverifiable
        with self.assertRaises(regions.CalibrationError):
            regions.expected_depth(0, 34, kmer_peak=34)

    def test_auto_depth_accepts_matching_mode(self):
        c, src, warn = regions.expected_depth("auto", 33, kmer_peak=34)
        self.assertEqual((c, warn), (33.0, ""))
        self.assertIn("checked", src)

    def test_configured_depth_is_shared_baseline(self):
        # the configured baseline is used as given, whatever the assembly's own mode is
        self.assertEqual(regions.expected_depth(34, 17, kmer_peak=34), (34.0, "configured", ""))
        c, src, warn = regions.expected_depth(20, 17, kmer_peak=34)
        self.assertEqual(c, 20.0)
        self.assertIn("differs", warn)

    def test_shared_baseline_classes_unpurged_and_purged(self):
        sci = {"regions": {"depth_classes": {"none": 0.15, "low": 0.3, "half_max": 0.7, "normal_max": 1.4,
                                             "dup_max": 2.5}}}

        def classes(depths, c):
            n = len(depths)
            w = {"depth_all": np.array(depths, float), "depth_hq": np.array(depths, float),
                 "tested_bp": np.ones(n), "mixed_sites": np.zeros(n)}
            return list(regions.classify_windows(w, c, sci)[1])

        # unpurged: separately kept alleles at ~17x, homozygous sequence at ~33x, a true collapse at 68x
        unpurged = [17, 16, 18, 17, 33, 34, 32, 68]
        # purged: alleles collapsed to one copy (both haplotypes' reads) at ~34x, one kept allele pair at 17x
        purged = [34, 31, 35, 33, 17, 68]
        self.assertEqual(classes(unpurged, 34.0), ["half"] * 4 + ["normal"] * 3 + ["dup"])
        self.assertEqual(classes(purged, 34.0), ["normal"] * 4 + ["half", "dup"])
        # the defect being prevented: the unpurged mode (17) as baseline turns homozygous sequence into 'dup'
        self.assertEqual(classes(unpurged, 17.0)[4:7], ["dup"] * 3)

    def test_best_copy_pieces(self):
        alns = [{"q": "x", "t": "T", "ts": 0, "te": 100, "identity": 1.0, "alen": 100},
                {"q": "y", "t": "T", "ts": 50, "te": 150, "identity": 0.98, "alen": 100}]
        p = pairalign.best_copy_pieces(alns)
        best = {(i, s, e) for i, s, e, b in p if b}
        red = {(i, s, e) for i, s, e, b in p if not b}
        self.assertEqual(best, {(0, 0, 50), (0, 50, 100), (1, 100, 150)})
        self.assertEqual(red, {(1, 50, 100)})

    def test_target_to_query_reverse(self):
        a = {"qs": 0, "qe": 100, "ts": 1000, "te": 1100, "strand": "-"}
        self.assertEqual(pairalign.target_to_query(a, 1000, 1010), (90, 100))


GFA_CTG = """S\tptg1\t*\tLN:i:300\trd:i:10
A\tptg1\t0\t+\tr1\t0\t100\tid:i:1\tHG:A:a
A\tptg1\t100\t+\tr2\t0\t100\tid:i:2\tHG:A:a
A\tptg1\t200\t+\tr4\t0\t100\tid:i:4\tHG:A:a
S\tptg2\t*\tLN:i:100\trd:i:5
A\tptg2\t0\t+\tr3\t0\t100\tid:i:3\tHG:A:a
"""
GFA_UTG = """S\tutg1\t*\tLN:i:100\trd:i:10
S\tutg2\t*\tLN:i:100\trd:i:5
S\tutg3\t*\tLN:i:100\trd:i:5
S\tutg4\t*\tLN:i:100\trd:i:10
L\tutg1\t+\tutg2\t+\t0M
L\tutg1\t+\tutg3\t+\t0M
L\tutg2\t+\tutg4\t+\t0M
L\tutg3\t+\tutg4\t+\t0M
A\tutg1\t0\t+\tr1\t0\t100\tid:i:1
A\tutg2\t0\t+\tr2\t0\t100\tid:i:2
A\tutg3\t0\t+\tr3\t0\t100\tid:i:3
A\tutg4\t0\t+\tr4\t0\t100\tid:i:4
"""
GFA_CTG_PURGED = """S\tptg1\t*\tLN:i:300\trd:i:10
A\tptg1\t0\t+\tr1\t0\t100\tid:i:1\tHG:A:a
A\tptg1\t100\t+\tr2\t0\t100\tid:i:2\tHG:A:a
A\tptg1\t200\t+\tr4\t0\t100\tid:i:4\tHG:A:a
"""


class Graph(unittest.TestCase):
    def setUp(self):
        self.d = tmpdir()
        for name, text in (("ctg.gfa", GFA_CTG), ("utg.gfa", GFA_UTG), ("ctg_b.gfa", GFA_CTG_PURGED)):
            with open(os.path.join(self.d, name), "w") as fh:
                fh.write(text)

    def test_topology_and_paths(self):
        utg = G.parse_gfa(os.path.join(self.d, "utg.gfa"))
        t = G.topology(utg)
        self.assertEqual((t["segments"], t["links"], t["simple_bubbles"]), (4, 4, 1))
        self.assertEqual(t["tips"], 2)
        ctg = G.parse_gfa(os.path.join(self.d, "ctg.gfa"))
        rows, found = G.contig_paths(ctg, utg)
        self.assertEqual(found, 1.0)
        self.assertEqual([r[3] for r in rows if r[0] == "ptg1"], ["utg1", "utg2", "utg4"])
        trans = G.annotate_transitions(rows, utg)
        self.assertTrue(all(x["linked"] for x in trans))
        self.assertEqual(trans[0]["exit_degree"], 2)  # utg1 + end branches

    def test_neighbourhood_render_and_retention(self):
        utg = G.parse_gfa(os.path.join(self.d, "utg.gfa"))
        nodes = G.neighborhood(utg, {"utg2"}, 1, 10)
        self.assertEqual(nodes, {"utg1", "utg2", "utg4"})
        G.write_subgraph(utg, nodes, os.path.join(self.d, "sub.gfa"))
        G.render(utg, nodes, {"utg2"}, {"utg1", "utg2", "utg4"}, 10, "test", os.path.join(self.d, "sub"))
        self.assertTrue(os.path.getsize(os.path.join(self.d, "sub.svg")) > 1000)
        r = G.retention(os.path.join(self.d, "ctg.gfa"), os.path.join(self.d, "ctg_b.gfa"),
                        os.path.join(self.d, "utg.gfa"), self.d)
        self.assertEqual(r["b_reads_also_in_a"], 1.0)
        self.assertEqual(r["dropped_bp_by_utg_rd_class"], {"half": 100})


def write_config(d, paths, extra=None):
    cfg = {
        "profile": "test_v1", "outdir": os.path.join(d, "out"), "tmpdir": os.path.join(d, "tmp"),
        "science_overrides": asmqc_data.SCIENCE_OVERRIDES,
        "tools": {"minimap2": {"cmd": MINIMAP2, "version_cmd": f"{MINIMAP2} --version"},
                  "samtools": {"cmd": SAMTOOLS, "version_cmd": f"{SAMTOOLS} --version | head -1"},
                  "seqkit": {"cmd": SEQKIT, "version_cmd": f"{SEQKIT} version"}},
        "resources": {"threads": 2, "sort_threads": 1, "sort_mem_per_thread": "200M", "reserve_disk_gb": 1,
                      "min_free_ram_gb": 1, "nice": 0},
        "references": {"REF": {"fasta": paths["ref"], "chromosomes_regex": r"^chr[A-Z]$"}},
        "samples": {"SAMPLE_A": {
            "hifi_reads": [paths["reads"]], "hifi_read_count": paths["n_reads"], "expected_depth": 30,
            "assemblies": {"unpurged": {"fasta": paths["asm_a"], "label": "unpurged"},
                           "purged": {"fasta": paths["asm_b"], "label": "purged"}},
            "pairs": [{"a": "unpurged", "b": "purged", "use_b_extra_targets": False}]}},
    }
    if extra:
        extra(cfg)
    path = os.path.join(d, "asmqc.yaml")
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh)
    return path


def run_cli(*args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = g2f_asmqc.main(list(args))
    return rc, buf.getvalue()


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


@unittest.skipUnless(HAVE_TOOLS, "minimap2/samtools/seqkit not available")
class EndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tmpdir()
        cls.paths = asmqc_data.build(os.path.join(cls.d, "data"))
        cls.cfg = write_config(cls.d, cls.paths)
        cls.rc, cls.out = run_cli("-c", cls.cfg, "run")
        cls.root = os.path.join(cls.d, "out", "test_v1", "samples", "SAMPLE_A")

    def test_all_stages_validated(self):
        self.assertEqual(self.rc, 0, self.out)
        rc, out = run_cli("-c", self.cfg, "status")
        self.assertNotIn("planned", out)
        self.assertNotIn("failed", out)

    def test_rerun_skips(self):
        rc, out = run_cli("-c", self.cfg, "run")
        self.assertEqual(rc, 0)
        self.assertNotIn("[run ]", out)

    def test_read_count_and_depth(self):
        v = json.load(open(os.path.join(self.root, "unpurged", "map_hifi", "VALIDATED.json")))
        self.assertEqual(v["validation"]["primary"], self.paths["n_reads"])
        s = json.load(open(os.path.join(self.root, "purged", "read_support", "summary.json")))
        self.assertTrue(25 <= s["depth_mode_all"] <= 35, s["depth_mode_all"])

    def test_planted_features_found(self):
        reg = tsv(os.path.join(self.root, "unpurged", "regions", "regions.tsv"))
        calls = {(r["contig"], r["call"]) for r in reg}
        self.assertIn(("ctg1", "collapsed_repeat_or_duplicate"), calls)
        self.assertIn(("ctg3", "haploid_depth_region"), calls)
        self.assertIn(("ctg4", "near_identical_duplicate"), calls)
        col = [r for r in reg if r["contig"] == "ctg1" and r["call"] == "collapsed_repeat_or_duplicate"]
        self.assertTrue(any(int(r["start"]) <= 52000 and int(r["end"]) >= 54000 for r in col), col)
        bps = tsv(os.path.join(self.root, "unpurged", "regions", "breakpoints.tsv"))
        mis = [b for b in bps if b["contig"] == "ctg2" and b["call"] == "suspected_misjoin"]
        self.assertTrue(any(abs(int(b["pos"]) - 150000) < 200 for b in mis), bps)
        # the only other misjoin is the junction the collapsed second repeat copy leaves at ctg1:120 kb
        other = [b for b in bps if b["call"] == "suspected_misjoin" and b["contig"] != "ctg2"]
        self.assertEqual([(b["contig"], abs(int(b["pos"]) - 120000) < 200) for b in other], [("ctg1", True)])

    def test_reference_junction_and_evidence(self):
        j = tsv(os.path.join(self.root, "unpurged", "ref_REF_analysis", "junctions.tsv"))
        self.assertTrue(any(x["contig"] == "ctg2" and x["kind"] == "inter_chromosomal" for x in j), j)
        cand = tsv(os.path.join(self.root, "unpurged", "evidence", "candidates.tsv"))
        top = cand[0]  # reads + reference agree: ranked first
        self.assertEqual((top["contig"], top["category"], top["confidence"], top["n_evidence_types"]),
                         ("ctg2", "suspected_misassembly", "high", "2"))
        c1 = [c for c in cand if c["contig"] == "ctg1" and c["category"] == "suspected_misassembly"]
        self.assertIn("missing repeat copy", c1[0]["read_evidence"])

    def test_pair_classes(self):
        s = json.load(open(os.path.join(self.root, "pairs", "unpurged_vs_purged", "analysis", "summary.json")))
        bp = s["class_bp"]
        self.assertGreater(bp["redundant_identical"], 25000)   # ctg4
        self.assertGreater(bp["redundant_divergent"], 35000)   # ctg3 (1.5% divergent)
        self.assertLess(bp["A_unique"], 2000)
        self.assertLess(s["b_only_bp"], 2000)
        x = s["depth_crosstab_bp"]
        # each retained allele copy is at half depth; the purged assembly's single copy at full depth
        self.assertGreater(x.get("redundant_divergent|A:half|B:normal", 0), 30000)
        self.assertGreater(x.get("redundant_identical|A:half|B:normal", 0), 25000)

    def test_report(self):
        rc, out = run_cli("-c", self.cfg, "report")
        rep = out.strip().splitlines()[-1]
        self.assertTrue(os.path.exists(os.path.join(rep, "report.html")))
        md = open(os.path.join(rep, "report.md")).read()
        self.assertIn("suspected_misassembly", md)
        m = tsv(os.path.join(rep, "metrics.tsv"))
        self.assertTrue(any(r["metric"] == "N50" for r in m))
        self.assertTrue(os.path.exists(os.path.join(rep, "figures", "depth_histograms.svg")))


@unittest.skipUnless(HAVE_TOOLS, "minimap2/samtools/seqkit not available")
class RunnerSafety(unittest.TestCase):
    def setUp(self):
        self.d = tmpdir()
        self.paths = asmqc_data.build(os.path.join(self.d, "data"), depth_per_hap=2)

    def test_plan_runs_nothing(self):
        cfg = write_config(self.d, self.paths)
        rc, out = run_cli("-c", cfg, "plan")
        self.assertIn("plan only", out)
        self.assertFalse(os.path.exists(os.path.join(self.d, "out", "test_v1", "samples")))

    def test_failure_is_recorded_and_preserved(self):
        def broken(c):
            c["tools"]["minimap2"] = {"cmd": "false", "version_cmd": "echo broken"}
        cfg = write_config(self.d, self.paths, broken)
        rc, out = run_cli("-c", cfg, "run", "--stages", "prepare,map_hifi", "--samples", "SAMPLE_A")
        self.assertEqual(rc, 1)
        # stages run in a fixed order (purged before unpurged) and the run stops at the first failure
        part = os.path.join(self.d, "out", "test_v1", "samples", "SAMPLE_A", "purged", "map_hifi.partial")
        self.assertTrue(os.path.exists(os.path.join(part, "FAILED.json")))
        rc, out = run_cli("-c", cfg, "status")
        self.assertIn("failed", out)
        rc, out = run_cli("-c", cfg, "run", "--stages", "map_hifi")
        attempts = [x for x in os.listdir(os.path.dirname(part)) if x.startswith("map_hifi.attempt-")]
        self.assertEqual(len(attempts), 1)  # first attempt kept, not deleted

    def test_stale_input_is_not_overwritten(self):
        cfg = write_config(self.d, self.paths)
        rc, _ = run_cli("-c", cfg, "run", "--stages", "prepare")
        self.assertEqual(rc, 0)
        os.utime(self.paths["asm_b"], (1, 1))
        rc, out = run_cli("-c", cfg, "status")
        self.assertIn("stale", out)
        rc, out = run_cli("-c", cfg, "run", "--stages", "prepare")
        self.assertEqual(rc, 1)
        rc, out = run_cli("-c", cfg, "run", "--stages", "prepare", "--rerun-stale")
        self.assertEqual(rc, 0)
        prep = os.path.join(self.d, "out", "test_v1", "samples", "SAMPLE_A", "purged")
        self.assertTrue([x for x in os.listdir(prep) if x.startswith("prepare.previous-")])

    def test_depth_calibration_fails_closed_then_reruns_only_downstream(self):
        root = os.path.join(self.d, "out", "test_v1", "samples", "SAMPLE_A")

        def setting(expected, kmer):
            def f(c):
                c["samples"]["SAMPLE_A"]["expected_depth"] = expected
                c["samples"]["SAMPLE_A"]["kmer_peak"] = kmer
            return write_config(self.d, self.paths, f)

        # auto without a k-mer peak: regions refuses to classify, the failure is recorded
        cfg = setting("auto", None)
        rc, out = run_cli("-c", cfg, "run", "--stages", "prepare,map_hifi,read_support,regions")
        self.assertEqual(rc, 1, out)
        with open(os.path.join(root, "purged", "regions.partial", "FAILED.json")) as fh:
            failed = json.load(fh)
        self.assertIn("CalibrationError", failed["error"])
        self.assertFalse(os.path.exists(os.path.join(root, "purged", "regions", "VALIDATED.json")))
        m = json.load(open(os.path.join(root, "unpurged", "read_support", "summary.json")))["depth_mode_all"]

        # auto with the read-depth mode at half the k-mer peak (the 17x vs 34x case)
        cfg = setting("auto", 2 * m)
        rc, out = run_cli("-c", cfg, "run", "--stages", "regions")
        self.assertEqual(rc, 1, out)
        self.assertIn("half the k-mer peak", out)

        # explicit shared baseline: everything validates; the failed attempts are kept
        cfg = setting(2 * m, 2 * m)
        rc, out = run_cli("-c", cfg, "run")
        self.assertEqual(rc, 0, out)
        self.assertNotIn("[run ] map_hifi", out)  # validated alignments are reused
        self.assertTrue([x for x in os.listdir(os.path.join(root, "purged")) if x.startswith("regions.attempt-")])
        for a in ("purged", "unpurged"):
            s = json.load(open(os.path.join(root, a, "regions", "summary.json")))
            self.assertEqual((s["expected_depth"], s["expected_depth_source"]), (2.0 * m, "configured"))

        # changing the baseline makes regions stale by parameters; its dependants follow by input signature
        # once regions is rerun (status does not propagate upstream staleness)
        cfg = setting(m, 2 * m)
        rc, out = run_cli("-c", cfg, "status")
        stale = {l.split()[1] for l in out.splitlines() if l.split()[:1] == ["stale"]}
        self.assertEqual({x.split(":")[0] for x in stale}, {"regions"}, out)
        rc, out = run_cli("-c", cfg, "run", "--rerun-stale")
        self.assertEqual(rc, 0, out)
        ran = {l.split()[2].split(":")[0] for l in out.splitlines() if l.startswith("[run ]")}
        self.assertEqual(ran, {"regions", "evidence", "pair_analysis"}, out)
        self.assertTrue([x for x in os.listdir(os.path.join(root, "purged")) if x.startswith("regions.previous-")])

    def test_wrong_read_count_fails_validation(self):
        def wrong(c):
            c["samples"]["SAMPLE_A"]["hifi_read_count"] = 1
        cfg = write_config(self.d, self.paths, wrong)
        rc, out = run_cli("-c", cfg, "run", "--stages", "prepare,map_hifi")
        self.assertEqual(rc, 1)
        self.assertIn("expected HiFi read count", out)


if __name__ == "__main__":
    unittest.main()
