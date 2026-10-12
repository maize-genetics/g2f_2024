"""Tests for ONT tar archives and the optional genotype fingerprinting module.

Unit tests use the standard library only. The end-to-end test sketches
synthetic reads with the real sourmash library (tool `fingerprint_python` in
config/tools/biohpc.yaml) and is skipped if that interpreter is unavailable.
"""
import gzip
import io
import json
import os
import random
import shlex
import shutil
import subprocess
import tarfile
import tempfile
import unittest

import yaml
from util import PIPELINE_DIR, SAMTOOLS, make_cfg, mkdirs, write_fastq

import fingerprint_data as fd
from g2f import archive, discovery
from g2f import fingerprint as fp

FP_SCRIPT = os.path.join(PIPELINE_DIR, "scripts", "g2f_fingerprint.py")


def fingerprint_python():
    with open(os.path.join(PIPELINE_DIR, "config", "tools", "biohpc.yaml")) as fh:
        cmd = yaml.safe_load(fh)["tools"]["fingerprint_python"]["cmd"]
    try:
        r = subprocess.run(shlex.split(cmd) + ["-c", "import sourmash"], capture_output=True)
    except OSError:
        return None
    return cmd if r.returncode == 0 else None


def _tar(path, members):
    with tarfile.open(path, "w") as tf:
        for name, data in members:
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return path


def _gz_fastq(records):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for name, seq in records:
            gz.write(f"@{name}\n{seq}\n+\n{'5' * len(seq)}\n".encode())
    return buf.getvalue()


ONT_HDR = ("{} runid=r1 read=1 ch=1 flow_cell_id=PAQ00000 sample_id=S1 "
           "basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.2.0")


class ArchiveTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_pass_members_streamed_in_natural_order_fail_excluded(self):
        p = _tar(os.path.join(self.tmp, "run_fastq.tar"), [
            ("./fastq_pass/x_pass_10.fastq.gz", _gz_fastq([("r10", "ACGT" * 5)])),
            ("./fastq_pass/x_pass_2.fastq.gz", _gz_fastq([("r2", "ACGT" * 5)])),
            ("./fastq_fail/x_fail_0.fastq.gz", _gz_fastq([("bad", "ACGT" * 5)])),
        ])
        lst = archive.list_members(p)
        self.assertEqual([m["name"] for m in lst["use"]],
                         ["./fastq_pass/x_pass_2.fastq.gz", "./fastq_pass/x_pass_10.fastq.gz"])
        self.assertEqual(len(lst["excluded"]), 1)
        out = io.BytesIO()
        n, nbytes = archive.stream_raw(p, out)
        self.assertEqual(n, 2)
        text = gzip.decompress(out.getvalue()).decode()
        self.assertEqual([line[1:] for line in text.splitlines()[::4]], ["r2", "r10"])
        self.assertEqual("".join(archive.open_text(p)), text)

    def test_unexpected_member_and_mixed_compression_are_errors(self):
        p = _tar(os.path.join(self.tmp, "a.tar"), [
            ("./fastq_pass/a.fastq.gz", _gz_fastq([("r", "ACGT")])),
            ("./sequencing_summary.txt", b"x\n")])
        self.assertTrue(archive.list_members(p)["errors"])
        p = _tar(os.path.join(self.tmp, "b.tar"), [
            ("./fastq_pass/a.fastq.gz", _gz_fastq([("r", "ACGT")])),
            ("./fastq_pass/b.fastq", b"@r\nACGT\n+\n5555\n")])
        self.assertTrue(archive.list_members(p)["errors"])

    def test_truncated_archive_detected(self):
        p = _tar(os.path.join(self.tmp, "c.tar"), [
            ("./fastq_pass/a.fastq.gz", _gz_fastq([(f"r{i}", "ACGT" * 500) for i in range(200)]))])
        lst = archive.list_members(p)
        with open(p, "rb") as fh:
            data = fh.read()
        with open(p, "wb") as fh:
            fh.write(data[: lst["use"][0]["offset"] + 1000])
        with self.assertRaises(archive.ArchiveError):
            archive.stream_raw(p, io.BytesIO(), members=lst)


class TarDiscoveryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_ont_tar_is_discovered_and_selects_hifi_ont_method(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(os.path.join(d["hifi"], "h.fastq.gz"),
                    [f"m84109_240101_000000_s1/{i}/ccs" for i in range(5)], [1000] * 5)
        rng = random.Random(1)
        _tar(os.path.join(d["ont"], "run1_fastq.tar"), [
            ("./fastq_pass/p_0.fastq.gz", _gz_fastq(
                [(ONT_HDR.format(f"0000000{i}-0000-0000-0000-000000000000"),
                  "".join(rng.choice("ACGT") for _ in range(3000))) for i in range(3)]))])
        inputs, issues = discovery.discover(make_cfg(self.tmp))
        self.assertFalse([i for i in issues if i.severity == "error"], issues)
        ont = [i for i in inputs["S1"] if i.tech == "ont"][0]
        self.assertEqual(ont.format, "tar")
        self.assertEqual(ont.identity["n_members"], 1)
        self.assertEqual(ont.identity["sampled_sample_ids"], ["S1"])
        method, _ = discovery.select_method("S1", inputs["S1"], make_cfg(self.tmp))
        self.assertIn("hifi-ont", method)

    def test_tar_in_hifi_is_rejected(self):
        d = mkdirs(self.tmp, "S1", ("hifi",))
        _tar(os.path.join(d["hifi"], "x.tar"), [("./a.fastq.gz", _gz_fastq([("r", "ACGT")]))])
        _, issues = discovery.discover(make_cfg(self.tmp))
        self.assertIn("UNSUPPORTED_FORMAT", [i.code for i in issues])


# --------------------------------------------------------------------------
# Fingerprint statistics on simulated hash sets (no sourmash needed)
# --------------------------------------------------------------------------

def _sim_counts(kmers, depth, rng, err_rate=0.0, next_err=[10 ** 12]):
    """Abundances of a genome's k-mer hashes sequenced at `depth` (Poisson-ish)."""
    out = {}
    for h in kmers:
        c = sum(1 for _ in range(int(depth * 2)) if rng.random() < 0.5)
        if c:
            out[h] = c
    for _ in range(int(len(kmers) * err_rate)):
        next_err[0] += 1
        out[next_err[0]] = 1
    return out


def _sim_query(kmers, sens, rng, others=(), n_err=0):
    q = {h for h in kmers if rng.random() < sens}
    for o, frac in others:
        q |= {h for h in o if rng.random() < sens * frac}
    q |= {rng.getrandbits(62) + 2 ** 62 for _ in range(n_err)}
    return q


class FingerprintStatsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(3)
        base = list(range(1, 20001))
        cls.geno = {}
        nxt = 100000
        for g in ("A", "B", "C"):
            ks = set(rng.sample(base, 17000))       # shared ancestry
            ks |= set(range(nxt, nxt + 2000))       # private sequence
            nxt += 10000
            cls.geno[g] = ks
        cls.rng = rng
        P = {"min_count_floor": 2, "max_copy": 3.0}
        refs = {}
        cls.runs = {}
        for g, ks in cls.geno.items():
            runs = [_sim_counts(ks, 6, rng, 0.2), _sim_counts(ks, 6, rng, 0.2)]
            cls.runs[g] = runs
            tot = {}
            for r in runs:
                fp.add_counts(tot, r)
            refs[g] = fp.Reference(g, tot, P, ["r1", "r2"])
        cls.P = P
        cls.panel = fp.Panel(refs)

    def test_solid_band_excludes_errors_and_repeats(self):
        band = self.panel.refs["A"].band
        self.assertGreaterEqual(band["lo"], 2)
        self.assertTrue(10 <= band["peak"] <= 14, band)
        self.assertNotIn(10 ** 12 + 1, self.panel.refs["A"].solid)

    def test_private_and_distinct_sets(self):
        for g in "ABC":
            self.assertGreater(len(self.panel.private[g]), 1500)
            for h in "ABC":
                if h != g:
                    self.assertGreaterEqual(self.panel.distinct_n[g][h], len(self.panel.private[g]))

    def test_peak_anchored_to_expected_depth(self):
        """At low depth a repeat peak must not be taken for the single-copy peak."""
        hist = {1: 50000, 2: 9000, 3: 4000, 4: 2500, 5: 1200, 6: 500, 7: 200}
        hist.update({c: 3000 for c in range(40, 50)})               # repeat peak
        self.assertGreaterEqual(fp.solid_band(hist)["peak"], 40)    # unanchored: repeat peak
        band = fp.solid_band(hist, expected_depth=3.5)
        self.assertLessEqual(band["peak"], 6)
        self.assertIn("low depth", band["note"])

    def test_with_variant_equals_fresh_panel(self):
        tot = dict(self.runs["B"][0])
        v = fp.Reference("B", tot, self.P, ["r1"], kind="loo")
        a = self.panel.with_variant("B", v)
        refs = dict(self.panel.refs, B=v)
        b = fp.Panel(refs)
        self.assertEqual(a.private, b.private)
        self.assertEqual(a.distinct_n, b.distinct_n)
        self.assertEqual(a.core, b.core)

    def _calib(self):
        controls = []
        for g in "ABC":
            v = fp.Reference(g, dict(self.runs[g][1]), self.P, ["r2"], kind="loo")
            vp = self.panel.with_variant(g, v)
            controls.append((g, fp.score(vp, set(self.runs[g][0]) & vp.universe, 3.29)))
        return fp.calibrate(controls, [], 0.0, 50)

    def test_classification_matches_mislabels_and_mixtures(self):
        cal = self._calib()
        self.assertTrue(cal["discriminating"], cal)
        rng = random.Random(9)
        G = self.geno

        def cls(q, exp):
            return fp.classify(fp.score(self.panel, q, 3.29), exp, cal, 50)

        r = cls(_sim_query(G["A"], 0.4, rng, n_err=30000), "A")
        self.assertEqual(r["classification"], "supported", r)
        r = cls(_sim_query(G["A"], 0.4, rng, n_err=30000), "C")     # mislabelled
        self.assertEqual(r["classification"], "potentially_mismatched", r)
        self.assertEqual(r["best"], "A")
        r = cls(_sim_query(G["B"], 0.4, rng, others=[(G["C"], 0.5)]), "B")   # mixture
        self.assertEqual(r["classification"], "ambiguous", r)
        self.assertIn("C", r["possible_admixture"])
        r = cls(_sim_query(G["B"], 0.4, rng), "")
        self.assertEqual((r["classification"], r["best"]), ("assigned", "B"))
        r = cls(_sim_query(G["B"], 0.0005, rng), "B")                # too little data
        self.assertEqual(r["classification"], "ambiguous", r)

    def test_overlapping_controls_mean_no_discrimination(self):
        cal = fp.calibrate([], [], 0.0, 50)
        self.assertFalse(cal["discriminating"])
        r = fp.classify(fp.score(self.panel, set(self.geno["A"]), 3.29), "A", cal, 50)
        self.assertEqual(r["classification"], "ambiguous")


# --------------------------------------------------------------------------
# End to end with real sourmash on synthetic reads
# --------------------------------------------------------------------------

FPY = fingerprint_python()


@unittest.skipUnless(FPY, "sourmash interpreter (tools.fingerprint_python) not available")
class FingerprintEndToEndTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base = os.environ.get("TMPDIR") or tempfile.gettempdir()
        cls.tmp = tempfile.mkdtemp(prefix="g2f_fp_", dir=base)
        cls.items = fd.build(os.path.join(cls.tmp, "data"),
                             samtools=SAMTOOLS if os.path.exists(SAMTOOLS) else None)
        recs = []
        for label, tech, run, path, _ in cls.items:
            stem = os.path.join(cls.tmp, "sk", f"{tech}.{label}.{run}")
            subprocess.run(shlex.split(FPY) + [
                FP_SCRIPT, "sketch", "--input", path, "--tech", tech, "--label", label,
                "--run-id", f"{label}_{run}", "--out-sig", stem + ".sig.gz",
                "--out-json", stem + ".json", "--ksizes", "21,31", "--scaled", "10",
                "--min-read-len", "1000", "--samtools", SAMTOOLS],
                check=True, capture_output=True)
            recs.append(stem + ".json")
        cls.out = os.path.join(cls.tmp, "out")
        subprocess.run(shlex.split(FPY) + [FP_SCRIPT, "analyze", "--records", *recs,
                                           "--outdir", cls.out, "--ksizes", "21,31",
                                           "--params", json.dumps({"genome_size": 300_000})],
                       check=True, capture_output=True)
        with open(os.path.join(cls.out, "assignments.tsv")) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh]
        cls.rows = {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_outputs_exist(self):
        for f in ("assignments.tsv", "similarity_matrix.tsv", "qc_summary.md",
                  "parameters_and_provenance.json", "k21/pairwise_long.tsv", "k31/controls.tsv",
                  "k21/calibration.json", "k21/panel.tsv"):
            self.assertTrue(os.path.exists(os.path.join(self.out, f)), f)
        with open(os.path.join(self.out, "parameters_and_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertTrue(all(i["sourmash_version"] for i in prov["inputs"]))
        self.assertTrue(prov["calibration"]["21"]["discriminating"])

    def test_correct_labels_supported(self):
        for q in ("A1_ontA1", "B1_ontB1", "C1_ontC1_ok"):
            self.assertEqual(self.rows[q]["classification"], "supported", self.rows[q])

    def test_mislabel_detected(self):
        r = self.rows["C1_ontC1_mislabel"]
        self.assertEqual((r["classification"], r["best"]), ("potentially_mismatched", "A1"))
        self.assertEqual(r["class_k31"], "potentially_mismatched")

    def test_mixture_off_panel_and_tiny_not_supported(self):
        self.assertEqual(self.rows["B2_ontB2_mix"]["classification"], "ambiguous")
        self.assertNotEqual(self.rows["A1_ontX"]["classification"], "supported")
        self.assertEqual(self.rows["B1_ontB1_tiny"]["classification"], "ambiguous")

    def test_pooled_runs_with_conflicting_genotypes_not_supported(self):
        self.assertNotEqual(self.rows["C1:all_ont_runs"]["classification"], "supported")
        self.assertEqual(self.rows["C1_ontC1_ok"]["consistent_best"], "False")

    def test_fail_members_excluded_from_sketch(self):
        with open(os.path.join(self.tmp, "sk", "ont.A1.ontA1.json")) as fh:
            rec = json.load(fh)
        with open(os.path.join(self.tmp, "sk", "ont.A1.ontX.json")) as fh:
            rec_x = json.load(fh)
        self.assertEqual(rec["format"], "tar")
        n_pass = sum(1 for line in archive.open_text(rec["input"]["path"])) // 4
        self.assertEqual(rec["sampling"]["reads_seen"], n_pass)
        self.assertGreater(rec_x["sampling"]["reads_used"], 0)

    def test_hifi_controls_supported(self):
        with open(os.path.join(self.out, "k21", "controls.tsv")) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh]
        cls = [dict(zip(rows[0], r))["classification"] for r in rows[1:]]
        self.assertTrue(cls and all(c == "supported" for c in cls), cls)


if __name__ == "__main__":
    unittest.main()
