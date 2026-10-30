"""Tests for scripts/g2f_allmaps.py (synthetic data only).

Unit tests always run. The end-to-end tests run the real HISAT2 and jcvi (ALLMAPS)
containers on the synthetic genome in allmaps_data.py (about a minute, 2 threads) and are
skipped unless ALLMAPS_HISAT2_IMAGE and ALLMAPS_JCVI_IMAGE name the .sif images.
"""
import contextlib
import gzip
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, HERE)

import g2f_allmaps as am  # noqa: E402
import allmaps_data  # noqa: E402

HISAT2 = os.environ.get("ALLMAPS_HISAT2_IMAGE")
JCVI = os.environ.get("ALLMAPS_JCVI_IMAGE")


def tmpdir():
    base = os.environ.get("TMPDIR")
    if not base or not os.path.isdir(base):
        raise unittest.SkipTest("TMPDIR must name an existing scratch directory (tests/run_tests.sh sets it)")
    return tempfile.mkdtemp(prefix="amtest_", dir=base)


def write_config(root, data, science_overrides=None):
    with open(os.path.join(HERE, "..", "config", "allmaps_science.yaml")) as fh:
        sci = yaml.safe_load(fh)
    sci["markers"]["lg_labels"] = {"1": "01", "2": "02", "3": "03"}
    sci["scaffold"].update(min_contig_bp=20000, min_unique_markers=10, min_distinct_positions=5,
                           conflict_min_markers=5, conflict_min_span_bp=5000)
    sci["allmaps"]["plot_format"] = None
    for sect, kv in (science_overrides or {}).items():
        sci[sect].update(kv)
    sp = os.path.join(root, "science.yaml")
    with open(sp, "w") as fh:
        yaml.safe_dump(sci, fh)
    cfg = {"profile": "amtest", "outdir": os.path.join(root, "out"), "tmpdir": os.path.join(root, "tmp"),
           "science": sp,
           "runtime": {"hisat2_image": HISAT2 or "/nonexistent_h.sif", "jcvi_image": JCVI or "/nonexistent_j.sif",
                       "bind": [root]},
           "resources": {"threads": 2, "cpus": "0-1", "nice": 10, "max_rss_gb": 16, "min_free_ram_gb": 1,
                         "reserve_disk_gb": 1, "max_cpu_busy_frac": 1.0},
           "markers": {"fasta": data["markers_fa"], "metadata": data["markers_txt"]},
           "samples": {"SAMPLE_A": {"fasta": data["asm_fa"], "label": "SAMPLE-A.v1"},  # hyphenated label: map name must stay parseable
                       "REF": {"fasta": data["ref_fa"], "label": "REF"}},
           "reference_positions": {"sample": "REF", "lg_to_seq": {"1": "R1", "2": "R2", "3": "R3"}}}
    p = os.path.join(root, "config.yaml")
    with open(p, "w") as fh:
        yaml.safe_dump(cfg, fh)
    return p


def load(path):
    with open(path) as fh:
        return json.load(fh)


def run(argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = am.main(argv)
    return rc, out.getvalue()


def sam_line(flag=0, mapq=60, cigar="100M", tags=("NH:i:1", "NM:i:0")):
    return ["um1", str(flag), "ctg1", "1001", str(mapq), cigar, "*", "0", "0", "A" * 100, "*"] + list(tags)


FLT = {"min_mapq": 60, "max_nh": 1, "reject_zs": True, "max_nm": 5}


class AlignmentFilterTest(unittest.TestCase):
    def test_unique_kept(self):
        c, h = am.classify_alignment(sam_line(), FLT)
        self.assertEqual(c, "kept")
        self.assertEqual((h["contig"], h["pos"], h["strand"]), ("ctg1", 1001, "+"))

    def test_reverse_strand(self):
        self.assertEqual(am.classify_alignment(sam_line(flag=16), FLT)[1]["strand"], "-")

    def test_rejections(self):
        cases = {"unmapped": sam_line(flag=4), "multi_nh": sam_line(tags=("NH:i:2", "NM:i:0")),
                 "low_mapq": sam_line(mapq=1), "second_best": sam_line(tags=("NH:i:1", "NM:i:0", "ZS:i:-3")),
                 "divergent": sam_line(tags=("NH:i:1", "NM:i:6")), "spliced": sam_line(cigar="50M200N50M"),
                 "clipped": sam_line(cigar="5S95M")}
        for want, f in cases.items():
            self.assertEqual(am.classify_alignment(f, FLT)[0], want, want)

    def test_xa_sa_alone_do_not_make_unique(self):
        # the reference script's grep -v XA/SA would keep this multi-mapper (HISAT2 never writes XA/SA)
        self.assertEqual(am.classify_alignment(sam_line(mapq=1, tags=("NH:i:5", "NM:i:0")), FLT)[0], "multi_nh")


class ReferencePositionsTest(unittest.TestCase):
    def test_remap_keeps_only_unique_same_lg(self):
        d = tmpdir()
        try:
            os.makedirs(os.path.join(d, "out", "x", "samples", "REF", "map"))
            with open(os.path.join(d, "out", "x", "samples", "REF", "map", "hits.tsv"), "w") as fh:
                fh.write("contig\tpos\tunique_id\tstrand\tnm\tmapq\tlg\tmap_pos\n")
                fh.write("R1\t2500000\tu1\t+\t0\t60\t1\t9.9\nR2\t10\tu2\t+\t0\t60\t1\t1.0\n")

            class C:
                outdir = os.path.join(d, "out", "x")

            class S:
                params = {"position_source": {"sample": "REF", "lg_to_seq": {"1": "R1"}}}
            hits = [{"contig": "c", "pos": 5, "uid": u, "lg": "1", "map_pos": 0.0} for u in ("u1", "u2", "u3")]
            out, rec = am.reference_positions(C, S, hits)
            self.assertEqual([h["uid"] for h in out], ["u1"])
            self.assertEqual(out[0]["map_pos"], 2.5)
            self.assertEqual(rec["dropped"], {"reference_sequence_not_its_lg": 1, "not_unique_in_reference": 1})
        finally:
            shutil.rmtree(d)


class SupportTest(unittest.TestCase):
    P = {"min_contig_bp": 1000, "min_unique_markers": 4, "min_distinct_positions": 3,
         "max_markers_per_contig": 3, "conflict_min_markers": 3, "conflict_min_span_bp": 100,
         "mixed_lg_min_frac": 0.1, "mixed_lg_min_markers": 3, "orientation_min_abs_rho": 0.5,
         "orientation_min_positions": 3}

    def hits(self, ctg, lgs_pos):
        return [{"contig": ctg, "pos": p, "uid": f"u{ctg}{i}", "lg": lg, "map_pos": mp}
                for i, (p, lg, mp) in enumerate(lgs_pos)]

    def test_overlap_flags(self):
        rows = [{"contig": "big", "object": "chr01", "order": 1, "length": 1000, "q05": 0.0, "q95": 10.0},
                {"contig": "copy", "object": "chr01", "order": 2, "length": 200, "q05": 2.0, "q95": 6.0},
                {"contig": "next", "object": "chr01", "order": 3, "length": 500, "q05": 11.0, "q95": 20.0},
                {"contig": "other", "object": "chr02", "order": 1, "length": 100, "q05": 2.0, "q95": 6.0}]
        fl = am.overlap_flags(rows)
        self.assertEqual(fl["copy"], {"OVERLAPS_LARGER", "ORDER_AMBIGUOUS"})
        self.assertEqual(fl["big"], {"ORDER_AMBIGUOUS"})
        self.assertNotIn("next", fl)
        self.assertNotIn("other", fl)   # another pseudomolecule
        self.assertEqual(am.quantile([1, 2, 3, 4, 5], 0.5), 3)

    def test_even_sample_deterministic_and_spread(self):
        items = list(range(100))
        s = am.even_sample(items, 4)
        self.assertEqual(s, am.even_sample(items, 4))
        self.assertEqual(s, [12, 37, 62, 87])
        self.assertEqual(am.even_sample(items[:3], 4), [0, 1, 2])

    def test_spearman(self):
        self.assertAlmostEqual(am.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(am.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertIsNone(am.spearman([1, 2, 3], [5, 5, 5]))

    def test_conflict_block_and_majority(self):
        h = self.hits("c1", [(100 * i, "1", i) for i in range(1, 9)] + [(1000 + 100 * i, "2", 50 + i) for i in range(4)])
        rows, sel = am.contig_support({"c1": 5000}, h, self.P)
        r = rows["c1"]
        self.assertEqual(r["majority_lg"], "1")
        self.assertEqual(r["status"], "ELIGIBLE")
        self.assertIn("CONFLICT_BLOCK", r["flags"])
        self.assertIn("MIXED_LG", r["flags"])
        self.assertEqual(len(r["conflict_blocks"]), 1)
        self.assertEqual(r["conflict_blocks"][0]["lg"], "2")
        self.assertTrue(all(m[1] == "1" for m in sel["c1"]))   # ALLMAPS input: majority LG only
        self.assertEqual(len(sel["c1"]), 3)                     # capped, evenly spaced

    def test_thresholds_use_unique_positions(self):
        # 6 unique markers but only 2 distinct map positions: too few informative markers
        h = self.hits("c2", [(100 * i, "1", 5.0 if i < 3 else 6.0) for i in range(6)])
        rows, sel = am.contig_support({"c2": 5000, "c3": 500}, h, self.P)
        self.assertEqual(rows["c2"]["status"], "TOO_FEW_MARKERS")
        self.assertEqual(rows["c3"]["status"], "NO_MARKERS")
        self.assertNotIn("c2", sel)

    def test_too_short(self):
        h = self.hits("c4", [(10 * i, "1", i) for i in range(10)])
        rows, _ = am.contig_support({"c4": 500}, h, self.P)
        self.assertEqual(rows["c4"]["status"], "TOO_SHORT")

    def test_map_csv_header_and_rows(self):
        d = tmpdir()
        try:
            p = os.path.join(d, "m.csv")
            n = am.write_map_csv(p, {"c1": [(10, "1", 0.5, "u1"), (20, "1", 0.7, "u2")]}, {"1": "01"})
            with open(p) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(n, 2)
            self.assertEqual(lines[0], "Scaffold ID,scaffold position,LG,genetic position")
            self.assertEqual(lines[1], "c1,10,01,0.5")
        finally:
            shutil.rmtree(d)


class AgpVerifyTest(unittest.TestCase):
    def setUp(self):
        self.d = tmpdir()
        self.src = os.path.join(self.d, "src.fa")
        allmaps_data.write_fasta(self.src, [("a", "ACGTACGTAA"), ("b", "GGGCCCTTTA"), ("u", "TTTT")])

    def tearDown(self):
        for f in getattr(self, "_fais", []):
            f.close()
        shutil.rmtree(self.d)

    def fai(self):
        f = am.Fai(self.src)
        self._fais = getattr(self, "_fais", []) + [f]
        return f

    def write(self, agp_rows, out_recs):
        agp = os.path.join(self.d, "o.agp")
        with open(agp, "w") as fh:
            fh.write("##agp-version 2.1\n")
            for r in agp_rows:
                fh.write("\t".join(map(str, r)) + "\n")
        out = os.path.join(self.d, "o.fasta")
        allmaps_data.write_fasta(out, out_recs)
        return agp, out

    GOOD = [("chr01", 1, 10, 1, "W", "a", 1, 10, "+"), ("chr01", 11, 13, 2, "U", 3, "scaffold", "yes", "map"),
            ("chr01", 14, 23, 3, "W", "b", 1, 10, "-"), ("u", 1, 4, 1, "W", "u", 1, 4, "?")]

    def test_good_reconstruction(self):
        agp, out = self.write(self.GOOD, [("chr01", "ACGTACGTAA" + "NNN" + "TAAAGGGCCC"), ("u", "TTTT")])
        v, _ = am.verify_agp(agp, out, self.fai())
        self.assertTrue(v["ok"], v)
        self.assertEqual((v["gaps"], v["gap_bp"], v["output_bp"], v["source_bp"]), (1, 3, 27, 24))

    def test_altered_base_detected(self):
        agp, out = self.write(self.GOOD, [("chr01", "ACGTACGTAA" + "NNN" + "TAAAGGGCCA"), ("u", "TTTT")])
        v, _ = am.verify_agp(agp, out, self.fai())
        self.assertFalse(v["ok"])
        self.assertEqual(v["n_records_mismatch"], 1)

    def test_wrong_orientation_detected(self):
        rows = list(self.GOOD)
        rows[2] = ("chr01", 14, 23, 3, "W", "b", 1, 10, "+")
        agp, out = self.write(rows, [("chr01", "ACGTACGTAA" + "NNN" + "TAAAGGGCCC"), ("u", "TTTT")])
        self.assertFalse(am.verify_agp(agp, out, self.fai())[0]["ok"])

    def test_missing_and_duplicated_contig_detected(self):
        rows = self.GOOD[:3] + [("u2", 1, 10, 1, "W", "a", 1, 10, "+")]
        agp, out = self.write(rows, [("chr01", "ACGTACGTAA" + "NNN" + "TAAAGGGCCC"), ("u2", "ACGTACGTAA")])
        v, _ = am.verify_agp(agp, out, self.fai())
        self.assertFalse(v["ok"])
        self.assertEqual((v["n_missing"], v["n_used_twice"]), (1, 1))

    def test_partial_component_detected(self):
        rows = [("chr01", 1, 5, 1, "W", "a", 1, 5, "+")] + self.GOOD[1:]
        agp, out = self.write(rows, [("chr01", "ACGTA"), ("u", "TTTT")])
        self.assertFalse(am.verify_agp(agp, out, self.fai())[0]["ok"])


class MarkersStageTest(unittest.TestCase):
    """The markers stage runs without containers."""

    def test_dedup_and_traceability(self):
        d = tmpdir()
        try:
            data = allmaps_data.make(d)
            cfgp = write_config(d, data)
            rc, _ = run(["-c", cfgp, "run", "--stages", "markers"])
            self.assertEqual(rc, 0)
            md = os.path.join(d, "out", "amtest", "markers")
            v = load(os.path.join(md, "VALIDATED.json"))["validation"]
            self.assertEqual(v["marker_records"], data["n_records"])
            self.assertEqual(v["distinct_sequences"], data["n_unique"])
            self.assertEqual(v["redundant_records"], data["n_records"] - data["n_unique"])
            self.assertEqual(v["inconsistent_duplicate_sequences"], 0)
            with gzip.open(os.path.join(md, "marker_to_unique.tsv.gz"), "rt") as fh:
                rows = [l.rstrip("\n").split("\t") for l in fh][1:]
            self.assertEqual(len(rows), data["n_records"])
            self.assertEqual(len({r[1] for r in rows}), data["n_unique"])
            # rerun is a no-op (validated), never an overwrite
            rc, _ = run(["-c", cfgp, "run", "--stages", "markers"])
            self.assertEqual(rc, 0)
        finally:
            shutil.rmtree(d)

    def test_inconsistent_duplicate_fails(self):
        d = tmpdir()
        try:
            data = allmaps_data.make(d)
            with gzip.open(data["markers_txt"], "rt") as fh:
                lines = fh.read().splitlines()
            # pg_7r1 has the same sequence as pg_7: give it another map position
            lines = [l.replace("\tpg_7r1\t", "\tpg_7r1\t9.9") if "\tpg_7r1\t" in l else l for l in lines]
            lines = [("\t".join(l.split("\t")[:2] + ["9.9"]) if "\tpg_7r1" in l else l) for l in lines]
            with gzip.open(data["markers_txt"], "wt") as fh:
                fh.write("\n".join(lines) + "\n")
            cfgp = write_config(d, data)
            rc, _ = run(["-c", cfgp, "run", "--stages", "markers"])
            self.assertEqual(rc, 1)
            self.assertTrue(os.path.exists(os.path.join(d, "out", "amtest", "markers.partial", "FAILED.json")))
        finally:
            shutil.rmtree(d)


@unittest.skipUnless(HISAT2 and JCVI, "set ALLMAPS_HISAT2_IMAGE and ALLMAPS_JCVI_IMAGE to run the container tests")
class EndToEndTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tmpdir()
        cls.data = allmaps_data.make(cls.d)
        cls.cfg = write_config(cls.d, cls.data)
        cls.rc, cls.out = run(["-c", cls.cfg, "run"])
        cls.base = os.path.join(cls.d, "out", "amtest", "samples", "SAMPLE_A")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def placement(self, variant="scaffold"):
        with open(os.path.join(self.base, variant, "placement.tsv")) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            return {r[0]: dict(zip(head, r)) for r in (l.rstrip("\n").split("\t") for l in fh)}

    def test_all_stages_validated(self):
        self.assertEqual(self.rc, 0, self.out)

    def test_filter_rejects_repeat_and_divergent(self):
        v = load(os.path.join(self.base, "map", "filter_stats.json"))
        self.assertFalse(v["xa_sa_tags_present"])
        self.assertGreaterEqual(v["primary_classes"].get("multi_nh", 0) + v["primary_classes"].get("low_mapq", 0), 10)
        self.assertEqual(v["primary_classes"].get("divergent", 0), 5)

    def check_truth(self, variant, inverted=()):
        pl = self.placement(variant)
        truth = self.data["truth"]
        for lg in ("1", "2", "3"):
            ctgs = sorted((int(pl[c]["order"]), c) for c in pl if pl[c]["object"] == f"chr0{lg}")
            names = [c for _, c in ctgs]
            want = [c for c in sorted(truth, key=lambda c: truth[c][1]) if truth[c][0] == lg]
            if lg == "3":
                want = ["ctg_chim"] + want  # placed by its LG3 majority (60 of 100 kb), before ctg07
            ok = names == want or names == want[::-1]
            self.assertTrue(ok, (lg, names, want))
            flip = names != want
            for c in names:
                if c == "ctg_chim":
                    continue
                o = truth[c][3]
                exp = o if not flip else {"+": "-", "-": "+"}[o]
                if c in inverted:  # the map itself is inverted here: ALLMAPS follows the map
                    exp = {"+": "-", "-": "+"}[exp]
                self.assertEqual(pl[c]["orientation"], exp, (variant, c))

    def test_supplied_positions_follow_the_map(self):
        """Order everywhere and orientation outside the inverted map segment match the truth; ctg05, which lies
        entirely in the segment, is oriented as the (inverted) map says."""
        self.check_truth("scaffold", inverted=("ctg05",))

    def test_reference_positions_recover_truth(self):
        self.check_truth("scaffold_refpos")
        v = load(os.path.join(self.base, "scaffold_refpos", "scaffold_stats.json"))
        self.assertEqual(v["position_source"]["reference_sample"], "REF")
        self.assertGreater(v["position_source"]["hits_kept"], 400)

    def test_no_false_overlap_flags(self):
        pl = self.placement("scaffold_refpos")
        self.assertFalse([c for c, r in pl.items() if "OVERLAPS_LARGER" in r["flags"] or "ORDER_AMBIGUOUS" in r["flags"]])

    def test_chimera_flagged_not_split(self):
        pl = self.placement()
        self.assertIn("CONFLICT_BLOCK", pl["ctg_chim"]["flags"])
        self.assertIn("MIXED_LG", pl["ctg_chim"]["flags"])
        self.assertEqual(pl["ctg_chim"]["object"], "chr03")
        self.assertEqual(pl["ctg_chim"]["length"], "100000")

    def test_unplaced_reasons(self):
        pl = self.placement()
        self.assertEqual(pl["ctg_short"]["status"], "UNPLACED_TOO_SHORT")
        self.assertEqual(pl["ctg_nomark"]["status"], "UNPLACED_NO_MARKERS")

    def test_reconstruction_verified(self):
        v = load(os.path.join(self.base, "scaffold", "agp_verification.json"))
        self.assertTrue(v["ok"], v)
        self.assertEqual(v["n_missing"], 0)
        self.assertEqual(v["n_used_twice"], 0)
        self.assertEqual(v["output_bp"], v["source_bp"] + v["gap_bp"])

    def test_deterministic_rerun(self):
        """A second, independent run (new output profile) gives a byte-identical AGP."""
        with open(self.cfg) as fh:
            c = yaml.safe_load(fh)
        c["profile"] = "amtest_rerun"
        p2 = os.path.join(self.d, "config2.yaml")
        with open(p2, "w") as fh:
            yaml.safe_dump(c, fh)
        rc, out = run(["-c", p2, "run"])
        self.assertEqual(rc, 0, out)
        a = os.path.join(self.base, "scaffold", "SAMPLE-A.v1.allmaps.agp")
        b = os.path.join(self.d, "out", "amtest_rerun", "samples", "SAMPLE_A", "scaffold", "SAMPLE-A.v1.allmaps.agp")
        def strip(p):
            with open(p) as fh:
                return [l for l in fh if not l.startswith("#")]
        self.assertEqual(strip(a), strip(b))


if __name__ == "__main__":
    unittest.main()
