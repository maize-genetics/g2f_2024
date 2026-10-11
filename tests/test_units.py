"""Config parsing, read statistics, validators, k-mer screen, summary."""
import gzip
import json
import math
import os
import shutil
import tempfile
import unittest

from util import BASE_TOP, PIPELINE_DIR, hifi_names, write_bytes

from g2f import gfa, kmer, readstats, summary, validate
from g2f.config import ConfigError, load_config

FIX = os.path.join(PIPELINE_DIR, "tests", "fixtures")


def jdump(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh)
    return path


class ConfigTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_tracked_config_loads(self):
        cfg = load_config(dict(BASE_TOP, data_dir=f"{self.tmp}/d", results_dir=f"{self.tmp}/r",
                               local=f"{self.tmp}/none.yaml"), PIPELINE_DIR)
        self.assertEqual(set(cfg["methods"]), {"hifiasm-0.25.0_hifi_l0_v1",
                                               "hifiasm-0.25.0_hifi-ont-ul50k_l0_v1"})
        ont = cfg["methods"]["hifiasm-0.25.0_hifi-ont-ul50k_l0_v1"]
        self.assertEqual(ont["ont"]["min_read_len"], 50000)
        self.assertEqual(ont["args"], ["-l0"])
        self.assertEqual(ont["expected_version"], "0.25.0-r726")
        self.assertEqual(cfg["science"]["qc"]["busco"]["lineages"], ["poales_odb12"])
        self.assertIn("PHBW8", cfg["samples"])
        self.assertEqual(cfg["samples"]["PHBW8"]["aliases"], ["EB2_PHBW8"])

    def test_results_inside_legacy_rejected(self):
        top = dict(BASE_TOP, data_dir=f"{self.tmp}/data", results_dir=f"{self.tmp}/legacy/results",
                   legacy_dirs=[f"{self.tmp}/legacy"], local=f"{self.tmp}/none.yaml")
        with self.assertRaises(ConfigError):
            load_config(top, PIPELINE_DIR)

    def test_bad_method_policy_rejected(self):
        sci = os.path.join(self.tmp, "science.yaml")
        with open(os.path.join(PIPELINE_DIR, "config", "science.yaml")) as fh:
            text = fh.read().replace("ont_absent: hifiasm-0.25.0_hifi_l0_v1",
                                     "ont_absent: hifiasm-0.25.0_hifi-ont-ul50k_l0_v1")
        with open(sci, "w") as fh:
            fh.write(text)
        top = dict(BASE_TOP, science=sci, data_dir=f"{self.tmp}/d", results_dir=f"{self.tmp}/r",
                   local=f"{self.tmp}/none.yaml")
        with self.assertRaisesRegex(ConfigError, "HiFi-only"):
            load_config(top, PIPELINE_DIR)

    def test_override_requires_reason(self):
        smp = os.path.join(self.tmp, "samples.yaml")
        with open(smp, "w") as fh:
            fh.write("samples:\n  S1:\n    method_override: hifiasm-0.25.0_hifi_l0_v1\n")
        top = dict(BASE_TOP, samples=smp, data_dir=f"{self.tmp}/d", results_dir=f"{self.tmp}/r",
                   local=f"{self.tmp}/none.yaml")
        with self.assertRaisesRegex(ConfigError, "override_reason"):
            load_config(top, PIPELINE_DIR)


class ReadStatsTest(unittest.TestCase):
    def test_hifi_stats_census_and_n50(self):
        movie = "m84109_240806_024441_s4"
        lines = [f"{n}\t{l}\t{q}\n" for n, l, q in
                 zip(hifi_names(movie, 4), [10, 20, 30, 40], [30, 30, 20, 10])]
        s = readstats.aggregate(lines, "hifi", True, genome_size=10)
        self.assertEqual(s["n_reads"], 4)
        self.assertEqual(s["total_bases"], 100)
        self.assertEqual(s["read_n50"], 30)      # 40+30 >= 50
        self.assertEqual(s["read_l50"], 2)
        self.assertEqual(s["census"], {movie: 4})
        self.assertEqual(s["approx_coverage"], 10.0)
        self.assertEqual(s["frac_q20"], 0.75)
        self.assertEqual(s["duplicate_read_ids"], 0)

    def test_hifi_duplicate_zmw_and_unparsed(self):
        lines = ["m1_000000_000000/5/ccs\t100\n", "m1_000000_000000/5/ccs\t100\n", "weird\t50\n"]
        s = readstats.aggregate(lines, "hifi", False)
        self.assertEqual(s["duplicate_read_ids"], 1)
        self.assertEqual(s["unparsed_names"], 1)

    def test_ont_names_with_tabs_and_ids(self):
        tmp = tempfile.mkdtemp()
        try:
            ids = os.path.join(tmp, "ids.gz")
            lines = ["uuid-1\tqs:f:12\tRG:Z:runA_model\t60000\t12.5\n",
                     "uuid-2 runid=runB ch=1\t70000\t14.0\n"]
            s = readstats.aggregate(lines, "ont", True, genome_size=1000, ids_out=ids)
            self.assertEqual(s["n_reads"], 2)
            self.assertEqual(s["census"], {"runA": 1, "runB": 1})
            with gzip.open(ids, "rt") as fh:
                self.assertEqual(fh.read().split(), ["uuid-1", "uuid-2"])
        finally:
            shutil.rmtree(tmp)


class ValidateTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_test_")
        self.inv = jdump(f"{self.tmp}/inv.json", {"inputs": [
            {"tech": "hifi", "input_id": "a", "notes": []},
            {"tech": "hifi", "input_id": "b", "notes": []},
            {"tech": "ont", "input_id": "o1", "notes": []},
            {"tech": "ont", "input_id": "o2", "notes": []}]})

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def hstats(self, name, movies, n=10, **kw):
        d = {"n_reads": n, "total_bases": n * 1000, "read_n50": 1000, "census": {m: n for m in movies},
             "name_kinds": {"ccs": n}, "unparsed_names": 0, "duplicate_read_ids": 0}
        d.update(kw)
        return jdump(f"{self.tmp}/{name}.json", d)

    def test_hifi_movie_overlap_between_inputs_fails(self):
        r = validate.validate_hifi({"a": self.hstats("a", ["m1"]), "b": self.hstats("b", ["m1", "m2"])},
                                   {}, self.inv, 1000)
        self.assertFalse(r["ok"])
        self.assertTrue(any("occurs in both" in e for e in r["errors"]))

    def test_hifi_ok_and_bam_count_mismatch(self):
        ok = validate.validate_hifi({"a": self.hstats("a", ["m1"]), "b": self.hstats("b", ["m2"])},
                                    {}, self.inv, 10000)
        self.assertTrue(ok["ok"], ok["errors"])
        self.assertEqual(ok["total_bases"], 20000)
        self.assertEqual(ok["approx_coverage"], 2.0)
        cnt = write_bytes(f"{self.tmp}/a.count", b"11\n")
        bad = validate.validate_hifi({"a": self.hstats("a", ["m1"])}, {"a": cnt}, self.inv, 10000)
        self.assertFalse(bad["ok"])

    def ont_case(self, filt_n=(3, 2), merged_n=5, ids=(["r1", "r2"], ["r3"])):
        raw, filt, idp = {}, {}, {}
        for iid, fn, rid in zip(("o1", "o2"), filt_n, ids):
            raw[iid] = jdump(f"{self.tmp}/{iid}.raw.json", {"n_reads": 10, "total_bases": 300000,
                                                            "duplicate_read_ids": 0})
            filt[iid] = jdump(f"{self.tmp}/{iid}.f.json", {"n_reads": fn, "total_bases": fn * 60000,
                                                           "min_len": 60000 if fn else 0})
            idp[iid] = f"{self.tmp}/{iid}.ids.gz"
            with gzip.open(idp[iid], "wt") as fh:
                fh.write("\n".join(rid) + "\n")
        merged = jdump(f"{self.tmp}/m.json", {"n_reads": merged_n, "total_bases": merged_n * 60000,
                                              "read_n50": 60000})
        return validate.validate_ont(raw, filt, merged, idp, 50000, 1, self.inv, 1000000)

    def test_ont_ok_reports_raw_and_retained(self):
        r = self.ont_case()
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual((r["raw_reads"], r["retained_reads"]), (20, 5))
        self.assertEqual(r["retained_approx_coverage"], 0.3)

    def test_ont_input_without_qualifying_reads_needs_review(self):
        r = self.ont_case(filt_n=(3, 0), merged_n=3)
        self.assertFalse(r["ok"])
        self.assertTrue(any("human review" in e for e in r["errors"]))

    def test_ont_duplicate_ids_and_merge_mismatch(self):
        r = self.ont_case(ids=(["r1", "r2"], ["r2"]), merged_n=4)
        self.assertFalse(r["ok"])
        self.assertTrue(any("occurs in both" in e for e in r["errors"]))
        self.assertTrue(any("merged" in e for e in r["errors"]))

    def test_busco_real_fixture(self):
        r = validate.validate_busco(os.path.join(FIX, "busco_short_summary_PHBW8_l0_liliopsida_odb10.json"),
                                    None)
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual(r["lineage"], "liliopsida_odb10")
        self.assertEqual((r["complete_pct"], r["single_pct"], r["duplicated_pct"]), (99.1, 81.9, 17.3))
        self.assertEqual(r["n_markers"], 3236)
        self.assertEqual(r["internal_stop_codon_count"], 127)
        self.assertEqual(r["gene_predictor"], "miniprot")

    def test_busco_inconsistent_counts_fail(self):
        with open(os.path.join(FIX, "busco_short_summary_PHBW8_l0_liliopsida_odb10.json")) as fh:
            d = json.load(fh)
        d["results"]["Missing BUSCOs"] += 5
        r = validate.validate_busco(jdump(f"{self.tmp}/b.json", d), None)
        self.assertFalse(r["ok"])

    def test_merqury_parse(self):
        qv = write_bytes(f"{self.tmp}/S.qv", b"S.m.p_ctg\t120\t1000000\t45.2\t3.0e-05\n")
        cp = write_bytes(f"{self.tmp}/S.completeness.stats", b"S.m.p_ctg\tall\t990\t1000\t99.0\n")
        r = validate.validate_merqury(qv, cp, "S.m.p_ctg")
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual((r["qv"], r["kmer_completeness_pct"]), (45.2, 99.0))
        bad = write_bytes(f"{self.tmp}/bad.qv", b"garbage\n")
        self.assertFalse(validate.validate_merqury(bad, cp, "S.m.p_ctg")["ok"])

    def test_seqkit_stats_real_fixture(self):
        st = validate.parse_seqkit_stats(os.path.join(FIX, "seqkit_stats_PHBW8_HiFi_Nano_l0.txt"))
        self.assertEqual(st["n_contigs"], 490)
        self.assertEqual(st["total_bp"], 2266382187)
        self.assertEqual(st["n50"], 121136088)
        self.assertEqual(st["l50"], 6)
        self.assertEqual(st["gc_pct"], 46.71)
        self.assertEqual(st["gap_bp"], 0)

    def test_hifiasm_log_real_fixture(self):
        h = validate.parse_hifiasm_log(os.path.join(FIX, "hifiasm_PHBW8_HiFi_Nano_l0.stderr.excerpt.log"))
        self.assertEqual(h["version"], "0.25.0-r726")
        self.assertIn("--ul", h["command"])
        self.assertEqual(h["peak_rss_gb"], 133.736)
        self.assertEqual(h["peaks"][0], {"stage": "ha_ft_gen", "peak_hom": 37, "peak_het": -1})
        self.assertEqual(h["peaks"][1], {"stage": "ha_pt_gen", "peak_hom": 72, "peak_het": 36})

    def test_gfa_fasta_consistency(self):
        g = write_bytes(f"{self.tmp}/a.gfa", b"H\tVN:Z:1.0\nS\tc1\tACGT\tLN:i:4\nS\tc2\tGGNN\tLN:i:4\n")
        f = write_bytes(f"{self.tmp}/a.fa", b">c1\nACGT\n>c2\nGGNN\n")
        self.assertTrue(gfa.compare(g, f)["ok"])
        f2 = write_bytes(f"{self.tmp}/b.fa", b">c1\nACGT\n>c2\nGGN\n")
        self.assertFalse(gfa.compare(g, f2)["ok"])
        f3 = write_bytes(f"{self.tmp}/c.fa", b">c1\nACGT\n")
        self.assertFalse(gfa.compare(g, f3)["ok"])
        g2 = write_bytes(f"{self.tmp}/bad.gfa", b"S\tc1\tACGT\tLN:i:5\n")
        with self.assertRaises(ValueError):
            gfa.compare(g2, f)


class KmerScreenTest(unittest.TestCase):
    @staticmethod
    def hist(peak, het_height=0.0):
        out = []
        for x in range(1, 150):
            err = 5e7 * math.exp(-x / 1.5)
            hom = 3e7 * math.exp(-((x - peak) ** 2) / (2 * (peak / 6) ** 2))
            het = het_height * 3e7 * math.exp(-((x - peak / 2) ** 2) / (2 * (peak / 10) ** 2))
            out.append((x, int(err + hom + het)))
        return out

    def test_inbred_like_spectrum_not_flagged(self):
        r = kmer.screen(self.hist(37))
        self.assertTrue(r["ok"])
        self.assertEqual(r["main_peak_coverage"], 37)
        self.assertFalse(r["flagged"], r)

    def test_heterozygous_spectrum_flagged(self):
        r = kmer.screen(self.hist(60, het_height=0.6))
        self.assertTrue(r["flagged"], r)

    def test_short_histogram_is_error(self):
        self.assertFalse(kmer.screen([(1, 5), (2, 3)])["ok"])


class SummaryTest(unittest.TestCase):
    def test_refuses_unvalidated_records(self):
        ok = {"ok": True}
        with self.assertRaises(ValueError):
            summary.build("S", "m", 1, {"ok": False}, None, ok, ok, ok, [], None, {})


if __name__ == "__main__":
    unittest.main()
