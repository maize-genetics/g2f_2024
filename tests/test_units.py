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
        self.assertEqual(list(cfg["science"]["qc"]["busco"]["lineages"]), ["poales_odb12.2"])
        self.assertEqual(cfg["science"]["qc"]["busco"]["lineages"]["poales_odb12.2"]
                         ["datasets_version"], "odb12.2")
        self.assertFalse(cfg["science"]["kmer"]["screen"]["block_assembly_on_flag"])
        ext = cfg["external_assemblies"]["SAMPLE_A"]
        self.assertEqual(set(ext), {"previous_hifi_l0"})
        self.assertTrue(ext["previous_hifi_l0"]["p_ctg_fa"].endswith(
            "previous_assemblies/SAMPLE_A/SAMPLE_A.asm.bp.p_ctg.fa"))
        self.assertIn("prefiltered_ont",
                      cfg["samples"]["SAMPLE_A"]["acknowledge"]["SAMPLE_A_ONT_gt50kb.fastq.gz"])
        self.assertEqual(cfg["samples"]["SAMPLE_A"]["aliases"], ["SampleA_old_name"])

    def test_missing_samples_file_is_optional(self):
        cfg = load_config(dict(BASE_TOP, samples=f"{self.tmp}/none.yaml", data_dir=f"{self.tmp}/d",
                               results_dir=f"{self.tmp}/r", local=f"{self.tmp}/none.yaml"), PIPELINE_DIR)
        self.assertEqual(cfg["samples"], {})

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
        movie = "m84000_240101_000001_s1"
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
                self.assertEqual(fh.read(), "uuid-1\t60000\nuuid-2\t70000\n")
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

    def ids_file(self, name, rows):
        path = f"{self.tmp}/{name}.ids.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("".join(f"{r}\t{n}\n" for r, n in rows))
        return path

    def ont_case(self, filt_n=(3, 2), ids=(["r1", "r2"], ["r3"]), prefiltered=()):
        with open(self.inv) as fh:
            inv = json.load(fh)
        for i in inv["inputs"]:
            if i["input_id"] in prefiltered:
                i["notes"] = ["ack:prefiltered_ont"]
        inv_path = jdump(f"{self.tmp}/inv2.json", inv)
        sup, ret, sup_ids, ret_ids = {}, {}, {}, {}
        for iid, fn, rid in zip(("o1", "o2"), filt_n, ids):
            sup[iid] = jdump(f"{self.tmp}/{iid}.in.json", {"n_reads": 10, "total_bases": 300000,
                                                           "duplicate_read_ids": 0, "min_len": 1000})
            ret[iid] = jdump(f"{self.tmp}/{iid}.ret.json", {"n_reads": fn, "total_bases": fn * 60000,
                                                            "min_len": 60000 if fn else 0})
            sup_ids[iid] = self.ids_file(f"{iid}.in", [(r, 30000) for r in rid])
            ret_ids[iid] = self.ids_file(f"{iid}.ret", [(f"{iid}_{k}", 60000) for k in range(fn)])
        return validate.validate_ont(sup, ret, sup_ids, ret_ids, 50000, 1, inv_path, 1000000)

    def test_ont_ok_reports_raw_and_retained(self):
        r = self.ont_case()
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual((r["raw_reads"], r["retained_reads"]), (20, 5))
        self.assertEqual(r["input_type"], "raw")
        self.assertEqual(r["retained_approx_coverage"], 0.3)
        self.assertEqual(r["retained_read_n50"], 60000)

    def test_ont_prefiltered_yield_is_not_reported_as_raw(self):
        r = self.ont_case(prefiltered=("o1", "o2"))
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual(r["input_type"], "prefiltered")
        self.assertIsNone(r["raw_reads"])
        self.assertIsNone(r["raw_bases"])
        self.assertEqual(r["supplied_reads"], 20)
        self.assertEqual(r["per_input"]["o1"]["yield_type"], "prefiltered_supplied")
        self.assertTrue(r["warnings"])
        mixed = self.ont_case(prefiltered=("o1",))
        self.assertEqual(mixed["input_type"], "mixed")
        self.assertIsNone(mixed["raw_reads"])

    def test_prefiltered_file_with_short_reads_fails(self):
        with open(self.inv) as fh:
            inv = json.load(fh)
        inv["inputs"][2]["notes"] = ["ack:prefiltered_ont"]
        inv_path = jdump(f"{self.tmp}/inv3.json", inv)
        st = jdump(f"{self.tmp}/pf.json", {"n_reads": 2, "total_bases": 100000, "min_len": 49999,
                                           "duplicate_read_ids": 0})
        ids = self.ids_file("pf", [("a", 49999), ("b", 50001)])
        r = validate.validate_ont({"o1": st}, {"o1": st}, {"o1": ids}, {"o1": ids}, 50000, 1,
                                  inv_path, 1000000)
        self.assertFalse(r["ok"])
        self.assertTrue(any("shorter than 50000" in e for e in r["errors"]))

    def test_ont_input_without_qualifying_reads_needs_review(self):
        r = self.ont_case(filt_n=(3, 0))
        self.assertFalse(r["ok"])
        self.assertTrue(any("human review" in e for e in r["errors"]))

    def test_ont_duplicate_ids_between_inputs(self):
        r = self.ont_case(ids=(["r1", "r2"], ["r2"]))
        self.assertFalse(r["ok"])
        self.assertTrue(any("occurs in both" in e for e in r["errors"]))

    def test_stats_from_other_aggregator_code_rejected(self):
        st = self.hstats("a", ["m1"], aggregator_id="old")
        r = validate.validate_hifi({"a": st}, {}, self.inv, 1000, aggregator="new")
        self.assertFalse(r["ok"])
        self.assertTrue(any("different read-statistics code" in e for e in r["errors"]))

    def test_busco_fixture(self):
        r = validate.validate_busco(os.path.join(FIX, "busco_short_summary_synthetic_liliopsida_odb10.json"),
                                    None)
        self.assertTrue(r["ok"], r["errors"])
        self.assertEqual(r["lineage"], "liliopsida_odb10")
        self.assertEqual((r["complete_pct"], r["single_pct"], r["duplicated_pct"]), (95.8, 80.3, 15.5))
        self.assertEqual(r["n_markers"], 3236)
        self.assertEqual(r["internal_stop_codon_count"], 96)
        self.assertEqual(r["gene_predictor"], "miniprot")

    def test_busco_pinned_dataset_enforced(self):
        fx = os.path.join(FIX, "busco_short_summary_synthetic_liliopsida_odb10.json")
        cfg = write_bytes(f"{self.tmp}/dataset.cfg", b"name=liliopsida_odb10\ncreation_date=2024-01-08\n"
                                                    b"number_of_BUSCOs=3236\nOrthoDB_version=10.1\n")
        ok = validate.validate_busco(fx, cfg, "6.1.0", {"name": "liliopsida_odb10",
                                                         "creation_date": "2024-01-08"})
        self.assertTrue(ok["ok"], ok["errors"])
        bad_ver = validate.validate_busco(fx, cfg, "6.0.0", None)
        self.assertFalse(bad_ver["ok"])
        wrong = validate.validate_busco(fx, cfg, "6.1.0", {"name": "poales_odb12.2",
                                                            "datasets_version": "odb12.2"})
        self.assertFalse(wrong["ok"])
        self.assertTrue(any("odb12.2" in e for e in wrong["errors"]))
        missing = validate.validate_busco(fx, f"{self.tmp}/nope.cfg", None, None)
        self.assertFalse(missing["ok"])

    def test_empty_alternate_contig_set_is_valid(self):
        st = write_bytes(f"{self.tmp}/a.tsv", b"file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len"
                         b"\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)"
                         b"\tsum_n\na.fa\t\t\t0\t0\t0\t0.0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0.00\t0.00\t0\n")
        chk = jdump(f"{self.tmp}/c.json", {"ok": True, "n_fasta_records": 0, "fasta_total_bp": 0})
        r = validate.validate_assembly_stats(st, chk)
        self.assertTrue(r["ok"], r["errors"])
        self.assertTrue(r["empty"])

    def test_busco_inconsistent_counts_fail(self):
        with open(os.path.join(FIX, "busco_short_summary_synthetic_liliopsida_odb10.json")) as fh:
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

    def test_seqkit_stats_fixture(self):
        st = validate.parse_seqkit_stats(os.path.join(FIX, "seqkit_stats_synthetic.txt"))
        self.assertEqual(st["n_contigs"], 312)
        self.assertEqual(st["total_bp"], 1234567890)
        self.assertEqual(st["n50"], 45678901)
        self.assertEqual(st["l50"], 9)
        self.assertEqual(st["gc_pct"], 44.12)
        self.assertEqual(st["gap_bp"], 0)

    def test_hifiasm_log_fixture(self):
        h = validate.parse_hifiasm_log(os.path.join(FIX, "hifiasm_synthetic.stderr.excerpt.log"))
        self.assertEqual(h["version"], "0.25.0-r726")
        self.assertIn("--ul", h["command"])
        self.assertEqual(h["peak_rss_gb"], 98.765)
        self.assertEqual(h["ul_passed_reads"], 9876)
        self.assertEqual(h["peaks"][0], {"stage": "ha_ft_gen", "peak_hom": 25, "peak_het": -1})
        self.assertEqual(h["peaks"][1], {"stage": "ha_pt_gen", "peak_hom": 48, "peak_het": 24})

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

    def test_real_meryl_spectrum_of_homozygous_data_not_flagged(self):
        h = kmer.read_histogram(os.path.join(FIX, "meryl_k21_synthetic_homozygous_300kb.hist"))
        r = kmer.screen(h)
        self.assertEqual(r["main_peak_coverage"], 26)
        self.assertFalse(r["flagged"], r)          # unsmoothed v1 flagged sampling noise
        self.assertTrue(kmer.screen(h, smooth_frac=0, min_local_max_ratio=0)["flagged"])

    def test_short_histogram_is_error(self):
        self.assertFalse(kmer.screen([(1, 5), (2, 3)])["ok"])


class SummaryTest(unittest.TestCase):
    def test_refuses_unvalidated_records(self):
        ok = {"ok": True}
        with self.assertRaises(ValueError):
            summary.build("S", "m", 1, {"ok": False}, None, ok, ok, ok, [], None, {})

    def test_pipeline_summary_requires_read_records(self):
        asm = {"ok": True, "n_contigs": 1, "total_bp": 10, "n50": 10, "l50": 1, "max_len": 10,
               "gc_pct": 50, "gap_bp": 0}
        with self.assertRaises(ValueError):
            summary.build("S", "m", 1, None, None, None, asm, None, [], None, {})
        rows = summary.build("S", "imported-x", 1, None, None, None, asm, None, [], None,
                             {"version": "0.25.0-r726"}, source="imported_legacy", validation_id="v")
        self.assertEqual(rows[0]["assembly_source"], "imported_legacy")
        self.assertEqual(rows[0]["hifi_reads"], "")


class IdentityTest(unittest.TestCase):
    def test_validation_identity_tracks_science_methods_and_code(self):
        from g2f.config import method_fingerprint, validation_identity
        tmp = tempfile.mkdtemp()
        try:
            code = write_bytes(f"{tmp}/v.py", b"x = 1\n")
            cfg = {"pipeline_dir": tmp, "science": {"a": 1},
                   "methods": {"m": {"method_id": "m", "args": ["-l0"], "description": "d"}}}
            v1, _ = validation_identity(cfg, [code])
            cfg2 = dict(cfg, science={"a": 2})
            self.assertNotEqual(v1, validation_identity(cfg2, [code])[0])
            cfg3 = dict(cfg, methods={"m": {"method_id": "m", "args": ["-l1"], "description": "d"}})
            self.assertNotEqual(v1, validation_identity(cfg3, [code])[0])
            write_bytes(code, b"x = 2\n")
            self.assertNotEqual(v1, validation_identity(cfg, [code])[0])
            # description edits and file paths do not change the scientific fingerprint
            m = {"method_id": "m", "args": ["-l0"], "description": "a", "_file": "/x"}
            self.assertEqual(method_fingerprint(m),
                             method_fingerprint(dict(m, description="b", _file="/y")))
        finally:
            shutil.rmtree(tmp)


if __name__ == "__main__":
    unittest.main()


class DiskPlanTest(unittest.TestCase):
    """Start gate and prepared-FASTQ retention decisions of the full-scale launcher."""

    def setUp(self):
        from g2f import diskplan
        self.dp = diskplan
        rec = lambda bam: {"inputs": [{"tech": "hifi", "status": "use", "format": "bam",  # noqa: E731
                                       "size": bam, "notes": []},
                                      {"tech": "ont", "status": "use", "format": "tar",
                                       "size": 10e9, "notes": []}]}
        self.est = {s: diskplan.estimate(rec(b)) for s, b in
                    (("A", 80e9), ("B", 50e9), ("C", 50e9), ("D", 40e9))}

    def test_estimate_components(self):
        e = self.est["A"]
        self.assertAlmostEqual(e["peak"], e["persist_keep"] + e["bins"], delta=2)
        self.assertAlmostEqual(e["persist_keep"] - e["persist_temp"], e["conv"], delta=2)

    def test_start_and_retention(self):
        d = self.dp.decide("A", self.est, 2000e9, {}, ["A", "B", "C", "D"], 150e9, 3)
        self.assertTrue(d["start"] and d["retain_prepared_hifi"])
        d = self.dp.decide("A", self.est, 600e9, {}, ["A", "B", "C", "D"], 150e9, 3)
        self.assertTrue(d["start"])
        self.assertFalse(d["retain_prepared_hifi"])                 # space needed later
        d = self.dp.decide("A", self.est, 600e9, {}, ["A", "B", "C", "D"], 150e9, 3, mode="keep")
        self.assertTrue(d["retain_prepared_hifi"])

    def test_running_genotypes_reserve_their_remaining_growth(self):
        run = {"B": {"used": 10e9, "assembled": False}, "C": {"used": 10e9, "assembled": False}}
        free = self.est["A"]["peak"] + 150e9 + 50e9
        d = self.dp.decide("A", self.est, free, run, ["A", "D"], 150e9, 3)
        self.assertFalse(d["start"])
        run = {k: dict(v, assembled=True, used=self.est[k]["persist_keep"]) for k, v in run.items()}
        d = self.dp.decide("A", self.est, free, run, ["A", "D"], 150e9, 3)
        self.assertTrue(d["start"])

    def test_resumed_genotype_needs_only_what_is_missing(self):
        own = {"used": self.est["A"]["persist_keep"], "assembled": True}
        d = self.dp.decide("A", self.est, 200e9, {}, ["A"], 150e9, 3, own=own)
        self.assertTrue(d["start"])


class VerifyPreparedTest(unittest.TestCase):
    """Run accounting of the pre-hifiasm gate against the raw-data manifest."""

    def _inv(self, *inputs):
        return {"inputs": [dict(tech=t, input_id=i, status=st, path=f"/d/G/{t}/{n}", format=fmt,
                                size=sz, fingerprint="f", notes=[])
                           for t, i, st, n, fmt, sz in inputs]}

    def test_excluded_run_is_accounted_and_counts_checked(self):
        from g2f import prepared
        inv = self._inv(("hifi", "h1", "use", "h1.bam", "bam", 10),
                        ("ont", "o1", "use", "o1_fastq.tar", "tar", 20),
                        ("ont", "o2", "excluded", "o2_fastq.tar", "", 0))
        man = [{"genotype": "G", "technology": t, "source_path": f"/src/{n}", "size_bytes": str(sz)}
               for t, n, sz in (("HiFi", "h1.bam", 10), ("ONT", "o1_fastq.tar", 20),
                                ("ONT", "o2_fastq.tar", 5))]
        chk = lambda p, n: {"path": p, "size": 1, "mtime": 1, "md5": "m", "gzip_ok": True,  # noqa: E731
                            "error": "", "records": n}
        checks = {("hifi", "h1"): chk("/r/h1.fastq.gz", 7), ("ont", "o1"): chk("/r/o1.fastq.gz", 3)}
        hv = {"per_input": {"h1": {"n_reads": 7, "total_bases": 70, "bam_primary_records": 7}}}
        ov = {"per_input": {"o1": {"retained": {"n_reads": 3, "total_bases": 30},
                                   "supplied": {"n_reads": 9}}}}
        res = prepared.verify_prepared("G", inv, hv, ov, checks, {}, ["/r/h1.fastq.gz", "/r/o1.fastq.gz"],
                                       min_len=50000, manifest_rows=man)
        self.assertTrue(res["ok"], res["errors"])
        bad = dict(checks)
        bad[("ont", "o1")] = chk("/r/o1.fastq.gz", 2)
        res = prepared.verify_prepared("G", inv, hv, ov, bad, {}, ["/r/h1.fastq.gz", "/r/o1.fastq.gz"],
                                       min_len=50000, manifest_rows=man)
        self.assertFalse(res["ok"])
        res = prepared.verify_prepared("G", inv, hv, ov, checks, {}, ["/r/h1.fastq.gz", "/r/o1.fastq.gz"],
                                       min_len=50000, manifest_rows=man[:2])   # o2 not in manifest
        self.assertIn("unexpected run", " ".join(res["errors"]))


class CpuSlotsTest(unittest.TestCase):
    def test_disjoint_blocks_and_refusal(self):
        from g2f import cpus
        self.assertEqual(cpus.compress([0, 1, 2, 5, 7, 8]), "0-2,5,7-8")
        self.assertEqual(cpus.slots(3, 2, allowed=range(8)), ["0-1", "2-3", "4-5"])
        self.assertEqual(cpus.slots(2, 2, allowed=[0, 1, 4, 5, 6]), ["0-1", "4-5"])
        with self.assertRaises(ValueError):
            cpus.slots(3, 40, allowed=range(119))                      # 120 > 119


class GenotypeWrapperPinTest(unittest.TestCase):
    """scripts/run_one_genotype.sh pins Snakemake (and so every tool) to its slot's CPU
    block; without a slot (an older launcher) it runs unpinned as before."""
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def run_wrapper(self, slot_args, cpu_slots):
        import subprocess
        import sys
        t = self.tmp
        os.makedirs(f"{t}/sandbox/scripts", exist_ok=True)
        os.makedirs(f"{t}/bin", exist_ok=True)
        os.makedirs(f"{t}/logs/state", exist_ok=True)
        shutil.copy(f"{PIPELINE_DIR}/scripts/run_one_genotype.sh", f"{t}/sandbox/scripts/")
        with open(f"{t}/sandbox/scripts/g2f_tool.py", "w") as fh:          # stub disk gate
            fh.write("import json, sys\n"
                     "if sys.argv[1] == 'disk-gate':\n"
                     f"    open('{t}/logs/state/G.running', 'w').close()\n"
                     "    print(json.dumps({'sample': 'G', 'start': True, 'retain_prepared_hifi': True}))\n")
        with open(f"{t}/bin/snakemake", "w") as fh:                          # stub: record CPUs
            fh.write(f"#!{sys.executable}\nimport os\n"
                     f"open('{t}/affinity', 'w').write(' '.join(map(str, sorted(os.sched_getaffinity(0)))))\n")
        os.chmod(f"{t}/bin/snakemake", 0o755)
        env = dict(os.environ, PATH=f"{t}/bin:" + os.environ["PATH"], LOGDIR=f"{t}/logs",
                   STATE=f"{t}/logs/state", RESDIR=f"{t}/results", RESOURCES="r.yaml", THREADS="1",
                   MEM_GB="1", RESERVE_GB="0", JOBS="2", RETAIN_MODE="keep", MAX_WAIT_H="0")
        if cpu_slots:
            env["CPU_SLOTS"] = cpu_slots
        r = subprocess.run(["bash", f"{t}/sandbox/scripts/run_one_genotype.sh", "G", *slot_args],
                           env=env, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + open(f"{t}/logs/launcher.log").read())
        with open(f"{t}/affinity") as fh:
            return set(map(int, fh.read().split()))

    def test_pinned_to_slot_block(self):
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < 2:
            self.skipTest("needs 2 CPUs")
        a, b = allowed[0], allowed[1]
        self.assertEqual(self.run_wrapper(["2"], f"{a} {b}"), {b})
        self.assertEqual(self.run_wrapper([], f"{a} {b}"), set(allowed))      # no slot: unpinned
