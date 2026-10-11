import os
import shutil
import tempfile
import unittest

from util import (HAVE_SAMTOOLS, hifi_names, make_cfg, mkdirs, ont_names, write_bam,
                  write_bytes, write_fasta, write_fastq)

from g2f import discovery

M1, M2 = "m84109_240806_024441_s4", "m84109_240611_191050_s2"
HIFI_ONLY = "hifiasm-0.25.0_hifi_l0_v1"
HIFI_ONT = "hifiasm-0.25.0_hifi-ont-ul50k_l0_v1"


def codes(issues, severity="error"):
    return sorted({i.code for i in issues if i.severity == severity})


class DiscoveryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def run_discover(self, samples=None, external=None):
        cfg = make_cfg(self.tmp, samples=samples, external=external)
        inputs, issues = discovery.discover(cfg)
        return cfg, inputs, issues

    # -- basic cases ------------------------------------------------------
    def test_multiple_hifi_runs_no_ont_selects_hifi_only(self):
        d = mkdirs(self.tmp, "S1")
        write_fastq(f"{d['hifi']}/run1.hifi.fastq.gz", hifi_names(M1, 20), [12000] * 20)
        write_fasta(f"{d['hifi']}/run2.hifi.fasta.gz", hifi_names(M2, 20), [15000] * 20)
        cfg, inputs, issues = self.run_discover()
        self.assertEqual(codes(issues), [])
        used = [i for i in inputs["S1"] if i.status == "use"]
        self.assertEqual(len(used), 2)
        self.assertEqual({i.format for i in used}, {"fastq", "fasta"})
        method, reason = discovery.select_method("S1", inputs["S1"], cfg)
        self.assertEqual(method, HIFI_ONLY)
        self.assertIn("no ONT", reason)

    def test_multiple_ont_runs_select_hifi_ont(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        write_fastq(f"{d['ont']}/flowcell1.fastq.gz", ont_names(30, seed=1), [5000, 80000] * 15)
        write_fastq(f"{d['ont']}/flowcell2.fastq.gz", ont_names(30, seed=2), [3000, 60000] * 15)
        cfg, inputs, issues = self.run_discover()
        self.assertEqual(codes(issues), [])
        self.assertEqual(len([i for i in inputs["S1"] if i.tech == "ont"]), 2)
        self.assertEqual(discovery.select_method("S1", inputs["S1"], cfg)[0], HIFI_ONT)

    def test_empty_ont_dir_means_absent(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        cfg, inputs, issues = self.run_discover()
        self.assertEqual(codes(issues), [])
        self.assertEqual(discovery.select_method("S1", inputs["S1"], cfg)[0], HIFI_ONLY)

    def test_no_hifi_is_error(self):
        d = mkdirs(self.tmp, "S1", ("ont",))
        write_fastq(f"{d['ont']}/a.fastq.gz", ont_names(5), [60000] * 5)
        _, _, issues = self.run_discover(
            samples={"S1": {"acknowledge": {"a.fastq.gz": {"prefiltered_ont": "test"}}}})
        self.assertIn("NO_HIFI", codes(issues))

    # -- invalid ONT never falls back --------------------------------------
    def test_present_but_invalid_ont_halts(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        with open(f"{d['ont']}/broken.fastq", "w") as fh:
            fh.write("@read1\nACGT\n+\nII\n")      # quality length mismatch
        _, _, issues = self.run_discover()
        self.assertIn("BAD_FORMAT", codes(issues))

    def test_pod5_requires_basecalling(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        write_bytes(f"{d['ont']}/run.pod5", b"\x8bPOD\r\n\x1a\n")
        _, _, issues = self.run_discover()
        self.assertIn("RAW_SIGNAL", codes(issues))

    def test_ont_fasta_rejected(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        write_fasta(f"{d['ont']}/reads.fasta.gz", [n.split()[0] for n in ont_names(5)], [9000] * 5)
        _, _, issues = self.run_discover()
        self.assertIn("UNSUPPORTED_ONT_FASTA", codes(issues))

    def test_prefiltered_ont_needs_acknowledgement(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        write_fastq(f"{d['ont']}/ont_gt50kb.fastq.gz", ont_names(10), [55000] * 10)
        _, _, issues = self.run_discover()
        self.assertIn("ONT_POSSIBLY_PREFILTERED", codes(issues))
        cfg, inputs, issues = self.run_discover(
            samples={"S1": {"acknowledge": {"ont_gt50kb.fastq.gz": {"prefiltered_ont": "legacy"}}}})
        self.assertEqual(codes(issues), [])
        self.assertEqual(discovery.select_method("S1", inputs["S1"], cfg)[0], HIFI_ONT)

    def test_ont_raw_and_filtered_copy_detected(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        names = ont_names(20)
        write_fastq(f"{d['ont']}/raw.fastq.gz", names, [60000] * 20, seed=3)
        write_fastq(f"{d['ont']}/filtered.fastq.gz", names[::2], [60000] * 10, seed=4)
        _, _, issues = self.run_discover(samples={"S1": {"acknowledge": {
            "raw.fastq.gz": {"prefiltered_ont": "t"}, "filtered.fastq.gz": {"prefiltered_ont": "t"}}}})
        self.assertIn("ONT_DUPLICATE_READS", codes(issues))

    def test_ont_excluded_by_human_selects_hifi_only_with_reason(self):
        d = mkdirs(self.tmp, "S1", ("hifi", "ont"))
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        write_fastq(f"{d['ont']}/bad.fastq.gz", ont_names(5), [1000] * 5)
        cfg, inputs, issues = self.run_discover(
            samples={"S1": {"exclude_inputs": {"bad.fastq.gz": "failed flow cell"}}})
        self.assertEqual(codes(issues), [])
        method, reason = discovery.select_method("S1", inputs["S1"], cfg)
        self.assertEqual(method, HIFI_ONLY)
        self.assertIn("excluded by human decision", reason)

    # -- duplicate / derived HiFi data ---------------------------------------
    def test_identical_copies_detected(self):
        d = mkdirs(self.tmp, "S1")
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 10), [12000] * 10)
        shutil.copy(f"{d['hifi']}/a.fastq.gz", f"{d['hifi']}/a_copy.fastq.gz")
        _, _, issues = self.run_discover()
        self.assertIn("DUPLICATE_COPY", codes(issues))

    def test_derived_fastx_with_shared_movie_detected(self):
        d = mkdirs(self.tmp, "S1")
        write_fastq(f"{d['hifi']}/run1.fastq.gz", hifi_names(M1, 10), [12000] * 10, seed=1)
        write_fasta(f"{d['hifi']}/combined.fasta.gz", hifi_names(M1, 10) + hifi_names(M2, 10),
                    [12000] * 20, seed=2)
        _, _, issues = self.run_discover()
        self.assertIn("HIFI_POSSIBLE_DERIVED_COPY", codes(issues))

    @unittest.skipUnless(HAVE_SAMTOOLS, "samtools not available")
    def test_bam_and_derived_fasta_require_explicit_decision(self):
        d = mkdirs(self.tmp, "S1")
        write_bam(f"{d['hifi']}/{M1}.hifi_reads.bc2027.bam", hifi_names(M1, 10), [12000] * 10, movie=M1)
        write_fasta(f"{d['hifi']}/S1.hifi.fasta.gz", hifi_names(M1, 10), [12000] * 10)
        _, _, issues = self.run_discover()
        err = [i for i in issues if i.code == "HIFI_POSSIBLE_DERIVED_COPY"]
        self.assertEqual(len(err), 1)
        self.assertIn("S1.hifi.fasta.gz", err[0].action)   # suggests excluding the derived file
        _, inputs, issues = self.run_discover(
            samples={"S1": {"exclude_inputs": {"S1.hifi.fasta.gz": "derived from the BAM"}}})
        self.assertEqual(codes(issues), [])
        status = {os.path.basename(i.path): i.status for i in inputs["S1"]}
        self.assertEqual(status["S1.hifi.fasta.gz"], "excluded")

    @unittest.skipUnless(HAVE_SAMTOOLS, "samtools not available")
    def test_two_bams_two_movies_ok_and_pbi_ignored(self):
        d = mkdirs(self.tmp, "S1")
        write_bam(f"{d['hifi']}/a.hifi_reads.bam", hifi_names(M1, 5), [12000] * 5, movie=M1)
        write_bam(f"{d['hifi']}/b.hifi_reads.bam", hifi_names(M2, 5), [12000] * 5, movie=M2)
        write_bytes(f"{d['hifi']}/a.hifi_reads.bam.pbi", b"PBI\x01")
        _, inputs, issues = self.run_discover()
        self.assertEqual(codes(issues), [])
        st = sorted((os.path.basename(i.path), i.status) for i in inputs["S1"])
        self.assertIn(("a.hifi_reads.bam.pbi", "ignored"), st)
        self.assertEqual(sum(1 for _, s in st if s == "use"), 2)

    @unittest.skipUnless(HAVE_SAMTOOLS, "samtools not available")
    def test_subread_bam_and_fail_reads_rejected(self):
        d = mkdirs(self.tmp, "S1")
        write_bam(f"{d['hifi']}/x.subreads.bam", [f"{M1}/1/0_100"], [100], readtype="SUBREAD", movie=M1)
        write_bam(f"{d['hifi']}/{M2}.fail_reads.bam", hifi_names(M2, 3), [9000] * 3, movie=M2)
        _, _, issues = self.run_discover()
        msgs = [i.message for i in issues if i.code == "NOT_HIFI"]
        self.assertEqual(len(msgs), 2)

    def test_non_pacbio_names_need_acknowledgement(self):
        d = mkdirs(self.tmp, "S1")
        write_fastq(f"{d['hifi']}/renamed.fastq.gz", [f"read{i}" for i in range(5)], [12000] * 5)
        _, _, issues = self.run_discover()
        self.assertIn("HIFI_IDENTITY_UNVERIFIED", codes(issues))
        _, _, issues = self.run_discover(
            samples={"S1": {"acknowledge": {"renamed.fastq.gz": {"unverified_hifi": "checked"}}}})
        self.assertEqual(codes(issues), [])

    def test_unknown_files_and_unexpected_dirs(self):
        d = mkdirs(self.tmp, "S1")
        write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 5), [12000] * 5)
        write_bytes(f"{d['hifi']}/notes.docx", b"x")
        os.makedirs(os.path.join(self.tmp, "data", "S1", "illumina"))
        _, _, issues = self.run_discover()
        self.assertIn("UNKNOWN_FILE", codes(issues))
        self.assertIn("UNEXPECTED_ENTRY", codes(issues))

    def test_cross_sample_duplicate(self):
        for s in ("S1", "S2"):
            d = mkdirs(self.tmp, s)
            write_fastq(f"{d['hifi']}/a.fastq.gz", hifi_names(M1, 5), [12000] * 5, seed=7)
        _, _, issues = self.run_discover()
        self.assertIn("CROSS_SAMPLE_DUPLICATE", codes(issues))

    def test_external_inputs_and_determinism(self):
        ext = os.path.join(self.tmp, "elsewhere")
        os.makedirs(ext)
        f = write_fasta(f"{ext}/EB2_X.hifi.fasta.gz", hifi_names(M1, 5), [12000] * 5)
        external = {"X": {"hifi": [f]}}
        cfg, inputs, issues = self.run_discover(external=external)
        self.assertEqual(codes(issues), [])
        self.assertEqual(inputs["X"][0].source, "external")
        a = discovery.inventory_records("X", inputs["X"], cfg)
        _, inputs2, _ = self.run_discover(external=external)
        b = discovery.inventory_records("X", inputs2["X"], cfg)
        self.assertEqual(a, b)
        self.assertEqual(discovery.stable_view(a), discovery.stable_view(b))

    def test_relative_external_path_rejected(self):
        _, _, issues = self.run_discover(external={"X": {"hifi": ["relative/path.fa.gz"]}})
        self.assertIn("BAD_EXTERNAL", codes(issues))


if __name__ == "__main__":
    unittest.main()
