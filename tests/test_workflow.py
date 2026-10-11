"""Workflow-level tests on tiny synthetic data.

* Dry-run DAG resolution (no tools executed).
* A complete end-to-end execution in which seqkit, samtools, md5sum and awk
  are the real tools, while hifiasm, pbtk bam2fastq, meryl, Merqury and BUSCO
  are replaced by stubs (tests/stubs/). This tests plumbing, validation gates,
  provenance and resumability; it does NOT validate the scientific tools.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml
from util import (HAVE_SAMTOOLS, HAVE_SEQKIT, PIPELINE_DIR, SAMTOOLS, SEQKIT, hifi_names,
                  mkdirs, ont_names, write_bam, write_bytes, write_fasta, write_fastq)

STUBS = os.path.join(PIPELINE_DIR, "tests", "stubs")
HIFI_ONT = "hifiasm-0.25.0_hifi-ont-ul50k_l0_v1"
HIFI_ONLY = "hifiasm-0.25.0_hifi_l0_v1"
M1, M2, M3 = "m84109_240806_024441_s4", "m84109_240611_191050_s2", "m64001_220101_000000"


def have_snakemake():
    r = subprocess.run([sys.executable, "-m", "snakemake", "--version"], capture_output=True)
    return r.returncode == 0


@unittest.skipUnless(HAVE_SAMTOOLS and HAVE_SEQKIT and have_snakemake(),
                     "needs samtools, seqkit and snakemake (module load python/3.12.7)")
class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_wf_")
        t = self.tmp
        s1 = mkdirs(t, "S1", ("hifi", "ont"))
        write_bam(f"{s1['hifi']}/{M1}.hifi_reads.bc2027.bam", hifi_names(M1, 10), [3000] * 10, movie=M1)
        write_fastq(f"{s1['hifi']}/{M2}.hifi.fastq.gz", hifi_names(M2, 10), [3500] * 10, seed=2)
        write_fastq(f"{s1['ont']}/fc1.fastq.gz", ont_names(12, seed=1), [4000, 52000, 61000] * 4, seed=3)
        write_fastq(f"{s1['ont']}/fc2.fastq.gz", ont_names(6, seed=2), [70000, 1000] * 3, seed=4)
        s2 = mkdirs(t, "S2", ("hifi",))
        write_fasta(f"{s2['hifi']}/S2.hifi.fasta.gz", hifi_names(M3, 12), [4000] * 12, seed=5)
        lin = os.path.join(t, "lineages", "poales_odb12")
        os.makedirs(lin)
        write_bytes(f"{lin}/dataset.cfg", b"name=poales_odb12\nOrthoDB_version=12.1\n"
                                          b"creation_date=2026-05-22\nnumber_of_BUSCOs=6282\n")
        tools = {"containers_dir": "", "tools": {
            "hifiasm": {"cmd": f"{STUBS}/hifiasm", "version_cmd": f"{STUBS}/hifiasm --version"},
            "seqkit": {"cmd": SEQKIT, "version_cmd": f"{SEQKIT} version"},
            "samtools": {"cmd": SAMTOOLS, "version_cmd": f"{SAMTOOLS} --version | head -1"},
            "bam2fastq": {"cmd": f"env SAMTOOLS={SAMTOOLS} {STUBS}/bam2fastq",
                          "version_cmd": f"{STUBS}/bam2fastq --version"},
            "meryl": {"cmd": f"{STUBS}/meryl", "version_cmd": f"{STUBS}/meryl --version"},
            "merqury": {"cmd": f"{STUBS}/merqury.sh", "version_cmd": "echo stub"},
            "busco": {"cmd": f"{STUBS}/busco", "version_cmd": f"{STUBS}/busco --version"},
            "md5sum": {"cmd": "md5sum", "version_cmd": "md5sum --version | head -1"},
        }}
        self.samples = {"samples": {"S1": {"aliases": ["EB2_S1"]}}}
        self.write_yaml("tools.yaml", tools)
        self.write_yaml("samples.yaml", self.samples)
        self.write_yaml("local.yaml", {"busco_lineage_dir": os.path.join(t, "lineages")})
        self.write_yaml("config.yaml", {
            "data_dir": f"{t}/data", "results_dir": f"{t}/results", "legacy_dirs": [],
            "samples": f"{t}/samples.yaml", "tools": f"{t}/tools.yaml", "local": f"{t}/local.yaml",
            "tmpdir": f"{t}/tmp"})

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def write_yaml(self, name, obj):
        with open(os.path.join(self.tmp, name), "w") as fh:
            yaml.safe_dump(obj, fh)

    def smk(self, *args):
        cmd = [sys.executable, "-m", "snakemake", "-s", f"{PIPELINE_DIR}/workflow/Snakefile",
               "--directory", self.tmp, "--configfile", f"{self.tmp}/config.yaml",
               "--cores", "2", *args, "--rerun-triggers", "mtime"]
        env = dict(os.environ, TMPDIR=self.tmp)
        return subprocess.run(cmd, capture_output=True, text=True, env=env)

    def test_dry_run_resolves_both_methods(self):
        r = self.smk("--dry-run", "--printshellcmds")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        out = r.stdout + r.stderr
        self.assertIn(f"results/S1/assembly/{HIFI_ONT}/", out)
        self.assertIn(f"results/S2/assembly/{HIFI_ONLY}/", out)
        self.assertIn("--ul", out)                       # ONT sample gets --ul
        self.assertIn("seqkit seq -m 50000", out.replace(SEQKIT, "seqkit"))
        self.assertIn("hifi_bam_to_fastq", out)          # BAM converted, FASTQ/FASTA used in place
        self.assertNotIn(f"results/S2/reads/ont", out)
        self.assertFalse(os.path.exists(f"{self.tmp}/results"), "dry run must not write results")

    def test_invalid_ont_halts_without_fallback(self):
        os.makedirs(f"{self.tmp}/data/S2/ont")
        write_bytes(f"{self.tmp}/data/S2/ont/run.pod5", b"x")
        r = self.smk("--dry-run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("RAW_SIGNAL", r.stdout + r.stderr)

    def test_end_to_end_with_stub_tools(self):
        r = self.smk()
        self.assertEqual(r.returncode, 0, r.stdout[-4000:] + r.stderr[-6000:])
        res = f"{self.tmp}/results"
        with open(f"{res}/summary/assembly_summary.tsv") as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh]
        hdr, body = rows[0], rows[1:]
        recs = {r_[0]: dict(zip(hdr, r_)) for r_ in body}
        self.assertEqual(set(recs), {"S1", "S2"})
        self.assertEqual(recs["S1"]["method_id"], HIFI_ONT)
        self.assertEqual(recs["S2"]["method_id"], HIFI_ONLY)
        self.assertEqual(recs["S1"]["hifi_reads"], "20")
        self.assertEqual(recs["S1"]["hifi_bases"], str(10 * 3000 + 10 * 3500))
        self.assertEqual(recs["S1"]["ont_raw_reads"], "18")
        self.assertEqual(recs["S1"]["ont_retained_reads"], str(8 + 3))   # >= 50 kb only
        self.assertEqual(recs["S2"]["ont_raw_reads"], "")
        self.assertEqual(recs["S1"]["contigs"], "3")
        self.assertEqual(recs["S1"]["busco_lineage"], "poales_odb12")
        self.assertIn("OrthoDB 12.1", recs["S1"]["busco_dataset"])
        self.assertEqual(recs["S1"]["merqury_qv"], "45.2")
        self.assertEqual(recs["S1"]["hifiasm_version"], "0.25.0-r726")
        # Provenance was recorded at execution time with the hifiasm version check.
        with open(f"{res}/S1/provenance/hifiasm/{HIFI_ONT}.json") as fh:
            prov = json.load(fh)
        self.assertEqual(prov["status"], "completed")
        self.assertEqual(prov["tools"]["hifiasm"]["reported"], "0.25.0-r726")
        self.assertIn("--ul", prov["command"])
        self.assertIn("-l0", prov["command"])
        self.assertTrue(glob.glob(f"{res}/S1/inputs/md5/hifi/*.md5"))
        # Assemblies are write-protected.
        fa = f"{res}/S1/assembly/{HIFI_ONT}/S1.{HIFI_ONT}.p_ctg.fa"
        self.assertFalse(os.access(fa, os.W_OK))
        # Resumable: a second invocation has nothing to do.
        r2 = self.smk("--dry-run")
        self.assertEqual(r2.returncode, 0, r2.stderr)
        self.assertIn("Nothing to be done", r2.stdout + r2.stderr)
        # Changing inputs under existing results is refused.
        write_fastq(f"{self.tmp}/data/S2/hifi/{M1}.extra.fastq.gz", hifi_names("m99999_250101_000000", 8),
                    [3000] * 8, seed=9)
        r3 = self.smk("--dry-run")
        self.assertNotEqual(r3.returncode, 0)
        self.assertIn("differ from the inventory", r3.stdout + r3.stderr)
        for p in glob.glob(f"{res}/**", recursive=True):
            if os.path.isfile(p):
                os.chmod(p, 0o644)   # allow tearDown to remove protected files

    def test_kmer_flag_blocks_assembly_until_reviewed(self):
        os.environ["STUB_KMER_HET"] = "1"
        try:
            r = self.smk("--keep-going", f"{self.tmp}/results/S2/kmer/hifi.k21.gate.json")
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("Assembly blocked pending human review", r.stdout + r.stderr)
            self.samples["samples"]["S2"] = {"kmer_review": {
                "status": "accepted", "reason": "test", "reviewer": "test", "date": "2026-10-11"}}
            self.write_yaml("samples.yaml", self.samples)
            r = self.smk(f"{self.tmp}/results/S2/kmer/hifi.k21.gate.json")
            self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        finally:
            del os.environ["STUB_KMER_HET"]


if __name__ == "__main__":
    unittest.main()
