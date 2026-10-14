"""Real-tool compatibility test (opt-in): the complete workflow on a tiny
synthetic genome (300 kb) with the REAL BioHPC tools from
config/tools/biohpc.yaml: pbtk 3.5.0 (container), hifiasm 0.25.0-r726,
meryl + Merqury 1.3 (container), BUSCO 6.1.0 (container), seqkit, samtools.

    G2F_REALTOOLS=1 bash tests/run_tests.sh test_realtools     # ~2-4 min, 2 threads

The pinned production lineage (poales_odb12.2) is not installed yet, so this
test pins the installed poales_odb12 release instead; it validates the BUSCO
interface, offline lineage access, dataset pinning and output parsing, not
the production lineage choice. Synthetic random sequence gives BUSCO 0%.
"""
import os
import shutil
import tempfile
import unittest

import yaml
from realtools_data import build
from test_workflow import HIFI_ONT, WorkflowCase, have_snakemake, unprotect
from util import PIPELINE_DIR, SAMTOOLS  # noqa: F401  (also puts scripts/ on sys.path)

LINEAGES = os.path.expanduser("~/busco_downloads/lineages")
CONTAINERS = os.path.join(os.path.dirname(PIPELINE_DIR), "containers")
ENABLED = os.environ.get("G2F_REALTOOLS") == "1"


def real_tools():
    with open(os.path.join(PIPELINE_DIR, "config", "tools", "biohpc.yaml")) as fh:
        t = yaml.safe_load(fh)
    t["containers_dir"] = CONTAINERS
    return t


@unittest.skipUnless(ENABLED and have_snakemake() and os.path.isdir(f"{LINEAGES}/poales_odb12")
                     and os.path.exists(f"{CONTAINERS}/pbtk-3.5.0.sif"),
                     "set G2F_REALTOOLS=1 (needs BioHPC tools, pbtk container, poales_odb12)")
class RealToolsTest(WorkflowCase):
    tools = staticmethod(real_tools)
    lineage = ("poales_odb12", None)

    def setUp(self):
        base = os.environ.get("TMPDIR") or tempfile.gettempdir()
        self._parent = tempfile.mkdtemp(prefix="g2f_real_", dir=base)
        tempfile.tempdir, old = self._parent, tempfile.tempdir
        try:
            super().setUp()
        finally:
            tempfile.tempdir = old

    def tearDown(self):
        unprotect(self._parent)
        shutil.rmtree(self._parent)

    def build_data(self):
        src = build(os.path.join(self.tmp, "synthetic"), samtools=SAMTOOLS)
        hifi = os.path.join(self.tmp, "data", "SYN", "hifi")
        ont = os.path.join(self.tmp, "data", "SYN", "ont")
        os.makedirs(hifi)
        os.makedirs(ont)
        os.rename(f"{src}/hifi.fastq.gz", f"{hifi}/SYN_A.hifi.fastq.gz")
        os.rename(f"{src}/hifi_subset.bam", f"{hifi}/SYN_B.hifi_reads.bam")
        os.rename(f"{src}/ont.fastq.gz", f"{ont}/SYN.ont.fastq.gz")

    def adjust_science(self):
        self.science["qc"]["busco"]["lineages"] = {"poales_odb12": {
            "datasets_version": "odb12", "creation_date": "2026-05-22", "number_of_BUSCOs": 6282}}

    def lineage_dir(self):
        return LINEAGES

    def test_full_workflow_with_real_tools(self):
        self.write_yaml("resources.yaml", {"tmpdir": f"{self.tmp}/tmp", "rules": {
            "default": {"threads": 2, "mem_mb": 4000}}})
        r = self.smk(cores="2")
        self.assertEqual(r.returncode, 0, r.stdout[-5000:] + r.stderr[-8000:])
        row = self.summary_rows()["SYN"]
        self.assertEqual(row["method_id"], HIFI_ONT)
        self.assertEqual(row["hifi_reads"], "600")             # 500 FASTQ + 100 converted by pbtk
        self.assertEqual(row["ont_input_type"], "raw")
        self.assertEqual(row["ont_raw_reads"], "35")
        self.assertEqual(row["ont_retained_reads"], "25")         # >= 50 kb
        self.assertEqual(row["hifiasm_version"], "0.25.0-r726")
        self.assertEqual(row["busco_lineage"], "poales_odb12")
        self.assertEqual(row["busco_datasets_version"], "odb12")
        self.assertEqual(row["busco_n"], "6282")
        self.assertGreater(int(row["total_bp"]), 290000)      # one ~300 kb contig expected
        self.assertLessEqual(int(row["contigs"]), 3)
        self.assertGreater(float(row["merqury_qv"]), 20)
        self.assertGreater(float(row["merqury_kmer_completeness_pct"]), 90)
        self.assertEqual(row["alt_contigs"], "0")              # homozygous: empty a_ctg
        print("\nreal-tool summary:", {k: row[k] for k in (
            "contigs", "total_bp", "contig_n50", "merqury_qv", "merqury_kmer_completeness_pct",
            "busco_C", "kmer_main_peak", "kmer_assessment", "ont_retained_reads")})


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(ENABLED and os.path.exists("/programs/hifiasm-0.25.0/hifiasm"),
                     "set G2F_REALTOOLS=1 (needs hifiasm 0.25.0)")
class HifiasmMultipleUlFilesTest(unittest.TestCase):
    """hifiasm 0.25.0 must consume every file of a comma-separated --ul list
    (the workflow passes one file per ONT input instead of a merged copy)."""

    def test_comma_separated_ul_files_are_all_consumed(self):
        import subprocess
        from realtools_data import genome, sample_reads, write_fastq
        from g2f.validate import parse_hifiasm_log
        hifiasm = "/programs/hifiasm-0.25.0/hifiasm"
        base = os.environ.get("TMPDIR") or tempfile.gettempdir()
        tmp = tempfile.mkdtemp(prefix="g2f_ul_", dir=base)
        try:
            g = genome()
            hifi = sample_reads(g, 600, 15000, seed=1)
            write_fastq(f"{tmp}/hifi.fastq.gz", [f"m84109_240101_000000_s1/{i}/ccs"
                                                 for i in range(1, 601)], hifi)
            ont = sample_reads(g, 25, 75000, seed=2)
            write_fastq(f"{tmp}/a.fastq.gz", [f"a{i}" for i in range(10)], ont[:10], q="5")
            write_fastq(f"{tmp}/b.fastq.gz", [f"b{i}" for i in range(15)], ont[10:], q="5")
            counts = {}
            for name, ul in (("A", "a"), ("B", "b"), ("AB", "a,b")):
                d = f"{tmp}/{name}"
                os.makedirs(d)
                files = ",".join(f"{tmp}/{x}.fastq.gz" for x in ul.split(","))
                with open(f"{d}/err.log", "w") as err:
                    subprocess.run([hifiasm, "-o", "T", "-t", "2", "-l0", "--ul", files,
                                    f"{tmp}/hifi.fastq.gz"], cwd=d, stderr=err,
                                   stdout=subprocess.DEVNULL, check=True)
                counts[name] = parse_hifiasm_log(f"{d}/err.log").get("ul_passed_reads")
            self.assertEqual(counts, {"A": 10, "B": 15, "AB": 25})
        finally:
            shutil.rmtree(tmp)
