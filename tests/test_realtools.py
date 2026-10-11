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
from util import PIPELINE_DIR, SAMTOOLS

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
        base = os.environ.get("TMPDIR") or "/workdir/" + os.environ.get("USER", "") + "/tmp"
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
