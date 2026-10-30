"""Synthetic tests for scripts/g2f_fasta_gz_verify.py."""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "g2f_fasta_gz_verify.py")
CONTIGS = b">a\nACGTACGT\n>b\nGGGTTT\n"
AGP_OK = ("##agp-version 2.1\n"
          "s1\t1\t8\t1\tW\ta\t1\t8\t+\n"
          "s1\t9\t18\t2\tU\t10\tscaffold\tyes\tmap\n"
          "s1\t19\t24\t3\tW\tb\t1\t6\t-\n")
SCAF = b">s1\nACGTACGT" + b"N" * 10 + b"AAACCC\n"


class T(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, name, data, gz=False):
        p = os.path.join(self.d, name)
        with (gzip.open if gz else open)(p, "wb") as fh:
            fh.write(data)
        return p

    def run_tool(self, *args):
        r = subprocess.run([sys.executable, TOOL, *args], capture_output=True, text=True)
        return r.returncode, json.loads(r.stdout)

    def test_source_pass(self):
        src, gz = self.put("c.fa", CONTIGS), self.put("c.fa.gz", CONTIGS, gz=True)
        rc, res = self.run_tool("source", "--gz", gz, "--source", src)
        self.assertEqual((rc, res["status"], res["records"], res["bp"]), (0, "PASS", 2, 14))

    def test_source_detects_changed_base(self):
        src, gz = self.put("c.fa", CONTIGS), self.put("c.fa.gz", CONTIGS.replace(b"GGGTTT", b"GGGTTA"), gz=True)
        rc, res = self.run_tool("source", "--gz", gz, "--source", src)
        self.assertEqual((rc, res["status"]), (1, "FAIL"))

    def test_source_detects_wrong_recorded_sha(self):
        src, gz = self.put("c.fa", CONTIGS), self.put("c.fa.gz", CONTIGS, gz=True)
        rc, _ = self.run_tool("source", "--gz", gz, "--source", src, "--source-sha256", "0" * 64)
        self.assertEqual(rc, 1)

    def test_source_detects_duplicate_ids(self):
        dup = CONTIGS + b">a\nAC\n"
        src, gz = self.put("c.fa", dup), self.put("c.fa.gz", dup, gz=True)
        rc, res = self.run_tool("source", "--gz", gz, "--source", src)
        self.assertEqual(rc, 1)
        self.assertIn("duplicate record IDs", res["failures"])

    def test_source_detects_stale_index(self):
        src, gz = self.put("c.fa", CONTIGS), self.put("c.fa.gz", CONTIGS, gz=True)
        self.put("c.fa.gz.fai", b"a\t8\t3\t8\t9\nb\t7\t15\t6\t7\n")
        rc, res = self.run_tool("source", "--gz", gz, "--source", src)
        self.assertEqual(rc, 1)
        self.assertIn(".fai does not match the records", res["failures"])

    def test_agp_pass_with_reverse_component(self):
        comp = self.put("c.fa", CONTIGS)
        agp = self.put("s.agp", AGP_OK.encode())
        gz = self.put("s.fa.gz", SCAF, gz=True)
        rc, res = self.run_tool("agp", "--gz", gz, "--agp", agp, "--components", comp)
        self.assertEqual((rc, res["status"], res["components_used"]), (0, "PASS", 2))

    def test_agp_detects_wrong_orientation(self):
        comp = self.put("c.fa", CONTIGS)
        agp = self.put("s.agp", AGP_OK.replace("1\t6\t-", "1\t6\t+").encode())
        gz = self.put("s.fa.gz", SCAF, gz=True)
        rc, _ = self.run_tool("agp", "--gz", gz, "--agp", agp, "--components", comp)
        self.assertEqual(rc, 1)

    def test_agp_detects_unused_component(self):
        comp = self.put("c.fa", CONTIGS + b">x\nAC\n")
        agp = self.put("s.agp", AGP_OK.encode())
        gz = self.put("s.fa.gz", SCAF, gz=True)
        rc, res = self.run_tool("agp", "--gz", gz, "--agp", agp, "--components", comp)
        self.assertEqual(rc, 1)
        rc, res = self.run_tool("agp", "--gz", gz, "--agp", agp, "--components", comp, "--allow-unused")
        self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main()
