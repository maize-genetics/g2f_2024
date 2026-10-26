"""Tests for g2f_release_fasta (record reordering + independent verification; synthetic FASTA only)."""
import json
import os
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
import g2f_release_fasta as rel

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_release_fasta.py")

# mixed wrapping (80, unwrapped, 3) and lower case must survive byte for byte
RECORDS = [
    ("LINEA_chrUn_ctg002", "ACGT" * 30, 80),
    ("LINEA_chr10_ctg001", "acgtNNacgt" * 5, 0),
    ("LINEA_chr2_ctg010", "GATTACA" * 9, 3),
    ("LINEA_chr2_ctg002", "TTTT", 80),
    ("LINEA_chrUn_ctg010", "CCCCGGGG", 80),
    ("LINEA_chr01_ctg001", "A" * 161, 80),
    ("LINEA_chrUn_ctg001", "G", 80),
]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq, w in records:
            fh.write(f">{rid}\n")
            lines = [seq] if not w else [seq[i:i + w] for i in range(0, len(seq), w)]
            fh.write("\n".join(lines) + "\n")


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


class SortKeyTest(unittest.TestCase):
    def test_natural_order(self):
        ids = ["X_chr10_ctg1", "X_chr2_ctg10", "X_chr2_ctg2", "X_chrUn_ctg3", "X_chr1_ctg1", "X_chrUn_ctg20"]
        self.assertEqual(sorted(ids, key=rel.sort_key),
                         ["X_chr1_ctg1", "X_chr2_ctg2", "X_chr2_ctg10", "X_chr10_ctg1", "X_chrUn_ctg3", "X_chrUn_ctg20"])

    def test_prefix_may_start_with_digit_or_contain_underscore(self):
        self.assertEqual(rel.sort_key("3LINE_chr01_ctg001"), (0, 1, 1))
        self.assertEqual(rel.sort_key("A_B_chrUn_ctg007"), (1, 0, 7))

    def test_unmatched_refused(self):
        with self.assertRaises(ValueError):
            rel.sort_key("ptg000001l")
        self.assertEqual(rel.sort_key("ptg000001l", unmatched_ok=True), (2, 0, 0))


class SortVerifyTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        self.src = os.path.join(self.d, "src.fa")
        write_fasta(self.src, RECORDS)
        self.dst = os.path.join(self.d, "dst.fa")

    def test_sort_then_verify(self):
        r = run("sort", "--in", self.src, "--out", self.dst, "--order-tsv", os.path.join(self.d, "o.tsv"))
        self.assertEqual(r.returncode, 0, r.stderr)
        got = [ln[1:].strip() for ln in open(self.dst) if ln.startswith(">")]
        self.assertEqual(got, ["LINEA_chr01_ctg001", "LINEA_chr2_ctg002", "LINEA_chr2_ctg010", "LINEA_chr10_ctg001",
                               "LINEA_chrUn_ctg001", "LINEA_chrUn_ctg002", "LINEA_chrUn_ctg010"])
        self.assertEqual(sorted(open(self.src).read().split(">")), sorted(open(self.dst).read().split(">")))
        j = os.path.join(self.d, "v.json")
        r = run("verify", "--src", self.src, "--dst", self.dst, "--json", j)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(json.load(open(j))["passed"])

    def test_verify_detects_changed_base(self):
        run("sort", "--in", self.src, "--out", self.dst)
        txt = open(self.dst).read()
        self.assertIn("acgtNNacgt", txt)
        open(self.dst, "w").write(txt.replace("acgtNNacgt", "acgtNAacgt", 1))
        r = run("verify", "--src", self.src, "--dst", self.dst)
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(json.loads(r.stdout)["checks"]["sequence_sha256_identical"])

    def test_verify_detects_rewrap_and_wrong_order(self):
        recs = [(i, s, 60 if w == 80 else w) for i, s, w in RECORDS]   # same bases, different wrapping
        write_fasta(self.dst, recs)                                     # and unsorted
        out = json.loads(run("verify", "--src", self.src, "--dst", self.dst).stdout)
        self.assertTrue(out["checks"]["sequence_sha256_identical"])
        self.assertFalse(out["checks"]["record_blocks_byte_identical"])
        self.assertFalse(out["checks"]["dst_in_expected_sort_order"])

    def test_verify_detects_duplicate_and_missing(self):
        write_fasta(self.dst, RECORDS[:-1] + [RECORDS[0]])
        out = json.loads(run("verify", "--src", self.src, "--dst", self.dst).stdout)
        self.assertFalse(out["checks"]["each_src_id_exactly_once_in_dst"])
        self.assertFalse(out["passed"])

    def test_sort_refuses_duplicates(self):
        write_fasta(self.src, RECORDS + [RECORDS[0]])
        self.assertNotEqual(run("sort", "--in", self.src, "--out", self.dst).returncode, 0)
        self.assertFalse(os.path.exists(self.dst))


if __name__ == "__main__":
    unittest.main()
