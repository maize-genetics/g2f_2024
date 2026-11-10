"""Tests for g2f_seqid_alias.py + g2f_seqid_alias_verify.py (synthetic sequences only; needs bgzip)."""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import util
import telomere_data as TD

ALI = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_seqid_alias.py")
VER = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_seqid_alias_verify.py")
IDS = ["SAMPLE_A-1_chr01", "SAMPLE_A-1_chr02", "contig_000017_unplaced", "ctg9"]
SEQ = {"SAMPLE_A-1_chr01": "ACGTNNNNNacgt" * 9, "SAMPLE_A-1_chr02": "GGGCCCAT" * 11, "contig_000017_unplaced": "ATATCG" * 7,
       "ctg9": "TTTTGA" * 5}


def run(script, *args):
    env = dict(os.environ, G2F_BGZIP=TD.BGZIP)
    return subprocess.run([sys.executable, script] + [str(x) for x in args], capture_output=True, text=True, env=env)


@unittest.skipUnless(TD.HAVE_TOOLS, "bgzip needed")
class Alias(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        raw = f"{cls.d}/S.fa"
        with open(raw, "w") as fh:
            for i in IDS:
                fh.write(f">{i} some description\n" + "\n".join(SEQ[i][k:k + 20] for k in range(0, len(SEQ[i]), 20)) + "\n")
        subprocess.run([TD.BGZIP, "-f", raw], check=True)
        cls.src = raw + ".gz"
        cls.fa, cls.map = f"{cls.d}/S.alias.fa.gz", f"{cls.d}/S.alias.tsv"
        r = run(ALI, "make", "--fasta", cls.src, "--code", "SA1", "--out-fasta", cls.fa, "--out-map", cls.map, "--threads", "1")
        assert r.returncode == 0, r.stderr

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def rows(self):
        with open(self.map) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_aliases(self):
        self.assertEqual([r["alias"] for r in self.rows()], ["SA1_chr01", "SA1_chr02", "SA1_u0001", "SA1_u0002"])
        with gzip.open(self.fa, "rt") as fh:
            self.assertEqual([l[1:].strip() for l in fh if l.startswith(">")], [r["alias"] for r in self.rows()])

    def test_verify_pass_and_tamper(self):
        r = run(VER, "--source", self.src, "--alias", self.fa, "--map", self.map, "--out", f"{self.d}/v.json")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        bad = f"{self.d}/bad.fa"
        with gzip.open(self.fa, "rt") as fh:
            t = fh.read()
        with open(bad, "w") as fh:
            fh.write(t.replace("GGGCCCAT", "GGGCCCAA", 1))
        r = run(VER, "--source", self.src, "--alias", bad, "--map", self.map, "--out", f"{self.d}/v2.json")
        self.assertEqual(r.returncode, 1)
        with open(f"{self.d}/v2.json") as fh:
            self.assertTrue(any(f.startswith("A2") for f in json.load(fh)["failures"]))

    def test_translate_gff3_and_rmout(self):
        g = f"{self.d}/e.gff3"
        with open(g, "w") as fh:
            fh.write("##gff-version 3\n##sequence-region SA1_chr02 1 88\n"
                     "SA1_chr02\tEDTA\tLTR\t5\t60\t.\t+\t.\tID=TE_1;Name=SA1_u0001_like\n"
                     "SA1_u0001\tEDTA\thelitron\t1\t20\t.\t-\t.\tID=TE_2\n")
        out = f"{self.d}/e.orig.gff3"
        r = run(ALI, "translate", "--map", self.map, "--in", g, "--out", out, "--format", "gff3")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out) as fh:
            t = fh.read()
        self.assertIn("##sequence-region SAMPLE_A-1_chr02 1 88", t)
        self.assertIn("SAMPLE_A-1_chr02\tEDTA\tLTR\t5\t60", t)
        self.assertIn("contig_000017_unplaced\tEDTA\thelitron", t)
        self.assertIn("Name=SA1_u0001_like", t)                      # attributes untouched
        rmo = f"{self.d}/e.out"
        with open(rmo, "w") as fh:
            fh.write("   SW  perc perc perc  query      position\n\n"
                     "  463  10.1  0.0  0.0  SA1_chr01        1   60 (57) + TE1  LTR/Gypsy  1 60 (0) 1\n")
        r = run(ALI, "translate", "--map", self.map, "--in", rmo, "--out", f"{self.d}/e.orig.out", "--format", "rmout")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(f"{self.d}/e.orig.out") as fh:
            line = fh.read().splitlines()[2].split()
        self.assertEqual(line[4], "SAMPLE_A-1_chr01")
        self.assertEqual(line[5:7], ["1", "60"])

    def test_refusals(self):
        g = f"{self.d}/u.gff3"
        with open(g, "w") as fh:
            fh.write("SA1_chr99\tEDTA\tLTR\t1\t2\t.\t+\t.\tID=x\n")
        r = run(ALI, "translate", "--map", self.map, "--in", g, "--out", f"{self.d}/u.out.gff3", "--format", "gff3")
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(os.path.exists(f"{self.d}/u.out.gff3"))
        r = run(ALI, "make", "--fasta", self.src, "--code", "TOOLONGCODE", "--out-fasta", f"{self.d}/x.fa.gz",
                "--out-map", f"{self.d}/x.tsv")
        self.assertNotEqual(r.returncode, 0)                            # TOOLONGCODE_chr01 > 13 characters
        r = run(ALI, "make", "--fasta", self.src, "--code", "SA1", "--out-fasta", self.fa, "--out-map", self.map)
        self.assertIn("refusing to overwrite", r.stderr)


if __name__ == "__main__":
    unittest.main()
