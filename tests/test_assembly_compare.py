"""Tests for g2f_assembly_compare (contiguity / completeness table; synthetic files only)."""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
import g2f_assembly_compare as cmp

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_assembly_compare.py")


class ContiguityTest(unittest.TestCase):
    def test_n50_ng50_and_size_classes(self):
        lens = sorted([60_000_000, 30_000_000, 5_000_000, 3_000_000, 2_000_000], reverse=True)  # 100 Mb
        c = cmp.contiguity(lens, 200_000_000)
        self.assertEqual((c["N50"], c["L50"]), (60_000_000, 1))
        self.assertEqual((c["NG50"], c["LG50"]), (2_000_000, 5))  # all 100 Mb needed to reach half of 200
        self.assertEqual((c["n_ge_10Mb"], c["frac_bp_ge_10Mb"], c["n_ge_50Mb"]), (2, 0.9, 1))

    def test_ng50_undefined_when_assembly_too_small(self):
        self.assertEqual(cmp.contiguity([10, 5], 100)["NG50"], 0)


class RunTest(unittest.TestCase):
    def test_table(self):
        d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        with open(os.path.join(d, "a.fai"), "w") as fh:
            fh.write("c1\t70\t0\t60\t61\nc2\t30\t0\t60\t61\n")
        with open(os.path.join(d, "busco.txt"), "w") as fh:
            fh.write("# The lineage dataset is: poales_odb12.2 (Creation date: x)\n"
                     "\tC:98.9%[S:73.9%,D:25.0%],F:0.1%,M:1.0%,n:6188,E:2.2%\n")
        with open(os.path.join(d, "m.qv"), "w") as fh:
            fh.write("asm.hap1\t10\t100\t55.5\t1e-6\nBoth\t20\t200\t54.0\t2e-6\n")
        with open(os.path.join(d, "m.cn"), "w") as fh:
            fh.write("asm.hap1\tall\t80\t100\t80.0\nboth\tall\t97\t100\t97.0\n")
        with open(os.path.join(d, "asm.tsv"), "w") as fh:
            fh.write("name\tline\tassembly_class\trole\tfai\tbusco\tmerqury_qv\tmerqury_qv_row\tmerqury_cn\t"
                     "merqury_cn_row\tmerqury_cn_pair_row\tnotes\n")
            fh.write(f"L.hap1\tL\tphased_haplotype\trelease_target\t{d}/a.fai\t{d}/busco.txt\t{d}/m.qv\tasm.hap1\t"
                     f"{d}/m.cn\tasm.hap1\tboth\tx\n")
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, SCRIPT, "run", "--assemblies", os.path.join(d, "asm.tsv"),
                            "--genome-size", "200", "--out", out, "--no-plots"], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(out, "assembly_compare.tsv")) as fh:
            r = next(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual((r["N50"], r["NG50"], r["busco_D"], r["merqury_qv"], r["merqury_completeness"],
                          r["merqury_pair_completeness"]), ("70", "30", "25.0", "55.5", "80.0", "97.0"))
        p = subprocess.run([sys.executable, SCRIPT, "run", "--assemblies", os.path.join(d, "asm.tsv"), "--out", out],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)  # never overwrites


if __name__ == "__main__":
    unittest.main()
