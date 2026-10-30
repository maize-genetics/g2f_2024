"""Tests for g2f_dup_pairs (same-haplotype duplicate candidates; synthetic sequences, real minimap2 + samtools)."""
import csv
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_dup_pairs.py")


def find(name, fallback):
    p = shutil.which(name)
    return p or (fallback if os.path.exists(fallback) else None)


MINIMAP2 = find("minimap2", "/programs/minimap2-2.30/minimap2")
SAMTOOLS = find("samtools", "/programs/samtools-1.20/bin/samtools")
HAVE_TOOLS = bool(MINIMAP2 and SAMTOOLS)


def mutate(seq, rate, rng):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for k, v in recs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 60):
                fh.write(v[i:i + 60] + "\n")
    with open(path + ".fai", "w") as fh:  # samtools-compatible index (60 bp lines)
        off = 0
        for k, v in recs.items():
            off += len(k) + 2
            fh.write(f"{k}\t{len(v)}\t{off}\t60\t61\n")
            off += len(v) + (len(v) + 59) // 60


def anchors(path, rows):
    with open(path, "w") as fh:
        fh.write("refChr\treferenceStart\treferenceEnd\tqueryChr\tqueryStart\tqueryEnd\tstrand\tgene\tblockIndex\t"
                 "score\n")
        for rec, pos, gene in rows:
            fh.write(f"chr1\t{1000 * int(gene[1:])}\t{1000 * int(gene[1:]) + 500}\t{rec}\t{pos}\t{pos + 500}\t+\t"
                     f"{gene}\t1\t1\n")


@unittest.skipUnless(HAVE_TOOLS, "minimap2/samtools not available")
class DupPairsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(7)
        d = cls.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        hap1_allele = "".join(rng.choice("ACGT") for _ in range(120_000))
        hap2_allele = mutate(hap1_allele, 0.01, rng)               # the other allele, 99 % identical
        flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        B = flank(80_000) + hap1_allele + flank(80_000)            # hap1's own copy, in a big record
        A = hap2_allele                                            # the other allele, also kept in hap1
        P = "".join(rng.choice("ACGT") for _ in range(100_000))    # unrelated sequence (paralogue-like)
        C = flank(50_000) + hap2_allele + flank(50_000)            # hap2's copy
        D = flank(40_000) + hap1_allele + flank(40_000)            # a hap2 copy of B's allele (other allele only in hap1)
        write_fasta(os.path.join(d, "hap1.fa"), {"B": B, "A": A, "P": P})
        write_fasta(os.path.join(d, "hap2.fa"), {"C": C, "D": D})
        g = [f"g{i}" for i in range(1, 11)]
        anchors(os.path.join(d, "hap1.anchors"),
                [("B", 80_000 + 10_000 * i, x) for i, x in enumerate(g)] + [("A", 10_000 * i, x) for i, x in enumerate(g)]
                + [("P", 10_000 * i, x) for i, x in enumerate(g[:5])])
        anchors(os.path.join(d, "hap2.anchors"), [("C", 50_000 + 10_000 * i, x) for i, x in enumerate(g)]
                + [("D", 40_000 + 10_000 * i, x) for i, x in enumerate(g)])
        with gzip.open(os.path.join(d, "win.tsv.gz"), "wt") as fh:
            fh.write("contig\tstart\tend\tdepth_ratio\tclass\n")
            for s in range(0, 120_000, 10_000):
                fh.write(f"A\t{s}\t{s + 10_000}\t0.5\thalf\n")
            for s in range(0, 280_000, 10_000):  # B: reads split with A over the shared allele, full depth outside
                half = 80_000 <= s < 200_000
                fh.write(f"B\t{s}\t{s + 10_000}\t{0.5 if half else 1.0}\t{'half' if half else 'normal'}\n")
        with open(os.path.join(d, "pairs.tsv"), "w") as fh:
            fh.write("pair_id\tline\tunit\tfasta\tanchors\trecord_a\trecord_b\twindows\tother_unit\tother_fasta\t"
                     "other_anchors\tother_records\n")
            fh.write(f"X1\tLINEA\tLINEA.hap1\t{d}/hap1.fa\t{d}/hap1.anchors\tA\tB\t{d}/win.tsv.gz\tLINEA.hap2\t"
                     f"{d}/hap2.fa\t{d}/hap2.anchors\tC\n")
            fh.write(f"X2\tLINEA\tLINEA.hap1\t{d}/hap1.fa\t{d}/hap1.anchors\tP\tB\t\t\t\t\t\n")
            fh.write(f"X3\tLINEA\tLINEA.hap1\t{d}/hap1.fa\t{d}/hap1.anchors\tA\tB\t{d}/win.tsv.gz\tLINEA.hap2\t"
                     f"{d}/hap2.fa\t{d}/hap2.anchors\tD\n")
        cls.out = os.path.join(d, "out")
        cls.proc = subprocess.run([sys.executable, SCRIPT, "run", "--pairs", os.path.join(d, "pairs.tsv"), "--out",
                                   cls.out, "--threads", "1", "--minimap2", MINIMAP2, "--samtools", SAMTOOLS,
                                   "--flank-bp", "20000", "--min-aln-bp", "1000"], capture_output=True, text=True)
        cls.rows = {}
        if cls.proc.returncode == 0:
            with open(os.path.join(cls.out, "dup_pairs.tsv")) as fh:
                cls.rows = {r["pair_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def test_ran(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "work")))  # extracted sequences removed

    def test_allelic_copy(self):
        r = self.rows["X1"]
        self.assertEqual((r["genes_shared"], r["shared_frac_of_a"]), ("10", "1.0"))
        self.assertGreater(float(r["AC_gc_identity"]), 0.999)
        self.assertLess(float(r["AB_gc_identity"]), 0.995)
        self.assertIn("other haplotype's allele", r["classification"])
        self.assertEqual((r["confidence"], r["depth_a_median"]), ("high", "0.5"))

    def test_both_alleles_in_one_haplotype_by_depth_split(self):
        r = self.rows["X3"]
        self.assertGreater(float(r["BC_gc_identity"]), 0.999)  # hap2 holds a copy of B's allele, not A's
        self.assertEqual((r["depth_b_median"], r["depth_b_flank_median"]), ("0.5", "1.0"))
        self.assertIn("both alleles of a heterozygous locus", r["classification"])
        self.assertEqual(r["confidence"], "high")
        self.assertTrue(r["other_haplotype_state"].startswith("carries a copy of B"))

    def test_reclassify_writes_a_new_directory(self):
        out = os.path.join(self.d, "re")
        p = subprocess.run([sys.executable, SCRIPT, "reclassify", "--from", self.out, "--out", out],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(out, "dup_pairs.tsv")) as fh:
            again = {r["pair_id"]: r["classification"] for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(again, {k: v["classification"] for k, v in self.rows.items()})

    def test_unrelated_sequence_is_paralogous_anchoring(self):
        r = self.rows["X2"]
        self.assertTrue(r["classification"].startswith("paralogous anchoring"), r)

    def test_plan_runs_nothing(self):
        p = subprocess.run([sys.executable, SCRIPT, "plan", "--pairs", os.path.join(self.d, "pairs.tsv")],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("X1\tAC", p.stdout)


if __name__ == "__main__":
    unittest.main()
