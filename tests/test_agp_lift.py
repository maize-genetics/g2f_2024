"""Tests for g2f_agp_lift (AGP-based coordinate conversion; synthetic sequences and features only)."""
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
import g2f_agp_lift as L

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_agp_lift.py")
RC = str.maketrans("ACGT", "TGCA")


def lines(path):
    with open(path) as fh:
        return fh.readlines()


def rc(s):
    return s.translate(RC)[::-1]


class LiftTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        rng = random.Random(7)
        self.ctg = {n: "".join(rng.choice("ACGT") for _ in range(ln))
                    for n, ln in (("ctgA", 3000), ("ctgB", 2000), ("ctgC", 1500), ("ctgU", 800))}
        # scaffold1 = ctgA(+) gap100 ctgB(-) gap100 ctgC(+); ctgU unplaced (its own object)
        s1 = self.ctg["ctgA"] + "N" * 100 + rc(self.ctg["ctgB"]) + "N" * 100 + self.ctg["ctgC"]
        self.scf = {"UNIT_chr01": s1, "ctgU": self.ctg["ctgU"]}
        a, b, c = 3000, 2000, 1500
        agp = ["##agp-version 2.1",
               f"UNIT_chr01\t1\t{a}\t1\tW\tctgA\t1\t{a}\t+",
               f"UNIT_chr01\t{a + 1}\t{a + 100}\t2\tU\t100\tscaffold\tyes\tmap",
               f"UNIT_chr01\t{a + 101}\t{a + 100 + b}\t3\tW\tctgB\t1\t{b}\t-",
               f"UNIT_chr01\t{a + b + 101}\t{a + b + 200}\t4\tU\t100\tscaffold\tyes\tmap",
               f"UNIT_chr01\t{a + b + 201}\t{a + b + 200 + c}\t5\tW\tctgC\t1\t{c}\t+",
               "ctgU\t1\t800\t1\tW\tctgU\t1\t800\t+"]
        self.agp = self.w("u.agp", "\n".join(agp) + "\n")
        self.cfa = self.w("c.fa", "".join(f">{k}\n{v}\n" for k, v in self.ctg.items()))
        self.sfa = self.w("s.fa", "".join(f">{k}\n{v}\n" for k, v in self.scf.items()))

    def w(self, name, text):
        p = os.path.join(self.d, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def run_cli(self, *args):
        return subprocess.run([sys.executable, SCRIPT] + list(args), capture_output=True, text=True)

    def test_transforms_match_sequence(self):
        A = L.AGP(self.agp)
        for comp, s, e in (("ctgA", 1, 50), ("ctgB", 1, 10), ("ctgB", 1990, 2000), ("ctgC", 700, 900), ("ctgU", 5, 9)):
            (obj, os_, oe, st), why = A.to_object(comp, s, e, "+")
            self.assertIsNone(why)
            piece = self.scf[obj][os_ - 1:oe]
            want = self.ctg[comp][s - 1:e]
            self.assertEqual(piece, rc(want) if st == "-" else want, (comp, s, e))
            self.assertEqual(A.to_component(obj, os_, oe, st)[0], (comp, s, e, "+"))
        self.assertEqual(A.to_object("ctgB", 1, 10, "+")[0][3], "-")   # strand flips on a '-' component

    def test_gap_and_multicomponent_refused(self):
        A = L.AGP(self.agp)
        self.assertEqual(A.to_component("UNIT_chr01", 2990, 3010)[1], "GAP_SPANNING")
        self.assertEqual(A.to_component("UNIT_chr01", 3005, 3050)[1], "IN_GAP")
        self.assertEqual(A.to_object("nope", 1, 2)[1], "NOT_IN_AGP")

    def test_gff_round_trip_and_hierarchy(self):
        gff = ["##gff-version 3",
               "ctgB\tt\tgene\t100\t900\t.\t+\t.\tID=g1",
               "ctgB\tt\tmRNA\t100\t900\t.\t+\t.\tID=m1;Parent=g1",
               "ctgB\tt\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1",
               "ctgB\tt\tCDS\t500\t900\t.\t+\t2\tID=c1;Parent=m1",
               "ctgU\tt\tgene\t10\t20\t.\t-\t.\tID=g2"]
        src = self.w("in.gff3", "\n".join(gff) + "\n")
        out, back, un = (os.path.join(self.d, x) for x in ("o.gff3", "b.gff3", "u.tsv"))
        r = self.run_cli("gff", "--agp", self.agp, "--to", "object", "-i", src, "-o", out, "--unmapped", un)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = [l.split("\t") for l in lines(out) if not l.startswith("#")]
        g1 = [x for x in rows if "ID=g1" in x[8]][0]
        self.assertEqual((g1[0], g1[6], int(g1[3]), int(g1[4])), ("UNIT_chr01", "-", 3000 + 100 + (2000 - 900) + 1,
                                                                   3000 + 100 + (2000 - 100) + 1))
        self.assertEqual({x[7] for x in rows if x[2] == "CDS"}, {"0", "2"})   # phase unchanged
        self.assertIn("##sequence-region UNIT_chr01 1 6700\n", lines(out))
        r = self.run_cli("gff", "--agp", self.agp, "--to", "component", "-i", out, "-o", back)
        self.assertEqual(r.returncode, 0, r.stderr)
        key = lambda f: sorted(tuple(l.rstrip("\n").split("\t")) for l in lines(f) if not l.startswith("#"))
        self.assertEqual(key(src), key(back))   # exact round trip, IDs unchanged

    def test_gap_spanning_gene_dropped_whole(self):
        gff = ["##gff-version 3",
               "UNIT_chr01\tt\tgene\t2900\t3200\t.\t+\t.\tID=gx",
               "UNIT_chr01\tt\texon\t2900\t2950\t.\t+\t.\tID=ex1;Parent=gx",
               "UNIT_chr01\tt\texon\t3150\t3200\t.\t+\t.\tID=ex2;Parent=gx",
               "UNIT_chr01\tt\tgene\t10\t20\t.\t+\t.\tID=ok"]
        src = self.w("s.gff3", "\n".join(gff) + "\n")
        out, un = os.path.join(self.d, "c.gff3"), os.path.join(self.d, "un.tsv")
        r = self.run_cli("gff", "--agp", self.agp, "--to", "component", "-i", src, "-o", out, "--unmapped", un)
        self.assertEqual(r.returncode, 0, r.stderr)
        kept = [l for l in lines(out) if not l.startswith("#")]
        self.assertEqual(len(kept), 1)
        self.assertIn("ID=ok", kept[0])
        reasons = [l.split("\t")[0] for l in lines(un)][1:]
        self.assertEqual(sorted(reasons), ["GAP_SPANNING", "GROUP_MEMBER_UNMAPPED", "GROUP_MEMBER_UNMAPPED"])

    def test_bed(self):
        src = self.w("i.bed", "ctgB\t0\t10\tf\t0\t+\nctgA\t0\t5\tg\t0\t-\n")
        out = os.path.join(self.d, "o.bed")
        self.assertEqual(self.run_cli("bed", "--agp", self.agp, "--to", "object", "-i", src, "-o", out).returncode, 0)
        rows = [l.split() for l in lines(out)]
        self.assertEqual(rows[0], ["UNIT_chr01", str(3100 + 2000 - 10), str(3100 + 2000), "f", "0", "-"])
        self.assertEqual(rows[1], ["UNIT_chr01", "0", "5", "g", "0", "-"])
        self.assertNotEqual(self.run_cli("bed", "--agp", self.agp, "--to", "object",
                                         "-i", self.w("x.bed", "a\t0\t1\tn\t0\t+\t0\t1\t0\t1\t1\t0\n"),
                                         "-o", out).returncode, 0)   # BED12 refused

    def test_check_with_chain(self):
        chain = os.path.join(self.d, "u.chain")
        self.assertEqual(self.run_cli("chain", "--agp", self.agp, "-o", chain).returncode, 0)
        r = self.run_cli("check", "--agp", self.agp, "--component-fasta", self.cfa, "--object-fasta", self.sfa,
                         "--chain", chain, "--n", "500")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        res = json.loads(r.stdout)
        self.assertEqual(res["status"], "PASS")
        self.assertGreater(res["tests"]["intervals_minus"], 0)
        self.assertEqual(res["tests"]["gap_spanning_refused"], 2)

    def test_check_detects_corruption(self):
        bad = dict(self.scf)
        s = bad["UNIT_chr01"]
        bad["UNIT_chr01"] = s[:3500] + ("A" if s[3500] != "A" else "C") + s[3501:]
        sfa = self.w("bad.fa", "".join(f">{k}\n{v}\n" for k, v in bad.items()))
        r = self.run_cli("check", "--agp", self.agp, "--component-fasta", self.cfa, "--object-fasta", sfa, "--n", "50")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("reconstruction differs", r.stdout)


if __name__ == "__main__":
    unittest.main()
