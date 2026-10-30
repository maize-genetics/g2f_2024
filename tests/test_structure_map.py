"""Tests for g2f_structure_map (structural-confidence classes from a riparian run + read evidence; synthetic only)."""
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
from test_riparian import HAVE_MPL, REF, collinear, write_anchors, write_fai

RIP = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_riparian.py")
SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_structure_map.py")


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def inversion_record(q):
    return (collinear("chr1", 100_000, q, 10_000, 20) + collinear("chr1", 300_000, q, 210_000, 6, strand="-")
            + collinear("chr1", 360_000, q, 270_000, 10))


class StructureMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        d = cls.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        write_fai(os.path.join(d, "ref.fa.fai"), REF)
        rows = inversion_record("P")                                                           # read-supported inv
        rows += (collinear("chr3", 10_000, "R", 10_000, 15) + collinear("chr3", 160_000, "R", 160_000, 5, strand="-")
                 + collinear("chr3", 210_000, "R", 210_000, 8))                                # contradicted inv
        rows += collinear("chr2", 10_000, "Q", 10_000, 15) + collinear("chr3", 500_000, "Q", 170_000, 4)  # join
        rows += collinear("chr2", 10_000, "S", 10_000, 5)                                      # redundant with Q
        write_anchors(os.path.join(d, "a.anchors"), rows)
        write_fai(os.path.join(d, "a.fai"), {"P": 400_000, "R": 300_000, "Q": 210_000, "S": 60_000, "T": 500_000})
        write_anchors(os.path.join(d, "b.anchors"), inversion_record("P2"))                    # same inversion
        write_fai(os.path.join(d, "b.fai"), {"P2": 400_000})
        with open(os.path.join(d, "units.tsv"), "w") as fh:
            fh.write("line\tunit\tgroup\tlabel\tanchors\tfai\n")
            fh.write(f"LINEA\tLINEA\tg1\tsynthetic\t{d}/a.anchors\t{d}/a.fai\n")
            fh.write(f"LINEB\tLINEB\tg1\tsynthetic\t{d}/b.anchors\t{d}/b.fai\n")
        cls.rip = os.path.join(d, "rip")
        p = subprocess.run([sys.executable, RIP, "run", "--units", os.path.join(d, "units.tsv"), "--ref-fai",
                            os.path.join(d, "ref.fa.fai"), "--out", cls.rip, "--max-gap-bp", "100000",
                            "--no-plots"], capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        with open(os.path.join(d, "ev.tsv"), "w") as fh:
            fh.write("unit\trecord\tstart\tend\tkind\tcall\tsource\tdetail\n")
            fh.write("LINEA\tP\t200900\t210000\tjunction\tSUPPORTED\ttest\tleft boundary\n")
            fh.write("LINEA\tP\t260900\t270000\tjunction\tCONTINUOUS\ttest\tright boundary\n")
            fh.write("LINEA\tR\t150900\t160000\tjunction\tCONTRADICTED\ttest\tsplit reads to a partner\n")
            fh.write("LINEA\tQ\t150900\t170000\tjunction\tBREAK\ttest\tno spanning read\n")
            fh.write("LINEA\tQ\t150900\t170000\tedit\tretained\ttest\tdocumented\n")
        with gzip.open(os.path.join(d, "win.tsv.gz"), "wt") as fh:
            fh.write("contig\tstart\tend\tdepth_ratio\tclass\tmultimap_frac\tmixed_per_kb\n")
            for s in range(0, 400_000, 10_000):
                fh.write(f"P\t{s}\t{s + 10_000}\t1.0\tnormal\t0\t0\n")
            for s in range(0, 60_000, 10_000):
                fh.write(f"S\t{s}\t{s + 10_000}\t0.5\thalf\t0\t0\n")
        with open(os.path.join(d, "windows.tsv"), "w") as fh:
            fh.write(f"unit\tpath\nLINEA\t{d}/win.tsv.gz\n")
        cls.out = os.path.join(d, "out")
        args = ["run", "--riparian", cls.rip, "--evidence", os.path.join(d, "ev.tsv"), "--windows",
                os.path.join(d, "windows.tsv"), "--out", cls.out, "--terminal-bp", "50000", "--gap-bp", "100000",
                "--slop-bp", "5000"]
        if not HAVE_MPL:
            args.append("--no-plots")
        cls.proc = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
        cls.reg = read(os.path.join(cls.out, "regions.tsv")) if cls.proc.returncode == 0 else []

    def one(self, unit, record, cls_):
        hits = [r for r in self.reg if r["unit"] == unit and r["record"] == record and r["region_class"] == cls_]
        self.assertEqual(len(hits), 1, (unit, record, cls_, hits))
        return hits[0]

    def test_ran(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr)
        self.assertFalse(os.path.exists(self.out + ".partial"))

    def test_read_supported_inversion_is_a_reference_difference(self):
        r = self.one("LINEA", "P", "inversion")
        self.assertEqual(r["read_evidence"], "read-supported")
        self.assertIn("reference difference", r["interpretation"])
        self.assertTrue(r["recurrence"].startswith("2 lines"))
        self.assertEqual(r["depth_median"], "1.0")

    def test_contradicted_inversion_is_an_error_candidate(self):
        r = self.one("LINEA", "R", "inversion")
        self.assertEqual(r["read_evidence"], "read-contradicted")
        self.assertIn("assembly error candidate", r["interpretation"])

    def test_unspanned_terminal_join_is_not_called_contradicted(self):
        r = self.one("LINEA", "Q", "off_chromosome")
        self.assertEqual((r["subtype"], r["read_evidence"]), ("terminal", "not-spanned"))
        self.assertIn("absence of support", r["interpretation"])
        self.assertIn("retained", r["evidence_detail"])

    def test_untested_recurrent_inversion(self):
        r = self.one("LINEB", "P2", "inversion")
        self.assertEqual(r["read_evidence"], "untested")
        self.assertIn("other lines", r["interpretation"])

    def test_redundant_and_no_anchor(self):
        r = self.one("LINEA", "S", "redundant")
        self.assertEqual((r["partner"], r["depth_median"]), ("Q", "0.5"))
        self.one("LINEA", "T", "no_anchor")

    def test_summary_and_run_record(self):
        s = {r["unit"]: r for r in read(os.path.join(self.out, "unit_summary.tsv"))}
        self.assertEqual((s["LINEA"]["inversion_n"], s["LINEA"]["off_chromosome_n"], s["LINEA"]["redundant_n"]),
                         ("2", "1", "1"))
        self.assertEqual((s["LINEA"]["junctions_contradicted"], s["LINEA"]["junctions_not_spanned"]), ("1", "1"))
        with open(os.path.join(self.out, "RUN.json")) as fh:
            self.assertEqual(json.load(fh)["schema"], "g2f-structure-map/1")

    @unittest.skipUnless(HAVE_MPL, "matplotlib not available")
    def test_figures(self):
        for f in ("structure_g1.svg", "structure_g1.png", "structure_by_unit.pdf"):
            self.assertGreater(os.path.getsize(os.path.join(self.out, "figures", f)), 1000, f)

    def test_unknown_record_in_evidence_refused(self):
        bad = os.path.join(self.d, "bad.tsv")
        with open(bad, "w") as fh:
            fh.write("unit\trecord\tstart\tend\tkind\tcall\nLINEA\tNOPE\t1\t2\tjunction\tSUPPORTED\n")
        p = subprocess.run([sys.executable, SCRIPT, "check", "--riparian", self.rip, "--evidence", bad],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("not a record", p.stderr)


if __name__ == "__main__":
    unittest.main()
