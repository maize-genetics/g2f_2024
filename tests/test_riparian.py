"""Tests for g2f_riparian (riparian synteny tables and plots from AnchorWave anchors; synthetic data only)."""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
import g2f_riparian as rip

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_riparian.py")
try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

REF = {"chr1": 1_000_000, "chr2": 800_000, "chr3": 600_000, "scaf_1": 5_000}
HEADER = "refChr\treferenceStart\treferenceEnd\tqueryChr\tqueryStart\tqueryEnd\tstrand\tgene\tblockIndex\tscore\n"


def anchor(rc, rs, q, qs, strand, gene, block=1, ln=900):
    return [rc, rs, rs + ln, q, qs, qs + ln, strand, gene, block]


def collinear(rc, r0, q, q0, n, strand="+", step=10_000, tag="g"):
    out = []
    for i in range(n):
        rs = r0 + i * step if strand == "+" else r0 + (n - 1 - i) * step
        out.append(anchor(rc, rs, q, q0 + i * step, strand, f"{tag}{rc}_{q}_{i}"))
    return out


def write_anchors(path, rows):
    with open(path, "w") as fh:
        fh.write(HEADER + "#block begin\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\t1\n")
            fh.write(f"{r[0]}\t{r[2] + 1}\t{r[2] + 50}\t{r[3]}\t{r[5] + 1}\t{r[5] + 50}\t{r[6]}\tinteranchor\t{r[8]}\tNA\n")


def write_fai(path, lens):
    with open(path, "w") as fh:
        for k, v in lens.items():
            fh.write(f"{k}\t{v}\t0\t80\t81\n")


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class SegmentTest(unittest.TestCase):
    P = dict(max_gap=100_000, min_n=3)

    def dicts(self, rows):
        return [dict(rc=r[0], rs=r[1], re=r[2], q=r[3], qs=r[4], qe=r[5], strand=r[6], gene=r[7], block=r[8])
                for r in rows]

    def test_collinear_is_one_ribbon(self):
        rb = rip.segment(self.dicts(collinear("chr1", 1000, "A", 1000, 10)), 100_000, 3)
        self.assertEqual([(r["rc"], r["n"], r["drawn"]) for r in rb], [("chr1", 10, True)])

    def test_single_stray_gene_does_not_cut_a_block(self):
        rows = collinear("chr1", 1000, "A", 1000, 10)
        rows.insert(5, anchor("chr2", 500_000, "A", 1000 + 4 * 10_000 + 5000, "+", "stray"))
        rb = rip.segment(self.dicts(rows), 100_000, 3)
        drawn = [r for r in rb if r["drawn"]]
        self.assertEqual([(r["rc"], r["n"]) for r in drawn], [("chr1", 10)])
        self.assertEqual([(r["rc"], r["n"], r["drawn"]) for r in rb if not r["drawn"]], [("chr2", 1, False)])

    def test_strand_change_splits(self):
        rows = collinear("chr1", 1000, "A", 1000, 6) + collinear("chr1", 61_000, "A", 61_000, 6, strand="-")
        rb = rip.segment(self.dicts(rows), 100_000, 3)
        self.assertEqual([(r["strand"], r["n"]) for r in rb], [("+", 6), ("-", 6)])

    def test_large_reference_jump_splits(self):
        rows = collinear("chr1", 1000, "A", 1000, 5) + collinear("chr1", 700_000, "A", 51_000, 5)
        self.assertEqual(len(rip.segment(self.dicts(rows), 100_000, 3)), 2)


class LiftTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))

    def cmap(self, rows):
        p = os.path.join(self.d, "map.tsv")
        with open(p, "w") as fh:
            fh.write("candidate_record\tlength\toriginal_contig\toriginal_start\toriginal_end\toffset_to_original\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        return p

    def test_anchor_moves_to_piece_and_is_rebased(self):
        m = rip.read_coordinate_map(self.cmap([("X_p1", 40_000, "X", 1, 40_000, 0),
                                                ("X_p2", 60_000, "X", 40_001, 100_000, 40_000)]))
        a = [dict(rc="chr1", rs=1, re=901, q="X", qs=50_001, qe=50_901, strand="+", gene="g", block="1"),
             dict(rc="chr1", rs=1, re=901, q="Y", qs=1, qe=901, strand="+", gene="h", block="1")]
        out, straddle = rip.lift(a, {"X_p1": 40_000, "X_p2": 60_000, "Y": 5_000}, m)
        self.assertEqual([(x["q"], x["qs"], x["qe"]) for x in out], [("X_p2", 10_001, 10_901), ("Y", 1, 901)])
        self.assertEqual(straddle, 0)

    def test_straddling_anchor_is_counted(self):
        m = rip.read_coordinate_map(self.cmap([("X_p1", 40_000, "X", 1, 40_000, 0),
                                                ("X_p2", 60_000, "X", 40_001, 100_000, 40_000)]))
        a = [dict(rc="chr1", rs=1, re=901, q="X", qs=39_500, qe=40_400, strand="+", gene="g", block="1")]
        out, straddle = rip.lift(a, {"X_p1": 40_000, "X_p2": 60_000}, m)
        self.assertEqual((out[0]["q"], straddle), ("X_p1", 1))

    def test_gapped_map_refused(self):
        with self.assertRaises(SystemExit):
            rip.read_coordinate_map(self.cmap([("X_p1", 40_000, "X", 1, 40_000, 0),
                                               ("X_p2", 59_000, "X", 41_001, 100_000, 41_000)]))

    def test_unknown_contig_refused(self):
        a = [dict(rc="chr1", rs=1, re=901, q="Z", qs=1, qe=901, strand="+", gene="g", block="1")]
        with self.assertRaises(SystemExit):
            rip.lift(a, {"Y": 5_000}, None)


class RunTest(unittest.TestCase):
    """A five-record assembly covering every record class, run through the CLI."""

    @classmethod
    def setUpClass(cls):
        d = cls.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        write_fai(os.path.join(d, "ref.fa.fai"), REF)
        lens = {"A": 300_000, "B": 300_000, "C": 500_000, "D": 400_000, "E": 50_000}
        rows = collinear("chr2", 10_000, "A", 10_000, 20)                                # chr2, forward
        rows += collinear("chr1", 600_000, "B", 10_000, 10, strand="-")                   # chr1, reverse majority
        rows += collinear("chr1", 200_000, "B", 110_000, 4, tag="inv")                    # ... with an inversion
        rows += collinear("chr1", 100_000, "C", 10_000, 15) + collinear("chr3", 50_000, "C", 260_000, 15)  # 2 chr
        rows += collinear("chr3", 300_000, "D", 10_000, 15) + collinear("chr1", 900_000, "D", 360_000, 3)  # terminal
        rows += [anchor("scaf_1", 10, "D", 200_000, "+", "scaf")]
        write_anchors(os.path.join(d, "anchors"), rows)
        write_fai(os.path.join(d, "asm.fa.fai"), lens)
        with open(os.path.join(d, "names.tsv"), "w") as fh:
            fh.write("v1_id\tv1_1_id\nA\tLINEA_chr02_ctg001\n")
        with open(os.path.join(d, "units.tsv"), "w") as fh:
            fh.write("line\tunit\tgroup\tlabel\tanchors\tfai\tname_map\textra\n")
            fh.write(f"LINEA\tLINEA\tg1\tsynthetic\t{d}/anchors\t{d}/asm.fa.fai\t{d}/names.tsv\tkept\n")
        cls.out = os.path.join(d, "out")
        args = ["run", "--units", os.path.join(d, "units.tsv"), "--ref-fai", os.path.join(d, "ref.fa.fai"),
                "--out", cls.out, "--max-gap-bp", "100000", "--sec-min-span-bp", "50000", "--terminal-bp", "50000"]
        if not HAVE_MPL:
            args.append("--no-plots")
        cls.proc = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
        cls.seq = {r["record"]: r for r in read(os.path.join(cls.out, "sequences.tsv"))} if cls.proc.returncode == 0 else {}

    def test_ran(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr)
        self.assertFalse(os.path.exists(self.out + ".partial"))

    def test_classes(self):
        s = self.seq
        self.assertEqual(s["A"]["status"], "one chromosome, collinear")
        self.assertEqual(s["A"]["display_name"], "LINEA_chr02_ctg001")
        self.assertEqual((s["B"]["dominant_chr"], s["B"]["display_reversed"], s["B"]["minority_strand_ribbons"]),
                         ("chr1", "yes", "1"))
        self.assertEqual(s["B"]["orientation_vs_ref"], "mixed")
        self.assertEqual((s["C"]["candidate_interchromosomal"], s["C"]["chr_substantial"]), ("yes", "chr1,chr3"))
        self.assertEqual((s["D"]["candidate_interchromosomal"], s["D"]["off_chr_terminal_blocks"]), ("minor", "1"))
        self.assertEqual((s["E"]["plotted"], s["E"]["status"]), ("no", "no anchors"))

    def test_dominant_counts_and_fraction(self):
        self.assertEqual((self.seq["D"]["dominant_chr"], self.seq["D"]["anchors_on_dominant"],
                          self.seq["D"]["frac_anchors_on_dominant"]), ("chr3", "15", "0.8333"))

    def test_order_by_dominant_chromosome_then_position(self):
        order = sorted((int(r["plot_order"]), k) for k, r in self.seq.items() if r["plot_order"])
        self.assertEqual([k for _, k in order], ["C", "B", "A", "D"])  # chr1 @100k-pos, chr1 @600k, chr2, chr3

    def test_summary_and_inputs(self):
        u = read(os.path.join(self.out, "unit_summary.tsv"))[0]
        self.assertEqual((u["records_plotted"], u["candidate_interchromosomal"], u["minor_interchromosomal"],
                          u["records_no_anchors"], u["anchors_on_other_ref_seqs"], u["n_chr1"]),
                         ("4", "1", "1", "1", "1", "2"))
        kinds = {r["kind"] for r in read(os.path.join(self.out, "INPUTS.tsv"))}
        self.assertEqual(kinds, {"anchors", "fai", "name_map", "ref_fai", "units"})
        with open(os.path.join(self.out, "RUN.json")) as fh:
            run = json.load(fh)
        self.assertEqual(run["schema"], "g2f-riparian/1")
        self.assertEqual(run["units"][0]["extra"], "kept")

    @unittest.skipUnless(HAVE_MPL, "matplotlib not available")
    def test_figures_written(self):
        for f in ("riparian_g1.pdf", "riparian_g1.svg", "riparian_g1.png", "riparian_by_line.pdf",
                  "lines/LINEA.svg"):
            self.assertGreater(os.path.getsize(os.path.join(self.out, "figures", f)), 1000, f)

    def test_existing_output_refused(self):
        p = subprocess.run([sys.executable, SCRIPT, "run", "--units", os.path.join(self.d, "units.tsv"), "--ref-fai",
                            os.path.join(self.d, "ref.fa.fai"), "--out", self.out, "--no-plots"],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("exists", p.stderr)


class AssemblyOrderTest(unittest.TestCase):
    """--order both: grouping by an assignment table, FASTA order, released orientation, whole-group mirroring."""

    @classmethod
    def setUpClass(cls):
        d = cls.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        write_fai(os.path.join(d, "ref.fa.fai"), REF)
        # FASTA order: Z (chr1, '+'), Y (chr3, '-'), X (chr1, '-' majority, synteny-sorted first), W (chr3, '-'),
        # V (no anchors, assigned chr2), U (chr2 anchors, labelled chr1 in the release)
        lens = {"Z": 300_000, "Y": 200_000, "X": 300_000, "W": 200_000, "V": 50_000, "U": 200_000}
        rows = collinear("chr1", 500_000, "Z", 10_000, 12)
        rows += collinear("chr3", 10_000, "Y", 10_000, 10, strand="-")
        rows += collinear("chr1", 10_000, "X", 10_000, 12, strand="-")
        rows += collinear("chr3", 300_000, "W", 10_000, 10, strand="-")
        rows += collinear("chr2", 10_000, "U", 10_000, 10)
        write_anchors(os.path.join(d, "anchors"), rows)
        write_fai(os.path.join(d, "asm.fa.fai"), lens)
        with open(os.path.join(d, "assign.tsv"), "w") as fh:
            fh.write("# provisional\nrec\tchrom\nZ\tchr1\nY\tchr3\nX\tchr1\nW\tchr3\nV\tchr2\nU\tchr1\n")
        with open(os.path.join(d, "units.tsv"), "w") as fh:
            fh.write("line\tunit\tgroup\tlabel\tanchors\tfai\tassignment\tassignment_record_col\tassignment_chr_col\n")
            fh.write(f"LINEB\tLINEB\tg1\tsynthetic\t{d}/anchors\t{d}/asm.fa.fai\t{d}/assign.tsv\trec\tchrom\n")
        with open(os.path.join(d, "units_noassign.tsv"), "w") as fh:
            fh.write("line\tunit\tgroup\tlabel\tanchors\tfai\n")
            fh.write(f"LINEB\tLINEB\tg1\tsynthetic\t{d}/anchors\t{d}/asm.fa.fai\n")
        cls.out = os.path.join(d, "out")
        args = ["run", "--units", os.path.join(d, "units.tsv"), "--ref-fai", os.path.join(d, "ref.fa.fai"),
                "--out", cls.out, "--max-gap-bp", "100000", "--order", "both", "--compare", "LINEB"]
        if not HAVE_MPL:
            args = [a for a in args if a not in ("--compare", "LINEB")] + ["--no-plots"]
        cls.proc = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
        cls.seq = {r["record"]: r for r in read(os.path.join(cls.out, "sequences.tsv"))} if cls.proc.returncode == 0 else {}

    def test_ran(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr)

    def test_assembly_order_is_fasta_order_within_assigned_groups(self):
        order = sorted((int(r["assembly_order"]), k) for k, r in self.seq.items() if r["assembly_order"])
        # chr1: Z, X, U (FASTA order; U is labelled chr1); chr2: V; chr3 mirrored (all '-'): W, Y
        self.assertEqual([k for _, k in order], ["Z", "X", "U", "V", "W", "Y"])
        # the synteny order is unchanged: X (chr1 @ 10 kb) before Z (chr1 @ 500 kb)
        self.assertLess(int(self.seq["X"]["plot_order"]), int(self.seq["Z"]["plot_order"]))

    def test_group_mirroring_only_for_a_reverse_group(self):
        self.assertEqual((self.seq["Y"]["group_reversed"], self.seq["W"]["group_reversed"]), ("yes", "yes"))
        self.assertEqual((self.seq["Z"]["group_reversed"], self.seq["X"]["group_reversed"]), ("no", "no"))
        self.assertEqual(self.seq["X"]["display_reversed"], "yes")  # synteny mode still flips X itself

    def test_assignment_disagreement_and_unanchored_assigned_record(self):
        self.assertEqual((self.seq["U"]["assigned_chr"], self.seq["U"]["assignment_agrees"]), ("chr1", "no"))
        self.assertEqual((self.seq["V"]["assembly_group"], self.seq["V"]["plotted"]), ("chr2", "no"))

    @unittest.skipUnless(HAVE_MPL, "matplotlib not available")
    def test_figures_written(self):
        for f in ("riparian_g1.svg", "asmorder_g1.svg", "asmorder_by_line.pdf", "lines_asmorder/LINEB.png",
                  "compare/compare_LINEB.svg"):
            self.assertGreater(os.path.getsize(os.path.join(self.out, "figures", f)), 1000, f)

    def test_assembly_order_needs_assignment(self):
        p = subprocess.run([sys.executable, SCRIPT, "check", "--units", os.path.join(self.d, "units_noassign.tsv"),
                            "--ref-fai", os.path.join(self.d, "ref.fa.fai"), "--order", "assembly"],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("assignment", p.stderr)


if __name__ == "__main__":
    unittest.main()
