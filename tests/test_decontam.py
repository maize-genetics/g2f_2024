"""Tests for g2f_decontam (synthetic sequences and synthetic FCS reports only).

FCS itself is replaced by stub scripts that write reports in the FCS 0.5.5 formats. The organelle screen runs the
real minimap2 (G2F_MINIMAP2, PATH, or the BioHPC path); those tests are skipped when minimap2 is not available."""
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

import util
import g2f_decontam as D

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_decontam.py")
EXCLUDE = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_release_exclude.py")
MINIMAP2 = os.environ.get("G2F_MINIMAP2") or shutil.which("minimap2") or "/programs/minimap2-2.30/minimap2"
HAVE_MINIMAP2 = os.path.exists(MINIMAP2)


def run(*args, env=None):
    e = dict(os.environ, G2F_MINIMAP2=MINIMAP2, **(env or {}))
    return subprocess.run([sys.executable, SCRIPT] + [str(x) for x in args], capture_output=True, text=True, env=e)


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def jload(p):
    with open(p) as fh:
        return json.load(fh)


def rtext(p):
    with open(p) as fh:
        return fh.read()


def wtext(p, t):
    with open(p, "w") as fh:
        fh.write(t)


def wfasta(p, recs):
    with open(p, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def agp_line(obj, ob, oe, pn, kind, *rest):
    return "\t".join(map(str, (obj, ob, oe, pn, kind) + rest)) + "\n"


rng = random.Random(7)
PT = util.rand_seq(30000, rng)          # synthetic "plastid"
MT = util.rand_seq(40000, rng)          # synthetic "mitochondrion"
NUC = util.rand_seq(120000, rng)        # synthetic nuclear reference
FILL = util.rand_seq(60000, rng)        # sequence in neither reference

# contigs: name -> sequence
CONTIGS = {
    "u_pt": PT[2000:14000],                              # whole record plastid -> CONFIDENT, plastid
    "u_mt": MT[5000:25000],                              # whole record mitochondrion -> CONFIDENT, mitochondrion
    "u_mix": PT[0:10000] + NUC[0:10000],                 # half organelle + nuclear -> MIXED
    "u_ins": NUC[20000:50000][:15000] + PT[15000:16000] + NUC[50000:65000],   # nuclear with a 1 kb insertion -> none
    "u_bact": FILL[0:8000],                              # foreign (FCS stub says EXCLUDE)
    "u_trim": FILL[8000:20000],                          # FCS-GX TRIM (partial) -> review
    "c_a": NUC[70000:100000],                            # placed in chromosome object
    "c_org": MT[0:9000],                                 # placed, organellar -> screened only with --all-records
    "u_odd": PT[20000:28000],                            # unplaced but reverse-oriented in the AGP -> review
}


def build_assembly(d):
    fa = os.path.join(d, "S.contigs.fa")
    wfasta(fa, CONTIGS.items())
    subprocess.run([util.SAMTOOLS, "faidx", fa], check=True)
    L = {k: len(v) for k, v in CONTIGS.items()}
    agp = os.path.join(d, "S.agp")
    with open(agp, "w") as fh:
        fh.write("##agp-version\t2.1\n")
        a, b = L["c_a"], L["c_org"]
        fh.write(agp_line("S_chr01", 1, a, 1, "W", "c_a", 1, a, "+"))
        fh.write(agp_line("S_chr01", a + 1, a + 100, 2, "U", 100, "scaffold", "yes", "align_genus"))
        fh.write(agp_line("S_chr01", a + 101, a + 100 + b, 3, "W", "c_org", 1, b, "+"))
        for n in ("u_pt", "u_mt", "u_mix", "u_ins", "u_bact", "u_trim"):
            fh.write(agp_line(n, 1, L[n], 1, "W", n, 1, L[n], "+"))
        fh.write(agp_line("u_odd", 1, L["u_odd"], 1, "W", "u_odd", 1, L["u_odd"], "-"))
    org = os.path.join(d, "organelles.fa")
    wfasta(org, [("PT1", PT), ("MT1", MT)])
    nuc = os.path.join(d, "nuclear.fa")
    wfasta(nuc, [("N1", NUC)])
    return fa, agp, org, nuc


GX_HEADER = ('##[["FCS genome report", 2, 1], {"git-rev": "v0.5.5-test", "db": {"build-date": "2023-01-24"}, '
             '"run-info": {"asserted-div": "plnt:plants", "inferred-primary-divs": ["plnt:plants"]}}]\n'
             "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n")
TAX_HEADER = ('##[["GX taxonomy analysis report", 3, 1], {}]\n'
              "#seq-id\tseq-len\t(xp,lc,co,n,mt,pt,pm)-len\tcvg-by-all\tsep1\ttax-name-1\n")


def write_reports(adir, gdir, L, extra_gx=""):
    os.makedirs(adir, exist_ok=True)
    os.makedirs(gdir, exist_ok=True)
    wtext(os.path.join(adir, "fcs_adaptor_report.txt"),
          "#accession\tlength\taction\trange\tname\n"
          f"c_a\t{L['c_a']}\tACTION_TRIM\t101..120\tCONTAMINATION_SOURCE_TYPE_ADAPTOR:X:test adaptor\n")
    wtext(os.path.join(gdir, "S.fa.4577.fcs_gx_report.txt"), GX_HEADER +
          f"u_bact\t1\t{L['u_bact']}\t{L['u_bact']}\tEXCLUDE\tprok:a-proteobacteria\t100\tTestbacter\n"
          f"u_trim\t1001\t3000\t{L['u_trim']}\tTRIM\tprok:a-proteobacteria\t80\tTestbacter\n"
          f"c_a\t1\t{L['c_a']}\t{L['c_a']}\tEXCLUDE\tprok:a-proteobacteria\t100\tTestbacter\n"
          f"u_pt\t1\t500\t{L['u_pt']}\tREVIEW\tprok:cyanobacteria\t10\tTestcyano\n" + extra_gx)
    wtext(os.path.join(gdir, "S.fa.4577.taxonomy.rpt"), TAX_HEADER +
          f"u_pt\t{L['u_pt']}\t0,0,0,0,0,{L['u_pt']},0\t1\t|\tZea mays\n"
          f"c_a\t{L['c_a']}\t{L['c_a']},0,0,0,0,0,0\t1\t|\tZea mays\n")


STUB_ADAPTOR = """#!/usr/bin/env bash
set -e
while [[ $# -gt 0 ]]; do case $1 in --fasta-input) F=$2; shift 2;; --output-dir) O=$2; shift 2;;
  --euk|--prok) echo "$1" > "$O/mode.txt"; shift;; *) shift;; esac; done
cp "$STUB_REPORTS/adaptor/fcs_adaptor_report.txt" "$O/"
"""
STUB_FCS = """import os, shutil, sys
a = sys.argv
out = a[a.index("--out-dir") + 1]
open(os.path.join(out, "args.txt"), "w").write(" ".join(a[1:]))
for f in os.listdir(os.path.join(os.environ["STUB_REPORTS"], "gx")):
    shutil.copy(os.path.join(os.environ["STUB_REPORTS"], "gx", f), out)
"""


class Defaults(unittest.TestCase):
    def test_yaml_equals_builtin_defaults(self):
        with open(os.path.join(util.PIPELINE_DIR, "config", "decontam_science.yaml")) as fh:
            self.assertEqual(yaml.safe_load(fh), D.DEFAULTS)


class Units(unittest.TestCase):
    def test_merge_and_class(self):
        self.assertEqual(D.iv_bp([(0, 10), (5, 20), (30, 40)]), 30)
        o = D.DEFAULTS["organelle"]
        self.assertEqual(D.organelle_class(10000, 9500, 0, 500, o), "CONFIDENT_ORGANELLE")
        self.assertEqual(D.organelle_class(100000, 95000, 5000, 5000, o), "ORGANELLE_DOMINANT_MIXED")
        self.assertEqual(D.organelle_class(100000, 95000, None, 5000, o), "ORGANELLE_DOMINANT_MIXED")
        self.assertEqual(D.organelle_class(10000, 6000, 0, 4000, o), "ORGANELLE_DOMINANT_MIXED")
        self.assertEqual(D.organelle_class(10000, 4000, 0, 6000, o), "")

    def test_placement(self):
        with tempfile.TemporaryDirectory() as d:
            fa, agp, _, _ = build_assembly(d)
            p = D.placement(agp, D.fasta_lengths(fa), D.DEFAULTS["chromosome_regex"])
            self.assertEqual(p["c_a"], ("placed", "S_chr01"))
            self.assertEqual(p["c_org"], ("placed", "S_chr01"))
            self.assertEqual(p["u_pt"], ("unplaced", "u_pt"))
            self.assertEqual(p["u_odd"][0], "unplaced_nonstandard")


@unittest.skipUnless(util.HAVE_SAMTOOLS, "samtools not available")
class Fcs(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        self.fa, self.agp, _, _ = build_assembly(self.d)
        self.L = {k: len(v) for k, v in CONTIGS.items()}

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_summarize_classes_placement_and_totals(self):
        r = os.path.join(self.d, "rep")
        write_reports(os.path.join(r, "adaptor"), os.path.join(r, "gx"), self.L)
        out = os.path.join(self.d, "fcs")
        p = run("fcs-summarize", "--fasta", self.fa, "--label", "S", "--agp", self.agp, "--adaptor-dir",
                os.path.join(r, "adaptor"), "--gx-dir", os.path.join(r, "gx"), "--out", out)
        self.assertEqual(p.returncode, 0, p.stderr)
        calls = {(c["seq_id"], c["tool"]): c for c in tsv(os.path.join(out, "fcs_calls.tsv"))}
        self.assertEqual(calls[("u_bact", "FCS-GX")]["call_class"], "foreign_contaminant")
        self.assertEqual(calls[("u_bact", "FCS-GX")]["placement"], "unplaced")
        self.assertEqual(calls[("c_a", "FCS-GX")]["placement"], "placed")
        self.assertEqual(calls[("u_pt", "FCS-GX")]["call_class"], "gx_review")
        self.assertEqual(calls[("c_a", "FCS-adaptor")]["call_class"], "adaptor_vector")
        s = jload(os.path.join(out, "summary.json"))
        # foreign: u_bact whole + u_trim 2000 bp + c_a whole; adaptor 20 bp inside c_a counted once in the union
        self.assertEqual(s["foreign_bp"], self.L["u_bact"] + 2000 + self.L["c_a"])
        self.assertEqual(s["adaptor_bp"], 20)
        self.assertEqual(s["union_flagged_bp"], s["foreign_bp"])
        self.assertEqual(s["gx_organelle_dominant_records"], 1)
        self.assertEqual(s["gx_db_build"], "2023-01-24")
        v = jload(os.path.join(out, "VALIDATED.json"))
        self.assertEqual(v["step"], "fcs")
        self.assertEqual(v["schema"], "g2f-decontam/1")
        # never overwrites
        p = run("fcs-summarize", "--fasta", self.fa, "--label", "S", "--adaptor-dir", os.path.join(r, "adaptor"),
                "--out", out)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("refusing to overwrite", p.stderr)

    def test_wrong_fasta_is_refused(self):
        r = os.path.join(self.d, "rep")
        write_reports(os.path.join(r, "adaptor"), os.path.join(r, "gx"), self.L,
                      extra_gx="ghost\t1\t10\t10\tEXCLUDE\tprok:x\t100\tX\n")
        out = os.path.join(self.d, "fcs")
        p = run("fcs-summarize", "--fasta", self.fa, "--label", "S", "--adaptor-dir", os.path.join(r, "adaptor"),
                "--gx-dir", os.path.join(r, "gx"), "--out", out)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("not in", p.stderr)
        self.assertFalse(os.path.exists(out))

    def test_fcs_runs_both_tools_with_species_and_mode(self):
        r = os.path.join(self.d, "rep")
        write_reports(os.path.join(r, "adaptor"), os.path.join(r, "gx"), self.L)
        fcs = os.path.join(self.d, "fcsdir")
        os.makedirs(fcs)
        wtext(os.path.join(fcs, "run_fcsadaptor.sh"), STUB_ADAPTOR)
        wtext(os.path.join(fcs, "fcs.py"), STUB_FCS)
        db = os.path.join(self.d, "db")
        os.makedirs(db)
        wtext(os.path.join(db, "all.gxi"), "")
        out = os.path.join(self.d, "fcs")
        p = run("fcs", "--fasta", self.fa, "--label", "S", "--agp", self.agp, "--tax-id", 4577, "--gx-db",
                os.path.join(db, "all"), "--fcs-dir", fcs, "--threads", 2, "--out", out, env={"STUB_REPORTS": r})
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("--tax-id 4577", rtext(os.path.join(out, "gx", "args.txt")))
        self.assertEqual(rtext(os.path.join(out, "adaptor", "mode.txt")).strip(), "--euk")
        v = jload(os.path.join(out, "VALIDATED.json"))
        self.assertEqual(v["tax_id"], 4577)
        self.assertEqual(len(v["commands"]), 2)

    def test_fcs_requires_database_unless_skipped(self):
        p = run("fcs", "--fasta", self.fa, "--label", "S", "--tax-id", 4577, "--fcs-dir", self.d,
                "--out", os.path.join(self.d, "x"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("--gx-db", p.stderr)


@unittest.skipUnless(HAVE_MINIMAP2 and util.HAVE_SAMTOOLS, "minimap2 or samtools not available")
class OrganelleAndPlan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.fa, cls.agp, cls.org, cls.nuc = build_assembly(cls.d)
        cls.L = {k: len(v) for k, v in CONTIGS.items()}
        cls.out = os.path.join(cls.d, "org")
        p = run("organelle", "--fasta", cls.fa, "--label", "S", "--agp", cls.agp, "--organelles", cls.org,
                "--plastid-ids", "PT1", "--nuclear", cls.nuc, "--threads", 2, "--out", cls.out)
        assert p.returncode == 0, p.stderr
        cls.rows = {r["record"]: r for r in tsv(os.path.join(cls.out, "organelle_screen.tsv"))}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def test_classes(self):
        r = self.rows
        self.assertEqual(r["u_pt"]["window_class"], "CONFIDENT_ORGANELLE")
        self.assertEqual(r["u_pt"]["organelle_kind"], "ORGANELLE_PLASTID")
        self.assertEqual(r["u_mt"]["window_class"], "CONFIDENT_ORGANELLE")
        self.assertEqual(r["u_mt"]["organelle_kind"], "ORGANELLE_MITOCHONDRION")
        self.assertEqual(r["u_mix"]["window_class"], "ORGANELLE_DOMINANT_MIXED")
        self.assertGreater(int(r["u_mix"]["nonorganelle_nuclear_bp"]), 9000)
        self.assertEqual(r["u_ins"]["window_class"], "")
        self.assertAlmostEqual(int(r["u_ins"]["organelle_bp"]), 1000, delta=100)
        self.assertEqual(r["u_bact"]["window_class"], "")

    def test_placed_records_are_not_screened_by_default(self):
        self.assertNotIn("c_org", self.rows)
        self.assertNotIn("c_a", self.rows)
        out = os.path.join(self.d, "org_all")
        p = run("organelle", "--fasta", self.fa, "--label", "S", "--agp", self.agp, "--organelles", self.org,
                "--all-records", "--threads", 2, "--out", out)
        self.assertEqual(p.returncode, 0, p.stderr)
        rows = {r["record"]: r for r in tsv(os.path.join(out, "organelle_screen.tsv"))}
        self.assertEqual(rows["c_org"]["placement"], "placed")
        self.assertEqual(rows["c_org"]["window_class"], "CONFIDENT_ORGANELLE")
        self.assertEqual(rows["c_org"]["organelle_kind"], "ORGANELLE")      # no --plastid-ids: kind not split
        # without the nuclear test the non-organelle part of u_mix is treated as possibly nuclear: still MIXED
        self.assertEqual(rows["u_mix"]["window_class"], "ORGANELLE_DOMINANT_MIXED")
        self.assertEqual(rows["u_mix"]["nonorganelle_nuclear_bp"], "")

    def test_plan_rows_review_and_release_exclude_accepts_it(self):
        r = os.path.join(self.d, "rep")
        write_reports(os.path.join(r, "adaptor"), os.path.join(r, "gx"), self.L)
        fcs = os.path.join(self.d, "fcs_p")
        p = run("fcs-summarize", "--fasta", self.fa, "--label", "S", "--agp", self.agp, "--adaptor-dir",
                os.path.join(r, "adaptor"), "--gx-dir", os.path.join(r, "gx"), "--out", fcs)
        self.assertEqual(p.returncode, 0, p.stderr)
        out = os.path.join(self.d, "plan")
        p = run("plan", "--fcs", fcs, "--organelle", self.out, "--out", out)
        self.assertEqual(p.returncode, 0, p.stderr)
        plan = {x["record"]: x for x in tsv(os.path.join(out, "EXCLUSION_PLAN.candidate.tsv"))}
        self.assertEqual(set(plan), {"u_pt", "u_mt", "u_bact"})
        self.assertEqual(plan["u_bact"]["reason_class"], "FOREIGN_CONTAMINANT")
        self.assertEqual(plan["u_pt"]["reason_class"], "ORGANELLE_PLASTID")
        self.assertEqual(plan["u_mt"]["reason_class"], "ORGANELLE_MITOCHONDRION")
        self.assertTrue(rtext(os.path.join(out, "EXCLUSION_PLAN.candidate.tsv")).startswith("# CANDIDATE"))
        rev = {(x["record"], x["call"]): x for x in tsv(os.path.join(out, "REVIEW.tsv"))}
        self.assertIn(("u_trim", "foreign_contaminant"), rev)                 # TRIM is partial
        self.assertIn(("u_mix", "ORGANELLE_DOMINANT_MIXED"), rev)
        self.assertIn(("c_a", "FOREIGN_CONTAMINANT"), rev)                    # placed: never a plan row
        self.assertIn("structural edit", rev[("c_a", "FOREIGN_CONTAMINANT")]["reason"])
        self.assertIn(("c_a", "adaptor_vector"), rev)
        self.assertIn(("u_pt", "gx_review"), rev)
        self.assertIn(("u_odd", "ORGANELLE_PLASTID"), rev)                    # reverse-oriented unplaced object
        # the candidate plan is valid input for the exclusion tool
        p = subprocess.run([sys.executable, EXCLUDE, "check-plan", "--plan", os.path.join(out, "EXCLUSION_PLAN.candidate.tsv"),
                            "--agp", self.agp, "--contigs", self.fa], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("OK", p.stdout)

    def test_conflicting_calls_go_to_review(self):
        r = os.path.join(self.d, "rep_c")
        write_reports(os.path.join(r, "adaptor"), os.path.join(r, "gx"), self.L,
                      extra_gx=f"u_mt\t1\t{self.L['u_mt']}\t{self.L['u_mt']}\tEXCLUDE\tprok:x\t100\tX\n")
        fcs = os.path.join(self.d, "fcs_c")
        self.assertEqual(run("fcs-summarize", "--fasta", self.fa, "--label", "S", "--agp", self.agp, "--adaptor-dir",
                             os.path.join(r, "adaptor"), "--gx-dir", os.path.join(r, "gx"), "--out", fcs).returncode, 0)
        out = os.path.join(self.d, "plan_c")
        self.assertEqual(run("plan", "--fcs", fcs, "--organelle", self.out, "--out", out).returncode, 0)
        plan = {x["record"] for x in tsv(os.path.join(out, "EXCLUSION_PLAN.candidate.tsv"))}
        self.assertNotIn("u_mt", plan)
        rev = [x for x in tsv(os.path.join(out, "REVIEW.tsv")) if x["record"] == "u_mt"]
        self.assertEqual(rev[0]["reason"], "conflicting whole-record calls")

    def test_plan_refuses_unvalidated_input(self):
        p = run("plan", "--organelle", self.d, "--out", os.path.join(self.d, "plan_x"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("VALIDATED.json", p.stderr)


if __name__ == "__main__":
    unittest.main()
