"""Tests for g2f_agp_revise + g2f_agp_revise_verify (synthetic sequences only; needs samtools and bgzip)."""
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
import g2f_agp_revise as R


def wtext(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rows_text(path):
    with open(path) as fh:
        return fh.readlines()

REV = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_agp_revise.py")
VER = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_agp_revise_verify.py")
TEL = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_telomere.py")
MINIMAP2 = os.environ.get("G2F_MINIMAP2") or shutil.which("minimap2")


def run(script, *args):
    env = dict(os.environ, G2F_SAMTOOLS=util.SAMTOOLS, G2F_BGZIP=TD.BGZIP)
    env["PATH"] = os.path.dirname(util.SAMTOOLS) + os.pathsep + os.path.dirname(TD.BGZIP) + os.pathsep + env["PATH"]
    return subprocess.run([sys.executable, script] + [str(x) for x in args], capture_output=True, text=True, env=env)


def fasta(path):
    out, name = {}, None
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                name = l[1:].strip()
                out[name] = []
            else:
                out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


@unittest.skipUnless(TD.HAVE_TOOLS, "samtools and bgzip needed")
class Revise(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.p = TD.build(cls.d, gz=True)
        cls.out = os.path.join(cls.d, "rev")
        r = run(REV, "apply", "--plan", cls.p["plan"], "--agp", cls.p["old_agp"], "--contigs", cls.p["contigs"],
                "--out", cls.out, "--prefix", "UNIT.v2", "--threads", "1")
        assert r.returncode == 0, r.stderr
        cls.vjson = os.path.join(cls.d, "verify.json")
        cls.v = run(VER, "--old-agp", cls.p["old_agp"], "--new-agp", f"{cls.out}/UNIT.v2.agp", "--plan", cls.p["plan"],
                    "--contigs", cls.p["contigs"], "--old-scaffolds", cls.p["old_scaffolds"],
                    "--new-scaffolds", f"{cls.out}/UNIT.v2.fa.gz", "--chain", f"{cls.out}/UNIT.v2.chain",
                    "--gaps-bed", f"{cls.out}/UNIT.v2.gaps.bed", "--map", f"{cls.out}/UNIT.v2.contig_to_scaffold.tsv",
                    "--out", cls.vjson, "--intervals", "60")

    def test_output_matches_expected_layout(self):
        self.assertEqual(fasta(f"{self.out}/UNIT.v2.fa.gz"), self.p["new_expected_seqs"])
        body = [l for l in rows_text(f"{self.out}/UNIT.v2.agp") if not l.startswith("#")]
        exp = [l for l in rows_text(self.p["new_expected_agp"]) if not l.startswith("#")]
        self.assertEqual(body, exp)                  # unplaced ctgX inserted in ID order after ctgU

    def test_build_record_and_manifest(self):
        b = jload(f"{self.out}/BUILD.json")
        self.assertEqual(b["state"], "built (not yet verified)")
        self.assertEqual(b["agp_line_changes"], 4)    # 1 flip + gap removed + ctgX removed + ctgX new object
        m = rows(f"{self.out}/UNIT.v2.change_manifest.tsv")
        self.assertEqual([r["change"] for r in m], ["orientation", "removed from object", "removed from object", "new unplaced object"])
        c2s = {r["contig"]: r for r in rows(f"{self.out}/UNIT.v2.contig_to_scaffold.tsv")}
        self.assertEqual(c2s["ctgX"]["status"], "UNPLACED_BEYOND_TELOMERE")
        self.assertEqual(c2s["ctgT"]["orientation"], "-")

    def test_independent_verification_passes(self):
        self.assertEqual(self.v.returncode, 0, self.v.stdout + self.v.stderr)
        res = jload(self.vjson)
        self.assertTrue(res["pass"])
        self.assertEqual(set(res["checks"]), {"agp_syntax", "agp_reconstructs", "diff_equals_plan", "sequence", "accounting",
                                              "derived_tables", "lift_roundtrip"})

    def test_verifier_catches_unplanned_change(self):
        plan2 = os.path.join(self.d, "plan_missing_row.tsv")
        with open(self.p["plan"]) as fi, open(plan2, "w") as fo:
            fo.writelines(l for l in fi if "ctgT" not in l)
        r = run(VER, "--old-agp", self.p["old_agp"], "--new-agp", f"{self.out}/UNIT.v2.agp", "--plan", plan2,
                "--contigs", self.p["contigs"], "--old-scaffolds", self.p["old_scaffolds"], "--new-scaffolds", f"{self.out}/UNIT.v2.fa.gz",
                "--out", os.path.join(self.d, "v_bad_plan.json"), "--intervals", "20")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("diff_equals_plan", jload(os.path.join(self.d, "v_bad_plan.json"))["fails"])

    def test_verifier_catches_tampered_fasta(self):
        seqs = fasta(f"{self.out}/UNIT.v2.fa.gz")
        s = seqs["UNIT_chr01"]
        seqs["UNIT_chr01"] = s[:5000] + ("A" if s[5000] != "A" else "C") + s[5001:]
        bad = TD.write_fasta(os.path.join(self.d, "tampered.fa.gz"), seqs, gz=True)
        r = run(VER, "--old-agp", self.p["old_agp"], "--new-agp", f"{self.out}/UNIT.v2.agp", "--plan", self.p["plan"],
                "--contigs", self.p["contigs"], "--old-scaffolds", self.p["old_scaffolds"], "--new-scaffolds", bad,
                "--out", os.path.join(self.d, "v_tampered.json"), "--intervals", "20")
        self.assertNotEqual(r.returncode, 0)
        f = jload(os.path.join(self.d, "v_tampered.json"))["fails"]
        self.assertIn("agp_reconstructs", f)
        self.assertIn("sequence", f)

    def test_check_plan_rejects_non_terminal_unplace_and_wrong_orientation(self):
        _, body = R.read_agp(self.p["old_agp"])
        bad = [dict(action="UNPLACE_TERMINAL_EXTENSION", object="UNIT_chr01", component="ctgA", old_orientation="",
                    new_orientation="", new_status="", evidence_class="", evidence=""),
               dict(action="REVERSE_ORIENTATION", object="UNIT_chr02", component="ctgD", old_orientation="-",
                    new_orientation="+", new_status="", evidence_class="", evidence="")]
        probs = R.check(body, bad)
        self.assertTrue(any("not the terminal run" in p for p in probs), probs)
        self.assertTrue(any("orientation" in p for p in probs), probs)
        self.assertEqual(R.check(body, rows(self.p["plan"])), [])

    def test_apply_refuses_overwrite(self):
        r = run(REV, "apply", "--plan", self.p["plan"], "--agp", self.p["old_agp"], "--contigs", self.p["contigs"],
                "--out", self.out, "--prefix", "UNIT.v2")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("refusing to overwrite", r.stderr)

    def test_telomere_placement_improves(self):
        """end to end: scans of the old and revised scaffolds -> compare: 2 ends fixed, nothing regressed, no internal arrays"""
        sci, tidk = self.p["science"], TD.STUB_TIDK
        for lab, fa in (("old", self.p["old_scaffolds"]), ("new", f"{self.out}/UNIT.v2.fa.gz")):
            r = run(TEL, "scan", "--fasta", fa, "--label", lab, "--out", os.path.join(self.d, f"s_{lab}"), "--science", sci, "--tidk", tidk)
            self.assertEqual(r.returncode, 0, r.stderr)
        r = run(TEL, "compare", "--before", os.path.join(self.d, "s_old"), "--after", os.path.join(self.d, "s_new"),
                "--label", "UNIT", "--out", os.path.join(self.d, "cmp"), "--science", sci)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        s = jload(os.path.join(self.d, "cmp", "summary.json"))
        self.assertEqual((s["improved"], s["regressed"], s["internal_kb_arrays_after"]), (2, 0, 0))

    @unittest.skipUnless(MINIMAP2, "minimap2 not available (set G2F_MINIMAP2)")
    def test_orient_evidence_proposes_the_flip(self):
        ref = os.path.join(self.d, "ref.fa")
        seqs = {k: v for k, v in self.p["new_expected_seqs"].items() if "_chr" in k}
        wtext(ref, "".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        out = os.path.join(self.d, "orient.tsv")
        r = run(REV, "orient-evidence", "--agp", self.p["old_agp"], "--scaffolds", self.p["old_scaffolds"], "--reference", ref,
                "--objects-regex", r"_chr0[12]$", "--end", "start", "--out", out, "--window", "20000", "--minimap2", MINIMAP2,
                "--threads", "2", "--min-mapq", "0", "--tmpdir", self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        plan = rows(out + ".proposed_plan.tsv")
        self.assertEqual([(p["action"], p["component"], p["new_orientation"]) for p in plan], [("REVERSE_ORIENTATION", "ctgT", "-")])


if __name__ == "__main__":
    unittest.main()
