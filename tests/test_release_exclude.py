"""Tests for g2f_release_exclude + g2f_release_exclude_verify (synthetic sequences only; needs samtools and bgzip)."""
import csv
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import util
import telomere_data as TD

EXC = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_release_exclude.py")
VER = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_release_exclude_verify.py")


def run(script, *args):
    env = dict(os.environ, G2F_SAMTOOLS=util.SAMTOOLS, G2F_BGZIP=TD.BGZIP)
    return subprocess.run([sys.executable, script] + [str(x) for x in args], capture_output=True, text=True, env=env)


def wtext(p, t):
    with open(p, "w") as fh:
        fh.write(t)


def write_fa(path, recs, width):
    raw = path[:-3]
    with open(raw, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            if width:
                for i in range(0, len(s), width):
                    fh.write(s[i:i + width] + "\n")
            else:
                fh.write(s + "\n")
    subprocess.run([TD.BGZIP, "-f", raw], check=True)
    subprocess.run([util.SAMTOOLS, "faidx", path], check=True)


def read_fa(path):
    out, name = {}, None
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                name = l[1:].strip()
                out[name] = []
            else:
                out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def build_unit(d):
    rnd = random.Random(7)
    seq = lambda n: "".join(rnd.choice("ACGTacgt") for _ in range(n))
    C = {"c1": seq(500), "c2": seq(430), "u1": seq(210), "bact": seq(333), "org": seq(170)}
    order = ["c1", "c2", "u1", "bact", "org"]
    write_fa(f"{d}/U.contigs.fa.gz", [(k, C[k]) for k in order], 0)
    chrA = C["c1"] + "N" * 100 + C["c2"]
    scaf = [("U_chrA", chrA), ("u1", C["u1"]), ("bact", C["bact"]), ("org", C["org"])]
    write_fa(f"{d}/U.scaffolds.fa.gz", scaf, 80)
    agp = ["##agp-version\t2.1\n",
           "U_chrA\t1\t500\t1\tW\tc1\t1\t500\t+\n", "U_chrA\t501\t600\t2\tU\t100\tscaffold\tyes\talign_genus\n",
           "U_chrA\t601\t1030\t3\tW\tc2\t1\t430\t+\n"]
    agp += [f"{k}\t1\t{len(C[k])}\t1\tW\t{k}\t1\t{len(C[k])}\t+\n" for k in ("u1", "bact", "org")]
    wtext(f"{d}/U.agp", "".join(agp))
    wtext(f"{d}/U.gaps.bed", "U_chrA\t500\t600\n")
    ch = [f"chain 500 c1 500 + 0 500 U_chrA 1030 + 0 500 1\n500\n\n", f"chain 430 c2 430 + 0 430 U_chrA 1030 + 600 1030 2\n430\n\n"]
    ch += [f"chain {len(C[k])} {k} {len(C[k])} + 0 {len(C[k])} {k} {len(C[k])} + 0 {len(C[k])} {i}\n{len(C[k])}\n\n"
           for i, k in enumerate(("u1", "bact", "org"), 3)]
    wtext(f"{d}/U.chain", "".join(ch))
    wtext(f"{d}/U.map.tsv", "contig\tscaffold\n" + "".join(f"{k}\t{'U_chrA' if k[0] == 'c' else k}\n" for k in order))
    return C


def plan_file(path, rows):
    wtext(path, "unit\trecord\treason_class\tevidence\n" + "".join(f"U\t{r}\t{c}\tsynthetic\n" for r, c in rows))


def apply_args(d, plan, out):
    return ["apply", "--plan", plan, "--unit", "U", "--contigs", f"{d}/U.contigs.fa.gz", "--scaffolds",
            f"{d}/U.scaffolds.fa.gz", "--agp", f"{d}/U.agp", "--chain", f"{d}/U.chain", "--map", f"{d}/U.map.tsv",
            "--gaps-bed", f"{d}/U.gaps.bed", "--out", out, "--contigs-name", "U.contigs.v2.fa.gz",
            "--scaffolds-name", "U.scaffolds.v2.fa.gz", "--prefix", "U", "--threads", "1"]


def verify_args(d, plan, out, vjson):
    return ["--plan", plan, "--unit", "U", "--build", out, "--src-contigs", f"{d}/U.contigs.fa.gz",
            "--src-scaffolds", f"{d}/U.scaffolds.fa.gz", "--src-agp", f"{d}/U.agp", "--src-chain", f"{d}/U.chain",
            "--src-map", f"{d}/U.map.tsv", "--src-gaps-bed", f"{d}/U.gaps.bed", "--out", vjson]


@unittest.skipUnless(TD.HAVE_TOOLS, "samtools and bgzip needed")
class Exclude(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.C = build_unit(cls.d)
        cls.plan = f"{cls.d}/plan.tsv"
        plan_file(cls.plan, [("bact", "FOREIGN_CONTAMINANT"), ("org", "ORGANELLE_PLASTID")])
        cls.out = f"{cls.d}/build"
        r = run(EXC, *apply_args(cls.d, cls.plan, cls.out))
        assert r.returncode == 0, r.stderr
        cls.vjson = f"{cls.d}/verify.json"
        cls.v = run(VER, *verify_args(cls.d, cls.plan, cls.out, cls.vjson))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def test_verifier_passes(self):
        self.assertEqual(self.v.returncode, 0, self.v.stdout + self.v.stderr)
        self.assertEqual(json.load(open(self.vjson))["result"], "PASS")

    def test_outputs(self):
        c = read_fa(f"{self.out}/contigs/U.contigs.v2.fa.gz")
        self.assertEqual(list(c), ["c1", "c2", "u1"])
        self.assertEqual(c["u1"], self.C["u1"])
        s = read_fa(f"{self.out}/scaffolds_rg/U.scaffolds.v2.fa.gz")
        self.assertEqual(list(s), ["U_chrA", "u1"])
        agp = open(f"{self.out}/provenance/U.agp").read()
        self.assertNotIn("bact", agp)
        self.assertIn("U_chrA\t601", agp)
        self.assertEqual(open(f"{self.out}/provenance/U.gaps.bed").read(), "U_chrA\t500\t600\n")
        self.assertEqual(open(f"{self.out}/provenance/U.chain").read().count("chain "), 3)

    def test_kept_records_byte_identical(self):
        src = gzip.open(f"{self.d}/U.scaffolds.fa.gz", "rt").read()
        out = gzip.open(f"{self.out}/scaffolds_rg/U.scaffolds.v2.fa.gz", "rt").read()
        self.assertTrue(src.startswith(out.split(">u1")[0]))
        self.assertIn(">u1\n" + self.C["u1"][:80] + "\n", out)

    def test_excluded_table(self):
        rows = list(csv.DictReader(open(f"{self.out}/U.EXCLUDED_RECORDS.tsv"), delimiter="\t"))
        self.assertEqual([r["record"] for r in rows], ["bact", "org"])
        self.assertEqual(rows[0]["seq_sha256_upper"], hashlib.sha256(self.C["bact"].upper().encode()).hexdigest())
        self.assertEqual(int(rows[1]["length"]), 170)
        self.assertEqual(rows[0]["reason_class"], "FOREIGN_CONTAMINANT")

    def test_refuses_placed_contig(self):
        p = f"{self.d}/bad.tsv"
        plan_file(p, [("c2", "FOREIGN_CONTAMINANT")])
        r = run(EXC, "check-plan", "--plan", p, "--agp", f"{self.d}/U.agp", "--contigs", f"{self.d}/U.contigs.fa.gz")
        self.assertEqual(r.returncode, 1)
        self.assertIn("not an unplaced singleton", r.stdout)
        r = run(EXC, *apply_args(self.d, p, f"{self.d}/bad_build"))
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(os.path.exists(f"{self.d}/bad_build"))

    def test_refuses_unknown_record_and_overwrite(self):
        p = f"{self.d}/unk.tsv"
        plan_file(p, [("nope", "X")])
        r = run(EXC, "check-plan", "--plan", p, "--agp", f"{self.d}/U.agp", "--contigs", f"{self.d}/U.contigs.fa.gz")
        self.assertEqual(r.returncode, 1)
        r = run(EXC, *apply_args(self.d, self.plan, self.out))
        self.assertIn("refusing to overwrite", r.stderr)

    def test_verifier_catches_tampering(self):
        t = f"{self.d}/tampered"
        shutil.copytree(self.out, t)
        fa = f"{t}/contigs/U.contigs.v2.fa.gz"
        os.remove(fa)
        os.remove(fa + ".fai")
        write_fa(fa, [("c1", self.C["c1"]), ("c2", self.C["c2"])], 0)   # u1 silently lost
        vj = f"{self.d}/tampered.json"
        r = run(VER, *verify_args(self.d, self.plan, t, vj))
        self.assertEqual(r.returncode, 1)
        res = json.load(open(vj))
        self.assertEqual(res["checks"]["X1_contigs"], "FAIL")
        self.assertEqual(res["checks"]["X7_index"], "FAIL")   # BUILD.json hash no longer matches


HEADER = ("##agp-version\t2.1\n# REPRESENTATION: scaffolds v2 ({excluded} records excluded)\n"
          "# COMPONENTS: {contigs_name} (uncompressed sha256 {contigs_sha256}); objects in {scaffolds_name}\n")


@unittest.skipUnless(TD.HAVE_TOOLS, "samtools and bgzip needed")
class ExcludeAgpMetadata(unittest.TestCase):
    """--agp-header and --relabel-linkage: declared AGP metadata changes, and nothing else."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.C = build_unit(cls.d)
        agp = open(f"{cls.d}/U.agp").read().replace("\tyes\talign_genus\n", "\tyes\tmap\n")
        wtext(f"{cls.d}/U.agp", "# old comment: v1\n".join(agp.split("\n", 1)[0:1]) + "\n# old comment: v1\n" + agp.split("\n", 1)[1])
        cls.plan = f"{cls.d}/plan.tsv"
        plan_file(cls.plan, [("org", "ORGANELLE_PLASTID")])
        cls.hdr = f"{cls.d}/header.txt"
        wtext(cls.hdr, HEADER)
        cls.out = f"{cls.d}/build"
        cls.meta = ["--agp-header", cls.hdr, "--relabel-linkage", "map=align_genus"]
        r = run(EXC, *apply_args(cls.d, cls.plan, cls.out), *cls.meta)
        assert r.returncode == 0, r.stderr

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def verify(self, build, *extra):
        vj = f"{self.d}/v_{len(os.listdir(self.d))}.json"
        r = run(VER, *verify_args(self.d, self.plan, build, vj), *extra)
        return r.returncode, json.load(open(vj))

    def test_header_and_relabel(self):
        agp = open(f"{self.out}/provenance/U.agp").read()
        dec = hashlib.sha256(gzip.open(f"{self.out}/contigs/U.contigs.v2.fa.gz", "rb").read()).hexdigest()
        self.assertIn(f"(uncompressed sha256 {dec})", agp)
        self.assertIn("1 records excluded", agp)
        self.assertNotIn("old comment", agp)
        self.assertIn("\tU\t100\tscaffold\tyes\talign_genus\n", agp)
        self.assertNotIn("\tmap\n", agp)
        b = json.load(open(f"{self.out}/BUILD.json"))
        self.assertEqual(b["linkage_relabel"], {"from": "map", "to": "align_genus", "gap_lines": 1})
        self.assertEqual(b["contigs_sha256_uncompressed"], dec)
        code, res = self.verify(self.out, *self.meta)
        self.assertEqual(code, 0, res["failures"])

    def test_undeclared_change_fails(self):
        code, res = self.verify(self.out)                       # verifier not told about header/relabel
        self.assertEqual(code, 1)
        self.assertEqual(res["checks"]["X4_agp"], "FAIL")

    def test_other_gap_column_change_fails(self):
        t = f"{self.d}/tampered"
        shutil.copytree(self.out, t)
        p = f"{t}/provenance/U.agp"
        wtext(p, open(p).read().replace("\tscaffold\tyes\talign_genus", "\tscaffold\tno\talign_genus"))
        code, res = self.verify(t, *self.meta)
        self.assertEqual(code, 1)
        self.assertEqual(res["checks"]["X4_agp"], "FAIL")

    def test_builder_refusals(self):
        r = run(EXC, *apply_args(self.d, self.plan, f"{self.d}/b2"), "--relabel-linkage", "map=align_species")
        self.assertNotEqual(r.returncode, 0)
        bad = f"{self.d}/bad_header.txt"
        wtext(bad, "# no version line\n")
        r = run(EXC, *apply_args(self.d, self.plan, f"{self.d}/b3"), "--agp-header", bad)
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(os.path.exists(f"{self.d}/b2") or os.path.exists(f"{self.d}/b3"))


if __name__ == "__main__":
    unittest.main()
