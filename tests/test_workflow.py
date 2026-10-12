"""Workflow-level tests on tiny synthetic data.

* Dry-run DAG resolution (no tools executed).
* Complete end-to-end executions in which seqkit, samtools, md5sum and awk
  are the real tools, while hifiasm, pbtk, meryl, Merqury and BUSCO are
  replaced by stubs (tests/stubs/). These test plumbing, validation gates,
  provenance, rerun behaviour and resumability; they do NOT validate the
  scientific tools (see tests/test_realtools.py for that).
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml
from util import (HAVE_SAMTOOLS, HAVE_SEQKIT, PIPELINE_DIR, SAMTOOLS, SEQKIT, hifi_names,
                  mkdirs, ont_names, write_bam, write_bytes, write_fasta, write_fastq)

STUBS = os.path.join(PIPELINE_DIR, "tests", "stubs")
HIFI_ONT = "hifiasm-0.25.0_hifi-ont-ul50k_l0_v1"
HIFI_ONLY = "hifiasm-0.25.0_hifi_l0_v1"
M1, M2, M3 = "m84109_240806_024441_s4", "m84109_240611_191050_s2", "m64001_220101_000000"
HEAVY = {"hifiasm", "checksum", "hifi_read_stats", "hifi_bam_to_fastq", "ont_input_stats",
         "ont_filter", "meryl_count", "busco", "merqury", "gfa_to_fasta", "assembly_stats"}


def have_snakemake():
    r = subprocess.run([sys.executable, "-m", "snakemake", "--version"], capture_output=True)
    return r.returncode == 0


def jobs(out):
    """Job names from a Snakemake dry-run job-stats table."""
    names, on = set(), False
    for line in out.splitlines():
        if line.startswith("job ") and "count" in line:
            on = True
            continue
        if on:
            if line.startswith("total"):
                break
            parts = line.split()
            if parts and not set(parts[0]) <= {"-"}:
                names.add(parts[0])
    return names


def unprotect(root):
    for p in glob.glob(f"{root}/**", recursive=True):
        if os.path.isfile(p) and not os.path.islink(p):
            os.chmod(p, 0o644)


def stub_tools():
    return {"containers_dir": ".", "tools": {
        "hifiasm": {"cmd": f"{STUBS}/hifiasm", "version_cmd": f"{STUBS}/hifiasm --version"},
        "seqkit": {"cmd": SEQKIT, "version_cmd": f"{SEQKIT} version"},
        "samtools": {"cmd": SAMTOOLS, "version_cmd": f"{SAMTOOLS} --version | head -1"},
        "bam2fastq": {"cmd": f"env SAMTOOLS={SAMTOOLS} {STUBS}/bam2fastq",
                      "version_cmd": f"{STUBS}/bam2fastq --version"},
        "pbindex": {"cmd": f"{STUBS}/pbindex", "version_cmd": f"{STUBS}/pbindex --version"},
        "meryl": {"cmd": f"{STUBS}/meryl", "version_cmd": f"{STUBS}/meryl --version"},
        "merqury": {"cmd": f"{STUBS}/merqury.sh", "version_cmd": "echo stub"},
        "busco": {"cmd": f"{STUBS}/busco", "version_cmd": f"{STUBS}/busco --version"},
        "md5sum": {"cmd": "md5sum", "version_cmd": "md5sum --version | head -1"},
    }}


class WorkflowCase(unittest.TestCase):
    """Temporary project with config files that tests may modify."""
    tools = staticmethod(stub_tools)
    lineage = ("poales_odb12.2", b"name=poales_odb12.2\ndomain=eukaryota\nOrthoDB_version=12.2\n"
                                 b"creation_date=2026-05-13\n")

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="g2f_wf_")
        t = self.tmp
        self.build_data()
        name, cfgtext = self.lineage
        if cfgtext is not None:
            os.makedirs(os.path.join(t, "lineages", name))
            write_bytes(os.path.join(t, "lineages", name, "dataset.cfg"), cfgtext)
        self.samples = {"samples": {"S1": {"aliases": ["EB2_S1"]}}}
        with open(os.path.join(PIPELINE_DIR, "config", "science.yaml")) as fh:
            self.science = yaml.safe_load(fh)
        self.adjust_science()
        self.write_yaml("tools.yaml", self.tools())
        self.write_yaml("samples.yaml", self.samples)
        self.write_yaml("science.yaml", self.science)
        self.write_yaml("resources.yaml", {"tmpdir": f"{t}/tmp", "rules": {"default": {
            "threads": 1, "mem_mb": 1000}}})
        self.write_yaml("local.yaml", {"busco_lineage_dir": self.lineage_dir()})
        self.write_yaml("external.yaml", {})
        self.config = {
            "data_dir": f"{t}/data", "results_dir": f"{t}/results", "legacy_dirs": [f"{t}/legacy"],
            "samples": f"{t}/samples.yaml", "tools": f"{t}/tools.yaml", "local": f"{t}/local.yaml",
            "science": f"{t}/science.yaml", "resources": f"{t}/resources.yaml",
            "external_assemblies": f"{t}/external.yaml", "tmpdir": f"{t}/tmp"}
        self.write_yaml("config.yaml", self.config)

    def build_data(self):
        t = self.tmp
        s1 = mkdirs(t, "S1", ("hifi", "ont"))
        write_bam(f"{s1['hifi']}/{M1}.hifi_reads.bc2027.bam", hifi_names(M1, 10), [3000] * 10, movie=M1)
        write_fastq(f"{s1['hifi']}/{M2}.hifi.fastq.gz", hifi_names(M2, 10), [3500] * 10, seed=2)
        write_fastq(f"{s1['ont']}/fc1.fastq.gz", ont_names(12, seed=1), [4000, 52000, 61000] * 4, seed=3)
        write_fastq(f"{s1['ont']}/fc2.fastq.gz", ont_names(6, seed=2), [70000, 1000] * 3, seed=4)
        s2 = mkdirs(t, "S2", ("hifi",))
        write_fasta(f"{s2['hifi']}/S2.hifi.fasta.gz", hifi_names(M3, 12), [4000] * 12, seed=5)

    def adjust_science(self):
        # The stub BUSCO reports a real odb10 result (3236 markers), so the
        # stub tests pin only the fields their fake dataset.cfg provides.
        self.science["qc"]["busco"]["lineages"] = {"poales_odb12.2": {
            "datasets_version": "odb12.2", "creation_date": "2026-05-13"}}

    def lineage_dir(self):
        return os.path.join(self.tmp, "lineages")

    def tearDown(self):
        unprotect(self.tmp)
        shutil.rmtree(self.tmp)

    def write_yaml(self, name, obj):
        with open(os.path.join(self.tmp, name), "w") as fh:
            yaml.safe_dump(obj, fh)

    def smk(self, *args, env=None, cores="2"):
        cmd = [sys.executable, "-m", "snakemake", "-s", f"{PIPELINE_DIR}/workflow/Snakefile",
               "--directory", self.tmp, "--configfile", f"{self.tmp}/config.yaml",
               "--cores", cores, *args, "--rerun-triggers", "mtime"]
        e = dict(os.environ, TMPDIR=self.tmp, **(env or {}))
        r = subprocess.run(cmd, capture_output=True, text=True, env=e)
        r.all = r.stdout + r.stderr
        return r

    def summary_rows(self, pattern="assembly_summary.*.tsv", vid=None):
        paths = sorted(glob.glob(f"{self.tmp}/results/summary/{pattern}"))
        if vid:
            paths = [p for p in paths if vid in p]
        self.assertEqual(len(paths), 1, paths)
        with open(paths[0]) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh]
        return {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}


@unittest.skipUnless(HAVE_SAMTOOLS and HAVE_SEQKIT and have_snakemake(),
                     "needs samtools, seqkit and snakemake (module load python/3.12.7)")
class WorkflowTest(WorkflowCase):
    # ------------------------------------------------------------------ dry runs
    def test_dry_run_resolves_both_methods(self):
        r = self.smk("--dry-run", "--printshellcmds")
        self.assertEqual(r.returncode, 0, r.all)
        self.assertIn(f"results/S1/assembly/{HIFI_ONT}/", r.all)
        self.assertIn(f"results/S2/assembly/{HIFI_ONLY}/", r.all)
        ul = (f"--ul {self.tmp}/results/S1/reads/ont/ul50000/fc1.fastq.gz,"
              f"{self.tmp}/results/S1/reads/ont/ul50000/fc2.fastq.gz")
        self.assertIn(ul, r.all)                           # comma-separated, no merged copy
        self.assertIn("seqkit seq -m 50000", r.all.replace(SEQKIT, "seqkit"))
        self.assertIn("hifi_bam_to_fastq", r.all)
        self.assertNotIn("results/S2/reads/ont", r.all)
        self.assertFalse(os.path.exists(f"{self.tmp}/results"), "dry run must not write results")

    def test_invalid_ont_halts_without_fallback(self):
        os.makedirs(f"{self.tmp}/data/S2/ont")
        write_bytes(f"{self.tmp}/data/S2/ont/run.pod5", b"x")
        r = self.smk("--dry-run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("RAW_SIGNAL", r.all)

    def test_prefiltered_ont_used_as_supplied(self):
        d = f"{self.tmp}/data/S2/ont"
        os.makedirs(d)
        write_fastq(f"{d}/S2_ONT_gt50kb.fastq.gz", ont_names(5, seed=9),
                    [55000, 60000, 70000, 51000, 50001])
        self.samples["samples"]["S2"] = {"acknowledge": {
            "S2_ONT_gt50kb.fastq.gz": {"prefiltered_ont": "test: already filtered"}}}
        self.write_yaml("samples.yaml", self.samples)
        r = self.smk("--dry-run", "--printshellcmds")
        self.assertEqual(r.returncode, 0, r.all)
        self.assertIn(f"results/S2/assembly/{HIFI_ONT}/", r.all)
        self.assertIn(f"--ul {d}/S2_ONT_gt50kb.fastq.gz ", r.all)     # used in place
        self.assertNotIn("results/S2/reads/ont/ul50000/S2_ONT_gt50kb.fastq.gz", r.all)
        r = self.smk()
        self.assertEqual(r.returncode, 0, r.stderr[-4000:])
        row = self.summary_rows()["S2"]
        self.assertEqual(row["ont_input_type"], "prefiltered")
        self.assertEqual(row["ont_raw_reads"], "")                  # never reported as raw
        self.assertEqual(row["ont_supplied_reads"], "5")
        self.assertEqual(row["ont_retained_reads"], "5")

    def test_ont_tar_streamed_in_place(self):
        """MinKNOW tar delivery: pass members used (never extracted), fail members excluded."""
        import gzip
        import io
        import random
        import tarfile
        d = f"{self.tmp}/data/S2/ont"
        os.makedirs(d)
        rng = random.Random(5)
        names = ont_names(7, seed=11)

        def member(recs):
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
                for n, L in recs:
                    seq = "".join(rng.choice("ACGT") for _ in range(L))
                    gz.write(f"@{n}\n{seq}\n+\n{'5' * L}\n".encode())
            return buf.getvalue()

        with tarfile.open(f"{d}/run1_PAQ00000_abcd1234_fastq.tar", "w") as tf:
            for name, data in (
                    ("./fastq_pass/p_pass_x_0.fastq.gz", member(zip(names[:3], [55000, 1000, 60000]))),
                    ("./fastq_pass/p_pass_x_1.fastq.gz", member(zip(names[3:5], [70000, 2000]))),
                    ("./fastq_fail/p_fail_x_0.fastq.gz", member(zip(names[5:], [80000, 90000])))):
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
        r = self.smk("--dry-run", "--printshellcmds")
        self.assertEqual(r.returncode, 0, r.all)
        self.assertIn(f"results/S2/assembly/{HIFI_ONT}/", r.all)
        self.assertIn("tar-fastq", r.all)
        r = self.smk()
        self.assertEqual(r.returncode, 0, r.stderr[-4000:])
        row = self.summary_rows()["S2"]
        self.assertEqual(row["ont_input_type"], "raw")
        self.assertEqual(row["ont_raw_reads"], "5")                 # pass members only
        self.assertEqual(row["ont_retained_reads"], "3")            # >= 50 kb
        extracted = glob.glob(f"{self.tmp}/results/**/p_pass_x_*", recursive=True)
        self.assertEqual(extracted, [])

    # -------------------------------------------------------------- end to end
    def test_end_to_end_and_rerun_behaviour(self):
        r = self.smk()
        self.assertEqual(r.returncode, 0, r.stdout[-4000:] + r.stderr[-6000:])
        res = f"{self.tmp}/results"
        recs = self.summary_rows()
        self.assertEqual(set(recs), {"S1", "S2"})
        s1, s2 = recs["S1"], recs["S2"]
        self.assertEqual((s1["method_id"], s2["method_id"]), (HIFI_ONT, HIFI_ONLY))
        self.assertEqual(s1["assembly_source"], "pipeline")
        self.assertEqual(s1["hifi_reads"], "20")
        self.assertEqual(s1["hifi_bases"], str(10 * 3000 + 10 * 3500))
        self.assertEqual(s1["ont_input_type"], "raw")
        self.assertEqual(s1["ont_raw_reads"], "18")
        self.assertEqual(s1["ont_retained_reads"], str(8 + 3))   # >= 50 kb only
        self.assertEqual(s2["ont_raw_reads"], "")
        self.assertEqual(s1["contigs"], "3")
        self.assertEqual(s1["busco_lineage"], "poales_odb12.2")
        self.assertEqual(s1["busco_datasets_version"], "odb12.2")
        self.assertEqual(s1["merqury_qv"], "45.2")
        self.assertEqual(s1["hifiasm_version"], "0.25.0-r726")
        self.assertEqual(s1["kmer_assessment"], "not_flagged")
        vid = s1["validation_id"]
        self.assertTrue(os.path.exists(f"{res}/.identity/validation-{vid}.json"))
        with open(f"{res}/S1/provenance/hifiasm/{HIFI_ONT}.json") as fh:
            prov = json.load(fh)
        self.assertEqual(prov["status"], "completed")
        self.assertEqual(prov["tools"]["hifiasm"]["reported"], "0.25.0-r726")
        self.assertIn("--ul", prov["command"])
        self.assertGreater(prov["resource_usage"]["max_rss_gb_largest_process"], 0)
        with open(f"{res}/S1/assembly/{HIFI_ONT}/S1.{HIFI_ONT}.method.json") as fh:
            self.assertEqual(json.load(fh)["fingerprint"], s1["method_fingerprint"])
        self.assertTrue(os.path.exists(
            f"{res}/S1/reads/hifi/converted/{M1}.hifi_reads.bc2027.pbtk/in.bam.pbi"))
        self.assertFalse(glob.glob(f"{self.tmp}/data/S1/hifi/*.pbi"), "data/ must not be modified")
        fa = f"{res}/S1/assembly/{HIFI_ONT}/S1.{HIFI_ONT}.p_ctg.fa"
        self.assertFalse(os.access(fa, os.W_OK))

        # 1. Resumable: nothing to do.
        self.assertIn("Nothing to be done", self.smk("--dry-run").all)

        # 2. Resource-only change: still nothing to do.
        self.write_yaml("resources.yaml", {"tmpdir": f"{self.tmp}/tmp", "rules": {
            "default": {"threads": 2, "mem_mb": 2000}, "hifiasm_hifi": {"threads": 8}}})
        self.assertIn("Nothing to be done", self.smk("--dry-run").all)

        # 3. Science change: only cheap validation/summary jobs rerun, under a new identity.
        self.science["genome_size_bp"] = 2_000_000_000
        self.write_yaml("science.yaml", self.science)
        r = self.smk("--dry-run")
        self.assertEqual(r.returncode, 0, r.all)
        ran = jobs(r.all)
        self.assertIn("validate_hifi_reads", ran)
        self.assertIn("sample_summary", ran)
        self.assertFalse(ran & HEAVY, f"heavy jobs would rerun: {ran & HEAVY}")
        r = self.smk()
        self.assertEqual(r.returncode, 0, r.stderr[-4000:])
        tables = glob.glob(f"{res}/summary/assembly_summary.*.tsv")
        self.assertEqual(len(tables), 2)                  # old identity kept, new one added
        new_vid = [os.path.basename(p).split(".")[1] for p in tables if vid not in p][0]
        self.assertTrue(glob.glob(f"{res}/.identity/validation-{new_vid}.json"))
        new = self.summary_rows(vid=new_vid)
        self.assertEqual(new["S1"]["hifi_approx_cov"], str(round(65000 / 2e9, 2)))
        self.assertEqual(new["S1"]["validation_id"], new_vid)
        self.assertTrue(glob.glob(f"{res}/S1/summary/S1.{HIFI_ONT}.validated.{vid}.*.json"))

        # 4. Editing a used method in place is refused.
        mdir = f"{self.tmp}/methods"
        shutil.copytree(f"{PIPELINE_DIR}/config/methods", mdir)
        mf = f"{mdir}/{HIFI_ONLY}.yaml"
        with open(mf) as fh:
            m = yaml.safe_load(fh)
        m["args"] = ["-l1"]
        with open(mf, "w") as fh:
            yaml.safe_dump(m, fh)
        self.write_yaml("config.yaml", dict(self.config, methods_dir=mdir))
        r = self.smk("--dry-run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("was changed after it was used", r.all)
        self.write_yaml("config.yaml", self.config)

        # 5. Changing inputs under existing results is refused.
        write_fastq(f"{self.tmp}/data/S2/hifi/extra.fastq.gz",
                    hifi_names("m99999_250101_000000", 8), [3000] * 8, seed=9)
        r = self.smk("--dry-run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("differ from the inventory", r.all)

    def test_failed_assembly_resumes_safely(self):
        target = f"{self.tmp}/results/S2/assembly/{HIFI_ONLY}/S2.asm.bp.p_ctg.gfa"
        r = self.smk(target, env={"STUB_HIFIASM_FAIL": "1"})
        self.assertNotEqual(r.returncode, 0)
        d = os.path.dirname(target)
        self.assertFalse(os.path.exists(target))
        self.assertFalse(os.path.exists(f"{d}/S2.{HIFI_ONLY}.method.json"))
        self.assertTrue(os.path.exists(f"{d}/S2.asm.ec.bin"))           # stale intermediate
        with open(f"{self.tmp}/results/S2/provenance/hifiasm/{HIFI_ONLY}.json") as fh:
            prov = json.load(fh)
        self.assertEqual((prov["status"], prov["failure"]), ("failed", "nonzero_exit"))
        r = self.smk(target)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Refusing to start: intermediate hifiasm .bin files", r.all)
        os.rename(f"{d}/S2.asm.ec.bin", f"{self.tmp}/moved_aside.ec.bin")    # human decision
        r = self.smk(target)
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])

    def test_cleanup_switches(self):
        """cleanup: converted HiFi FASTQ is temporary; hifiasm .bin files are removed only
        after a successful run (a failed run keeps them for diagnosis)."""
        self.write_yaml("resources.yaml", {"tmpdir": f"{self.tmp}/tmp", "rules": {"default": {
            "threads": 1, "mem_mb": 1000}}, "cleanup": {
            "temp_converted_hifi": True, "remove_hifiasm_bins_after_success": True}})
        d1 = f"{self.tmp}/results/S1/assembly/{HIFI_ONT}"
        d2 = f"{self.tmp}/results/S2/assembly/{HIFI_ONLY}"
        r = self.smk(f"{d2}/S2.asm.bp.p_ctg.gfa", env={"STUB_HIFIASM_FAIL": "1"})
        self.assertNotEqual(r.returncode, 0)
        self.assertTrue(os.path.exists(f"{d2}/S2.asm.ec.bin"))          # kept after failure
        os.remove(f"{d2}/S2.asm.ec.bin")
        r = self.smk("per_sample")
        self.assertEqual(r.returncode, 0, r.stderr[-4000:])
        self.assertFalse(glob.glob(f"{d1}/*.bin"))
        with open(f"{d1}/S1.asm.removed_bin_files.txt") as fh:
            self.assertIn("S1.asm.ec.bin", fh.read())
        self.assertFalse(glob.glob(f"{self.tmp}/results/S1/reads/hifi/converted/*.fastq.gz"))
        self.assertTrue(glob.glob(f"{self.tmp}/results/*/summary/*.validated.*.json"))
        self.assertFalse(glob.glob(f"{self.tmp}/results/summary/assembly_summary.*.tsv"))
        r = self.smk("--dry-run", "per_sample")                        # nothing to redo
        self.assertEqual(r.returncode, 0, r.all)
        self.assertFalse(jobs(r.all) & HEAVY, r.all[-3000:])

    def test_concurrent_per_sample_processes(self):
        """The launcher runs one Snakemake process per sample at the same time in the
        same directory; their DAGs must not conflict (shared identity stamp made first)."""
        r = self.smk("identity_stamp")
        self.assertEqual(r.returncode, 0, r.all)
        procs = []
        for s in ("S1", "S2"):
            cmd = [sys.executable, "-m", "snakemake", "-s", f"{PIPELINE_DIR}/workflow/Snakefile",
                   "--directory", self.tmp, "--configfile", f"{self.tmp}/config.yaml",
                   "--cores", "1", "per_sample", "--rerun-triggers", "mtime",
                   "--config", f"only_samples={s}"]
            procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                          text=True, env=dict(os.environ, TMPDIR=self.tmp)))
        outs = [p.communicate()[0] for p in procs]
        for p, out in zip(procs, outs):
            self.assertEqual(p.returncode, 0, out[-4000:])
        self.assertEqual(len(glob.glob(f"{self.tmp}/results/*/summary/*.validated.*.json")), 2)
        r = self.smk()                                                 # project table only
        self.assertEqual(r.returncode, 0, r.all[-3000:])
        self.assertEqual(set(self.summary_rows()), {"S1", "S2"})

    def test_kmer_flag_warns_by_default_and_review_is_recorded(self):
        r = self.smk(env={"STUB_KMER_HET": "1"})
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        self.assertIn("k-mer screen flagged", r.all)
        recs = self.summary_rows()
        self.assertEqual(recs["S2"]["kmer_assessment"], "flagged_unreviewed")
        self.assertEqual(recs["S2"]["kmer_screen_flagged"], "True")
        self.samples["samples"]["S2"] = {"kmer_review": {
            "status": "accepted", "reason": "test", "reviewer": "t", "date": "2026-10-11"}}
        self.write_yaml("samples.yaml", self.samples)
        r = self.smk("--dry-run", env={"STUB_KMER_HET": "1"})
        ran = jobs(r.all)
        self.assertIn("kmer_assess", ran, r.all[-3000:])
        self.assertFalse(ran & HEAVY, ran)
        r = self.smk(env={"STUB_KMER_HET": "1"})
        self.assertEqual(r.returncode, 0, r.stderr[-3000:])
        tables = sorted(glob.glob(f"{self.tmp}/results/summary/assembly_summary.*.tsv"),
                        key=os.path.getmtime)
        self.assertEqual(len(tables), 2)
        with open(tables[-1]) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh]
        cur = {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}
        self.assertEqual(cur["S2"]["kmer_assessment"], "flagged_reviewed")
        self.assertEqual(cur["S1"]["kmer_assessment"], "flagged_unreviewed")

    def test_kmer_blocking_mode(self):
        self.science["kmer"]["screen"]["block_assembly_on_flag"] = True
        self.write_yaml("science.yaml", self.science)
        target = f"{self.tmp}/results/S2/assembly/{HIFI_ONLY}/S2.asm.bp.p_ctg.gfa"
        r = self.smk(target, env={"STUB_KMER_HET": "1"})
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Assembly blocked pending human review", r.all)
        self.assertFalse(os.path.exists(target))

    def test_wrong_busco_dataset_is_rejected(self):
        # The installed dataset.cfg does not match the pinned odb12.2 identity.
        write_bytes(f"{self.tmp}/lineages/poales_odb12.2/dataset.cfg",
                    b"name=poales_odb12.2\ncreation_date=2026-05-22\n")
        r = self.smk(f"{self.tmp}/results/S2/qc/{HIFI_ONLY}/busco/poales_odb12.2.validated.json")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("pinned dataset requires", r.all)

    # ------------------------------------------------------- imported assemblies
    def test_imported_legacy_assembly_qc_is_read_only(self):
        legacy = f"{self.tmp}/legacy/run1"
        os.makedirs(legacy)
        seqs = {"ptg000001l": "ACGT" * 500, "ptg000002l": "GGCCAATT" * 100}
        with open(f"{legacy}/X.bp.p_ctg.gfa", "w") as fh:
            for n, s in seqs.items():
                fh.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\n")
        with open(f"{legacy}/X.bp.p_ctg.fa", "w") as fh:
            for n, s in seqs.items():
                fh.write(f">{n}\n{s}\n")
        write_bytes(f"{legacy}/X.stderr.log", b"[M::main] Version: 0.25.0-r726\n")
        before = {p: (os.path.getsize(p), os.path.getmtime(p)) for p in glob.glob(f"{legacy}/*")}
        self.write_yaml("external.yaml", {"root": f"{self.tmp}/legacy", "assemblies": {"S2": {
            "run1": {"p_ctg_gfa": "run1/X.bp.p_ctg.gfa", "p_ctg_fa": "run1/X.bp.p_ctg.fa",
                     "hifiasm_log": "run1/X.stderr.log", "qc": ["stats", "busco"]}}}})
        r = self.smk("imported_qc")
        self.assertEqual(r.returncode, 0, r.stderr[-4000:])
        self.assertNotIn("hifi_read_stats", r.all)                    # no read processing
        row = self.summary_rows("imported_assemblies_summary.*.tsv")["S2"]
        self.assertEqual(row["method_id"], "imported-run1")
        self.assertEqual(row["assembly_source"], "imported_legacy")
        self.assertEqual(row["contigs"], "2")
        self.assertEqual(row["hifiasm_version"], "0.25.0-r726")
        self.assertEqual(row["hifi_reads"], "")
        self.assertFalse(glob.glob(f"{self.tmp}/results/summary/assembly_summary.*.tsv"))
        after = {p: (os.path.getsize(p), os.path.getmtime(p)) for p in glob.glob(f"{legacy}/*")}
        self.assertEqual(before, after, "legacy directory must be untouched")
        self.assertEqual(sorted(os.listdir(legacy)), ["X.bp.p_ctg.fa", "X.bp.p_ctg.gfa", "X.stderr.log"])
        # A registered historical file that changes afterwards is refused.
        with open(f"{legacy}/X.bp.p_ctg.fa", "a") as fh:
            fh.write(">extra\nACGT\n")
        r = self.smk("--dry-run", "imported_qc")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("changed or moved since import", r.all)


if __name__ == "__main__":
    unittest.main()
