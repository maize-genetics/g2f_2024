"""Tests for scripts/g2f_anchorwave.py (synthetic data only).

Unit tests always run. The end-to-end tests run the real AnchorWave container
on the synthetic genome in anchorwave_data.py (seconds, 2 threads) and are
skipped unless ANCHORWAVE_IMAGE names the .sif image.
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, HERE)

import g2f_anchorwave as aw  # noqa: E402
import anchorwave_data  # noqa: E402

IMAGE = os.environ.get("ANCHORWAVE_IMAGE")
ANCHORS = """#anchorwave proali -i ref.gff3 ...
refChr\treferenceStart\treferenceEnd\tqueryChr\tqueryStart\tqueryEnd\tstrand\tgene\tblockIndex\tscore
#block begin
chrA\t100\t200\tctg1\t100\t200\t+\tG1_T001\t1\t1
chrA\t201\t500\tctg1\t201\t500\t+\tinteranchor\t1\tNA
chrA\t501\t600\tctg1\t501\t600\t+\tG2_T001\t1\t1
#block begin
chrA\t100\t200\tctg9\t900\t1000\t-\tG1_T001\t2\t1
"""


def tmpdir():
    base = os.environ.get("TMPDIR")
    if not base or not os.path.isdir(base):
        raise unittest.SkipTest("TMPDIR must name an existing scratch directory (tests/run_tests.sh sets it)")
    return tempfile.mkdtemp(prefix="awtest_", dir=base)


def write_config(root, data, science_overrides=None):
    sci_path = os.path.join(HERE, "..", "config", "anchorwave_science.yaml")
    with open(sci_path) as fh:
        sci = yaml.safe_load(fh)
    for k, v in (science_overrides or {}).items():
        sci[k] = v
    sp = os.path.join(root, "science.yaml")
    with open(sp, "w") as fh:
        yaml.safe_dump(sci, fh)
    cfg = {"profile": "awtest", "outdir": os.path.join(root, "out"), "tmpdir": os.path.join(root, "tmp"),
           "science": sp,
           "runtime": {"image": IMAGE or "/nonexistent.sif", "bind": [root]},
           "resources": {"threads": 2, "cpus": "0-1", "nice": 10, "max_rss_gb": 16, "min_free_ram_gb": 1,
                         "reserve_disk_gb": 1, "max_cpu_busy_frac": 1.0},
           "reference": {"id": "REF", "fasta": data["ref_fa"], "gff3": data["gff"]},
           "samples": {"SAMPLE_A": {"assemblies": {"unpurged": {"fasta": data["asm_fa"]}}}}}
    p = os.path.join(root, "config.yaml")
    with open(p, "w") as fh:
        yaml.safe_dump(cfg, fh)
    return p


class TestParsers(unittest.TestCase):
    def setUp(self):
        self.d = tmpdir()

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_anchors(self):
        p = os.path.join(self.d, "anchors")
        with open(p, "w") as fh:
            fh.write(ANCHORS)
        rows = aw.parse_anchors(p)
        s = aw.anchor_summary(rows)
        self.assertEqual(s["rows"], 4)
        self.assertEqual(s["gene_anchors"], 3)
        self.assertEqual(s["distinct_genes"], 2)
        self.assertEqual(s["query_copies_per_gene"], {"1": 1, "2": 1})
        self.assertEqual(s["strand"], {"+": 2, "-": 1})
        self.assertEqual(s["blocks"], 2)

    def test_anchors_bad_header(self):
        p = os.path.join(self.d, "anchors")
        with open(p, "w") as fh:
            fh.write("a\tb\n1\t2\n")
        with self.assertRaises(aw.StageError):
            aw.parse_anchors(p)

    def test_maf(self):
        p = os.path.join(self.d, "x.maf")
        with open(p, "w") as fh:
            fh.write("##maf version=1\na score=1\ns chrA 0 10 + 100 ACGTACGTAC\ns ctg1 5 9 + 50 ACGTAC-TAC\n\n"
                     "a score=2\ns chrA 20 4 + 100 ACGT\ns ctg2 0 4 - 9 ACGT\n\n")
        m = aw.maf_summary(p)
        self.assertEqual(m, {"blocks": 2, "ref_aligned_bp": 14, "query_aligned_bp": 13, "query_sequences": 2})

    def test_maf_empty_fails(self):
        p = os.path.join(self.d, "x.maf")
        open(p, "w").close()
        with self.assertRaises(aw.StageError):
            aw.maf_summary(p)

    def test_sam_summary(self):
        p = os.path.join(self.d, "x.sam")
        with open(p, "w") as fh:
            fh.write("@SQ\tSN:c1\tLN:100\n@SQ\tSN:c2\tLN:100\n"
                     "q1\t0\tc1\t5\t60\t*\n"
                     "q1\t256\tc2\t9\t0\t*\n"
                     "q1\t2048\tc2\t50\t0\t*\n"
                     "q2\t4\t*\t0\t0\t*\n"
                     "q3\t16\tc2\t7\t60\t*\n")
        s = aw.sam_summary(p)
        self.assertEqual(s["sq_lines"], 2)
        self.assertEqual(s["queries"], 3)
        self.assertEqual(s["mapped_queries"], 2)
        self.assertEqual(s["hits_per_query"], {"0": 1, "1": 1, "2": 1})

    def test_index_split_detected(self):
        p = os.path.join(self.d, "mm.log")
        with open(p, "w") as fh:
            fh.write("[M::main::1.0*1.00] loaded/built the index for 3 target sequence(s)\n" * 2)
        with self.assertRaises(aw.StageError):
            aw.check_single_index_part(p)

    def test_config_rejects_compressed_query(self):
        data = {"ref_fa": "/x/ref.fa", "gff": "/x/ref.gff3", "asm_fa": "/x/asm.fa.gz"}
        p = write_config(self.d, data)
        with self.assertRaises(SystemExit):
            aw.Config(p)


class _PassThrough:
    """Runtime stand-in: runs the command directly instead of in the container."""

    def wrap(self, argv, tmp):
        return argv


def _res(**kw):
    r = {"threads": 1, "cpus": "0", "nice": 10, "max_rss_gb": 8, "min_free_ram_gb": 1, "reserve_disk_gb": 1}
    r.update(kw)
    return r


class TestRunCmdSafety(unittest.TestCase):
    """Memory cap and signal handling of run_cmd, with a plain Python child (no container)."""

    def setUp(self):
        self.d = tmpdir()
        os.makedirs(os.path.join(self.d, "tmp"))

    def tearDown(self):
        shutil.rmtree(self.d)

    def cfg(self, **kw):
        import types
        return types.SimpleNamespace(res=_res(**kw), tmpdir=os.path.join(self.d, "tmp"))

    def alloc(self, mb):
        return [sys.executable, "-c", f"b = b'\\x01' * ({mb} << 20); print(len(b))"]

    def need_scope(self, cfg):
        try:
            aw.check_memory_cap(cfg.res)
        except SystemExit as e:
            self.skipTest(f"no systemd user scope here: {e}")

    def test_memory_cap_kills_command_over_limit(self):
        cfg = self.cfg(memory_cap="systemd-scope", max_rss_gb=0.15)
        self.need_scope(cfg)
        with self.assertRaises(aw.StageError) as cm:
            aw.run_cmd(cfg, _PassThrough(), self.alloc(400), self.d, "big")
        self.assertIn("memory cap", str(cm.exception.args[0]))
        info = aw.run_cmd(cfg, _PassThrough(), self.alloc(20), self.d, "small")  # under the cap: unaffected
        self.assertEqual(info["exit_code"], 0)
        with open(os.path.join(self.d, "commands.sh")) as fh:
            self.assertIn("MemoryMax=150000000", fh.read())

    def test_log_survives_closed_output(self):
        class Closed:
            def write(self, s):
                raise BrokenPipeError(32, "Broken pipe")

            def flush(self):
                raise BrokenPipeError(32, "Broken pipe")
        with contextlib.redirect_stdout(Closed()):
            aw.log("the session is gone")  # must not raise

    def test_unknown_memory_cap_rejected(self):
        data = {"ref_fa": "/x/ref.fa", "gff": "/x/ref.gff3", "asm_fa": "/x/asm.fa"}
        p = write_config(self.d, data)
        with open(p) as fh:
            c = yaml.safe_load(fh)
        c["resources"]["memory_cap"] = "ulimit"
        with open(p, "w") as fh:
            yaml.safe_dump(c, fh)
        with self.assertRaises(SystemExit):
            aw.Config(p)

    def start_runner(self, workdir, res, child_code="import time; time.sleep(600)"):
        """A runner process (as under a detached screen session) running a long child command via
        run_cmd in `workdir`; returns (runner Popen, child pid)."""
        import subprocess
        import time
        os.makedirs(workdir, exist_ok=True)
        pidfile = os.path.join(self.d, "child.pid")
        child = [sys.executable, "-c",
                 f"import os; open({pidfile!r}, 'w').write(str(os.getpid()))\n{child_code}"]
        runner = ("import signal, sys, types\n"
                  f"sys.path.insert(0, {os.path.join(HERE, '..', 'scripts')!r})\n"
                  "import g2f_anchorwave as aw\n"
                  "for s in (signal.SIGTERM, signal.SIGHUP): signal.signal(s, aw._on_signal)\n"
                  f"cfg = types.SimpleNamespace(res={res!r}, tmpdir={os.path.join(self.d, 'tmp')!r})\n"
                  "class P:\n    def wrap(self, a, t): return a\n"
                  "try:\n"
                  f"    aw.run_cmd(cfg, P(), {child!r}, {workdir!r}, 'long')\n"
                  "except aw.Interrupted:\n    sys.exit(143)\n")
        r = subprocess.Popen([sys.executable, "-c", runner])
        for _ in range(200):
            if os.path.exists(pidfile) and os.path.getsize(pidfile):
                break
            time.sleep(0.1)
        with open(pidfile) as fh:
            return r, int(fh.read())

    def test_sigterm_to_runner_stops_the_command(self):
        # SIGTERM to the runner must stop the child's process group, not orphan it
        import signal
        r, cpid = self.start_runner(self.d, _res())
        r.send_signal(signal.SIGTERM)
        self.assertEqual(r.wait(timeout=60), 143)
        with self.assertRaises(ProcessLookupError):
            os.kill(cpid, 0)

    def test_memory_cap_covers_subprocesses_and_disables_swap(self):
        import subprocess
        cfg = self.cfg(memory_cap="systemd-scope", max_rss_gb=0.15)
        self.need_scope(cfg)
        self.assertIn("MemorySwapMax=0 (verified)", aw.check_memory_cap(cfg.res))
        # the limit is on the scope's cgroup: a grandchild's allocation counts
        grandchild = [sys.executable, "-c",
                      "import subprocess, sys; sys.exit(subprocess.call([sys.executable, '-c', "
                      "\"b = b'\\\\x01' * (400 << 20)\"]))"]
        with self.assertRaises(aw.StageError) as cm:
            aw.run_cmd(cfg, _PassThrough(), grandchild, self.d, "tree")
        # 247 = the parent's sys.exit(-9): the grandchild itself was SIGKILLed (uncapped it exits 0)
        self.assertEqual(cm.exception.args[1]["exit_code"], 247)
        # logs and the process record are kept for diagnosis
        for f in ("tree.log", "tree.time.txt", "process.json", "commands.sh"):
            self.assertTrue(os.path.exists(os.path.join(self.d, f)), f)
        # the cgroup of a capped command really carries the limits
        out = subprocess.run(aw.memory_cap_argv(cfg.res) + ["sh", "-c", aw.CGROUP_PROBE],
                             capture_output=True, text=True, check=True).stdout.split()
        self.assertLess(abs(int(out[0]) - 150_000_000), 4096)
        self.assertEqual(out[1], "0")

    def test_memory_cap_fails_closed(self):
        res = _res(memory_cap="systemd-scope", max_rss_gb=1)
        orig = aw.memory_cap_argv
        try:
            for fake in (["/nonexistent/systemd-run"], ["false"], ["env"]):  # missing / fails / no cgroup limits
                aw.memory_cap_argv = lambda r, unit=None, fake=fake: fake
                with self.assertRaises(SystemExit):
                    aw.check_memory_cap(res)
        finally:
            aw.memory_cap_argv = orig

    def check_sigkill_orphan(self, res):
        """SIGKILL to the runner cannot be handled: the command survives. The next status/run must detect
        it by its recorded identity and refuse to start; nothing is signalled automatically."""
        import signal
        data = {"ref_fa": "/x/ref.fa", "gff": "/x/ref.gff3", "asm_fa": "/x/asm.fa"}
        cfgp = write_config(self.d, data)
        cfg = aw.Config(cfgp)
        st = {s.sid: s for s in aw.build_stages(cfg)}["ref_cds"]
        r, cpid = self.start_runner(st.partial, res)
        try:
            r.send_signal(signal.SIGKILL)
            r.wait(timeout=60)
            os.kill(cpid, 0)  # the command is still alive: an orphan
            pids, how = aw.orphan_processes(st.partial)
            self.assertIn(cpid, pids)
            self.assertEqual(st.state(cfg), "orphaned")
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                aw.main(["-c", cfgp, "status"])
            self.assertIn("ORPHANED ref_cds", out.getvalue())
            with self.assertRaises(SystemExit) as cm:
                aw.main(["-c", cfgp, "run"])
            self.assertIn("refusing to run", str(cm.exception))
            self.assertIn(str(cpid), str(cm.exception))
            self.assertTrue(os.path.isdir(st.partial))  # the attempt directory was not moved
        finally:
            # cleanup of the test's own child only (pid recorded by the child itself)
            rec = aw.read_json(os.path.join(st.partial, "process.json"))
            os.killpg(rec["pgid"], signal.SIGKILL)
            for _ in range(100):
                if not aw.orphan_processes(st.partial)[0]:
                    break
                import time
                time.sleep(0.1)
        self.assertEqual(aw.orphan_processes(st.partial)[0], [])
        self.assertEqual(st.state(cfg), "interrupted")  # gone: a rerun may now move the attempt aside
        return how

    def test_sigkill_orphan_detected_by_process_group(self):
        self.assertIn("process group", self.check_sigkill_orphan(_res()))

    def test_sigkill_orphan_detected_by_scope(self):
        res = _res(memory_cap="systemd-scope", max_rss_gb=1)
        self.need_scope(self.cfg(**res))
        self.assertIn(".scope", self.check_sigkill_orphan(res))


@unittest.skipUnless(IMAGE and os.path.exists(IMAGE or ""), "ANCHORWAVE_IMAGE not set")
class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = tmpdir()
        cls.data = anchorwave_data.build(cls.root)
        cls.cfgp = write_config(cls.root, cls.data)
        cls.cfg = aw.Config(cls.cfgp)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.rc = aw.main(["-c", cls.cfgp, "run"])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root)

    def rec(self, sid):
        st = {s.sid: s for s in aw.build_stages(self.cfg)}[sid]
        return aw.read_json(os.path.join(st.outdir, "VALIDATED.json")), st

    def anchors(self, variant):
        _, st = self.rec(f"anchors:SAMPLE_A:unpurged:{variant}")
        return [r for r in aw.parse_anchors(os.path.join(st.outdir, "anchors")) if r["gene"] != "interanchor"]

    def test_all_validated(self):
        self.assertEqual(self.rc, 0)
        for st in aw.build_stages(self.cfg):
            self.assertEqual(st.state(self.cfg), "validated", st.sid)

    def test_run_record_provenance(self):
        runs = os.path.join(self.cfg.outdir, "runs")
        rids = sorted(os.listdir(runs))  # other tests rerun the pipeline: the first run is setUpClass's
        rid = rids[0]
        run = aw.read_json(os.path.join(runs, rid, "RUN.json"))
        self.assertEqual(run["run_id"], rid)
        self.assertEqual(run["exit_code"], 0)
        self.assertEqual(set(run["stage_states"].values()), {"validated"})
        for k in ("commit", "tracked_files_modified", "script_sha256"):
            self.assertIn(k, run["pipeline_git"])
        for k in ("image_sha256", "anchorwave", "minimap2", "apptainer"):
            self.assertTrue(run["runtime"][k], k)
        self.assertEqual(run["science_sha256"], aw.file_digest(os.path.join(runs, rid, "science.yaml")))
        self.assertIn("SAMPLE_A:unpurged", run["assemblies"])
        for st in aw.build_stages(self.cfg):
            rec = aw.read_json(os.path.join(st.outdir, "VALIDATED.json"))
            self.assertIn(rec["run_id"], rids, st.sid)
            self.assertTrue(os.path.exists(os.path.join(st.outdir, "process.json")), st.sid)

    def test_cds_and_self_placement(self):
        r, _ = self.rec("ref_cds")
        self.assertEqual(r["validation"]["cds_records"], self.data["genes"])
        r, _ = self.rec("ref_map")
        self.assertEqual(r["validation"]["primary_at_own_locus_frac"], 1.0)
        self.assertEqual(r["runtime"]["anchorwave"], "anchorwave v1.3.1")

    def test_retained_allele_needs_R2(self):
        r1 = {a["queryChr"] for a in self.anchors("R1Q1e1")}
        r2 = {a["queryChr"] for a in self.anchors("R2Q1e1")}
        self.assertNotIn("ctg4", r1)
        self.assertIn("ctg4", r2)

    def test_chimeric_contig_two_chromosomes(self):
        refs = {a["refChr"] for a in self.anchors("R2Q1e1") if a["queryChr"] == "ctg5"}
        self.assertEqual(refs, {"chrA", "chrB"})

    def test_small_inversion_needs_small_D(self):
        inv = lambda v: [a for a in self.anchors(v) if a["queryChr"] == "ctg2" and a["strand"] == "-"]
        self.assertEqual(inv("R2Q1e1"), [])           # default -D 25 bridges a 10-gene inversion
        self.assertEqual(len(inv("R2Q1e1D5")), 10)    # all 10 inverted genes with -D 5

    def test_alignment_maf(self):
        r, st = self.rec("align:SAMPLE_A:unpurged:R2Q1e1")
        self.assertGreater(r["validation"]["maf"]["blocks"], 0)
        self.assertGreater(r["validation"]["maf"]["query_aligned_bp"], 300_000)
        self.assertTrue(os.path.exists(os.path.join(st.outdir, "alignment.f.maf")))

    def test_synteny_summary(self):
        import g2f_anchorwave_synteny as syn
        a = lambda v: self.rec(f"anchors:SAMPLE_A:unpurged:{v}")[1].outdir + "/anchors"
        out = os.path.join(self.root, "synteny")
        with contextlib.redirect_stdout(io.StringIO()):
            syn.main(["--anchors", a("R2Q1e1"), "--variant", f"D5={a('R2Q1e1D5')}",
                      "--query-fai", self.data["asm_fa"] + ".fai", "--ref-fai", self.data["ref_fa"] + ".fai",
                      "--maf", self.rec("align:SAMPLE_A:unpurged:R2Q1e1")[1].outdir + "/alignment.maf",
                      "--zoom", "chrA:150000-300000", "--out", out])
        with open(os.path.join(out, "summary.json")) as fh:
            s = json.load(fh)
        self.assertEqual(s["contigs"]["contigs_multi_chrom_ge3"], 1)       # ctg5, the chimera
        def rows(f):
            with open(os.path.join(out, f)) as fh:
                return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
        j = rows("junctions.tsv")
        self.assertIn(["ctg5", "inter_chromosomal"], [[r[0], r[4]] for r in j[1:]])
        jd5 = rows("junctions_D5.tsv")
        self.assertEqual(sorted(r[4] for r in jd5[1:] if r[0] == "ctg2"), ["strand_switch", "strand_switch"])
        d = rows("duplicated_regions.tsv")
        h = d[0]
        pairs = {(r[h.index("copy1")].split(":")[0], r[h.index("copy2")].split(":")[0]) for r in d[1:]}
        self.assertTrue(any("ctg4" in p for p in pairs))
        self.assertTrue(os.path.exists(os.path.join(out, "dotplot_genome.png")))
        self.assertGreater(s["maf"]["aligned_column_identity"], 0.98)
        with self.assertRaises(SystemExit):     # never overwrites an existing output directory
            syn.main(["--anchors", a("R2Q1e1"), "--query-fai", self.data["asm_fa"] + ".fai",
                      "--ref-fai", self.data["ref_fa"] + ".fai", "--out", out])

    def test_rerun_skips_and_stale(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(aw.main(["-c", self.cfgp, "run"]), 0)
        self.assertEqual(out.getvalue().count("validated, skipped"), len(aw.build_stages(self.cfg)))
        # a science change makes dependent stages stale; they are not rerun without --rerun-stale
        root2 = os.path.join(self.root, "stale")
        os.makedirs(root2)
        sci = dict(self.cfg.science)
        sci["minimap2"] = {"args": sci["minimap2"]["args"] + " -G 20000"}
        cp = write_config(root2, self.data, {"minimap2": sci["minimap2"]})
        with open(cp) as fh:
            c = yaml.safe_load(fh)
        c["outdir"] = os.path.dirname(self.cfg.outdir)
        c["profile"] = self.cfg.profile
        with open(cp, "w") as fh:
            yaml.safe_dump(c, fh)
        cfg2 = aw.Config(cp)
        states = {s.sid: s.state(cfg2) for s in aw.build_stages(cfg2)}
        self.assertEqual(states["ref_cds"], "validated")
        self.assertEqual(states["ref_map"], "stale")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(aw.main(["-c", cp, "run", "--stages", "ref_map"]), 1)
        self.assertIn("STALE", out.getvalue())


@unittest.skipUnless(IMAGE and os.path.exists(IMAGE or ""), "ANCHORWAVE_IMAGE not set")
class TestFailurePath(unittest.TestCase):
    def test_failed_stage_kept_then_attempt_moved(self):
        root = tmpdir()
        try:
            data = anchorwave_data.build(root)
            cp = write_config(root, data, {"validation": {"min_ref_self_placement": 0.9,
                                                          "min_query_mapped": 1.01}})
            with contextlib.redirect_stdout(io.StringIO()):
                rc = aw.main(["-c", cp, "run", "--stages", "ref_cds,query_map,anchors"])
            self.assertEqual(rc, 1)
            cfg = aw.Config(cp)
            st = {s.sid: s for s in aw.build_stages(cfg)}
            q = st["query_map:SAMPLE_A:unpurged"]
            self.assertEqual(q.state(cfg), "failed")
            rec = aw.read_json(os.path.join(q.partial, "FAILED.json"))
            self.assertIn("map to the assembly", rec["error"])
            self.assertEqual(st["anchors:SAMPLE_A:unpurged:R2Q1e1"].state(cfg), "planned")
            with contextlib.redirect_stdout(io.StringIO()):
                aw.main(["-c", cp, "run", "--stages", "query_map"])
            kept = [d for d in os.listdir(os.path.dirname(q.outdir)) if d.startswith("query_map.attempt-")]
            self.assertEqual(len(kept), 1)
            self.assertTrue(os.path.exists(os.path.join(os.path.dirname(q.outdir), kept[0], "FAILED.json")))
        finally:
            shutil.rmtree(root)

    def test_invalid_final_dir_refused(self):
        root = tmpdir()
        try:
            data = anchorwave_data.build(root)
            cp = write_config(root, data)
            cfg = aw.Config(cp)
            st = {s.sid: s for s in aw.build_stages(cfg)}["ref_cds"]
            os.makedirs(st.outdir)
            open(os.path.join(st.outdir, "something"), "w").close()
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                aw.main(["-c", cp, "run", "--stages", "ref_cds"])
            self.assertTrue(os.path.exists(os.path.join(st.outdir, "something")))
        finally:
            shutil.rmtree(root)


if __name__ == "__main__":
    unittest.main()
