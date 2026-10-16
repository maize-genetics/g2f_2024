#!/usr/bin/env python3
"""Assembly structural validation workflow.

  g2f_asmqc.py -c CONFIG check-env            tools, Python modules, inputs
  g2f_asmqc.py -c CONFIG plan [opts]           states, estimates; runs nothing
  g2f_asmqc.py -c CONFIG run  [opts]           run planned stages in order
  g2f_asmqc.py -c CONFIG status [opts]         state of every stage
  g2f_asmqc.py -c CONFIG report [opts]         write a new report per sample

Options for plan/run/status/report:
  --samples S1,S2      restrict to samples
  --stages T1,T2       restrict to stage types (prepare, graph, ref_index,
                       ref_align, ref_analysis, pair_align, map_hifi,
                       read_support, regions, evidence, pair_analysis)
  --rerun-stale        (run) keep stale results as <dir>.previous-<stamp>
                       and rerun them
  --keep-going         (run) continue with independent stages after a failure

Run it detached so the session closes when the run ends, e.g.
  run_detached asmqc LOGDIR -- python3 scripts/g2f_asmqc.py -c CONFIG run
See docs/assembly_validation.md.
"""
import argparse
import fcntl
import importlib
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from asmqc import config as C  # noqa: E402
from asmqc import report, runner, stages as S  # noqa: E402


def select(st, samples, types):
    out = []
    for sid in S.order(st):
        parts = sid.split(":")
        if types and parts[0] not in types:
            continue
        if samples and parts[0] != "ref_index" and parts[1] not in samples:
            continue
        out.append(sid)
    return out


def cmd_status(cfg, st, sids, show_est=False, ctx=None):
    tot = {"disk_final_gb": 0.0, "disk_tmp_gb": 0.0}
    print(f"profile {cfg['profile']}  root {cfg['root']}  science {cfg['science_hash']}")
    for sid in sids:
        state, why = runner.stage_state(st[sid])
        line = f"  {state:10s} {sid}"
        if why and state != "validated":
            line += f"  ({why})"
        if show_est and state != "validated":
            e = st[sid].estimate(ctx) or {}
            if e:
                line += "  est: " + ", ".join(f"{k}={round(v, 1) if isinstance(v, float) else v}" for k, v in e.items())
                tot["disk_final_gb"] += e.get("disk_final_gb", 0)
                tot["disk_tmp_gb"] = max(tot["disk_tmp_gb"], e.get("disk_tmp_gb", 0))
        print(line)
    if show_est:
        print(f"planned stages: final outputs ~{tot['disk_final_gb']:.0f} GB, peak temporary ~{tot['disk_tmp_gb']:.0f} GB "
              f"(stages run one at a time); free now {runner.free_disk_gb(cfg['root']):.0f} GB, "
              f"RAM available {runner.available_ram_gb():.0f} GB")


def cmd_check_env(cfg, ctx):
    ok = True
    print("tools:")
    for name, v in ctx.tool_versions().items():
        flag = "OK " if v["exit_code"] == 0 else "ERR"
        ok &= v["exit_code"] == 0
        print(f"  {flag} {name:10s} {v['reported']}   ({v['cmd']})")
    print("python modules:")
    for m in ("numpy", "pysam", "matplotlib", "networkx", "yaml"):
        try:
            mod = importlib.import_module(m)
            print(f"  OK  {m:10s} {getattr(mod, '__version__', '')}")
        except ImportError:
            ok = False
            print(f"  ERR {m:10s} missing")
    print(f"  python {sys.version.split()[0]} ({sys.executable})")
    print("inputs:")
    paths = [r["fasta"] for r in cfg["references"].values()]
    for s in cfg["samples"].values():
        paths += s["hifi_reads"]
        for a in s["assemblies"].values():
            paths += [p for p in (a["fasta"], a["ctg_gfa"], a["utg_gfa"], a["assembler_log"]) if p]
            paths += a["extra_targets"]
            if a["fasta"].endswith(".gz"):
                paths += [a["fasta"] + ".fai", a["fasta"] + ".gzi"]
    for p in paths:
        e = os.path.exists(p)
        ok &= e
        print(f"  {'OK ' if e else 'MISSING'} {p}")
    print(f"free disk at output root: {runner.free_disk_gb(cfg['root']):.0f} GB; "
          f"tmpdir {cfg['tmpdir']}: {runner.free_disk_gb(cfg['tmpdir']):.0f} GB; "
          f"RAM available {runner.available_ram_gb():.0f} GB; CPUs {os.cpu_count()}")
    print("environment check:", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


def cmd_run(cfg, st, sids, ctx, rerun_stale, keep_going):
    os.makedirs(cfg["root"], exist_ok=True)
    lockfh = open(os.path.join(cfg["root"], ".runner.lock"), "w")
    try:
        fcntl.flock(lockfh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit(f"another g2f_asmqc run holds {cfg['root']}/.runner.lock")
    print(f"run start {runner.now()} pid {os.getpid()} git {runner.git_commit(ctx.pipeline_dir)}", flush=True)
    failed, n_done = set(), 0
    for sid in sids:
        stage = st[sid]
        state, why = runner.stage_state(stage)
        if state == "validated":
            print(f"[skip] {sid}: validated", flush=True)
            continue
        if state == "stale" and rerun_stale:
            runner.move_stale(ctx, stage)
        bad = [d for d in stage.deps if d in failed or (d in st and runner.stage_state(st[d])[0] != "validated")]
        if bad:
            print(f"[wait] {sid}: dependencies not validated: {', '.join(bad)}", flush=True)
            failed.add(sid)
            continue
        print(f"[run ] {sid} {runner.now()}", flush=True)
        try:
            runner.execute(ctx, stage)
            n_done += 1
            print(f"[ok  ] {sid} {runner.now()}", flush=True)
        except runner.StageError as e:
            print(f"[FAIL] {e}", flush=True)
            failed.add(sid)
            if not keep_going:
                break
    print(f"run end {runner.now()}: {n_done} stage(s) validated, {len(failed)} not completed", flush=True)
    lockfh.close()
    return 1 if failed else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True)
    ap.add_argument("command", choices=["check-env", "plan", "run", "status", "report"])
    ap.add_argument("--samples", default="")
    ap.add_argument("--stages", default="")
    ap.add_argument("--rerun-stale", action="store_true")
    ap.add_argument("--keep-going", action="store_true")
    args = ap.parse_args(argv)
    try:
        cfg = C.load(args.config)
    except C.ConfigError as e:
        sys.exit(f"configuration error: {e}")
    samples = [x for x in args.samples.split(",") if x]
    types = [x for x in args.stages.split(",") if x]
    for t in types:
        if t not in S.TYPE_ORDER:
            sys.exit(f"unknown stage type {t}")
    for s in samples:
        if s not in cfg["samples"]:
            sys.exit(f"unknown sample {s}")
    ctx = runner.Context(cfg, C.PIPELINE_DIR)
    st = S.build(cfg, samples or None)
    sids = select(st, samples, types)
    if args.command == "check-env":
        return cmd_check_env(cfg, ctx)
    if args.command in ("status", "plan"):
        cmd_status(cfg, st, sids, show_est=args.command == "plan", ctx=ctx)
        if args.command == "plan":
            print("plan only: nothing was run")
        return 0
    if args.command == "run":
        return cmd_run(cfg, st, sids, ctx, args.rerun_stale, args.keep_going)
    if args.command == "report":
        stamp = runner.stamp()
        for s in samples or cfg["samples"]:
            print(report.run(cfg, st, s, stamp))
        return 0


if __name__ == "__main__":
    sys.exit(main())
