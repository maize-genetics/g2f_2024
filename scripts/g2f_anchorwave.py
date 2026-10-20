#!/usr/bin/env python3
"""Gene-anchored genome alignment of contig assemblies against a reference
with AnchorWave (proali), as resumable stages with explicit states.

Stages (in order):
  ref_cds                   anchorwave gff2seq: longest full-length CDS per gene
  ref_map                   minimap2 splice mapping of those CDS to the reference
  query_map:S:A             the same CDS mapped to assembly A of sample S
  anchors:S:A:V             proali collinear anchors only (-n, -ns), variant V
  align:S:A:V               proali anchors + base-level alignment (MAF), variant V

States (from files in a stage's output directory):
  planned    no output directory
  running    <dir>.partial exists and the run lock is held
  failed     <dir>.partial holds FAILED.json (kept for diagnosis)
  validated  <dir>/VALIDATED.json matches the current inputs and parameters
  stale      VALIDATED.json exists but inputs or parameters changed

A stage runs in <dir>.partial and becomes <dir> only after its outputs
validate. Nothing is deleted: an earlier failed or interrupted attempt is
moved to <dir>.attempt-<stamp>, and --rerun-stale moves an old result to
<dir>.previous-<stamp>.

  g2f_anchorwave.py -c CONFIG check-env
  g2f_anchorwave.py -c CONFIG plan   [--stages ...] [--samples ...]
  g2f_anchorwave.py -c CONFIG run    [--stages ...] [--samples ...] [--rerun-stale]
  g2f_anchorwave.py -c CONFIG status

--stages takes stage kinds (ref_cds,ref_map,query_map,anchors,align) or full
stage ids; dependencies must already be validated or be selected too.
"""
import argparse
import datetime
import fcntl
import getpass
import gzip
import hashlib
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time

import yaml

SCHEMA = "g2f-anchorwave-stage/1"
CODE_VERSION = 1  # bump when stage outputs change meaning
KINDS = ("ref_cds", "ref_map", "query_map", "anchors", "align")
# sha256 of this file as loaded by this process (stages record the code that actually ran,
# even if the file is edited while a run is in progress)
with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()
ANCHOR_COLUMNS = ("refChr", "referenceStart", "referenceEnd", "queryChr", "queryStart",
                  "queryEnd", "strand", "gene", "blockIndex")


class StageError(Exception):
    pass


class Interrupted(BaseException):
    """SIGTERM/SIGHUP received by the runner: stop the running stage, record it, end the run."""


def _on_signal(signum, frame):
    raise Interrupted(signal.Signals(signum).name)


# ---------------------------------------------------------------- utilities

def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def stamp():
    return datetime.datetime.now().strftime("%Y%m%dT%H%M%S")


def log(msg):
    # the terminal or log pipe may be gone (e.g. its screen session was closed); records in the
    # stage and run directories, not this stream, are what must survive
    try:
        print(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)
    except OSError:
        pass


def write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def short_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()[:12]


def file_digest(path, algo="sha256"):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def file_signature(paths):
    sig = []
    for p in paths:
        try:
            st = os.stat(p)
            sig.append([p, st.st_size, int(st.st_mtime)])
        except OSError:
            sig.append([p, None, None])
    return sig


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def free_disk_gb(path):
    while not os.path.exists(path):
        path = os.path.dirname(path)
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def available_ram_gb():
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1e6
    return float("nan")


def parse_cpus(spec):
    cpus = set()
    for part in str(spec).split(","):
        a, _, b = part.partition("-")
        cpus.update(range(int(a), int(b or a) + 1))
    return sorted(cpus)


def cpu_busy_fraction(cpus, seconds=2.0):
    """Mean busy fraction of the given CPUs over a short interval (/proc/stat)."""
    def sample():
        out = {}
        with open("/proc/stat") as fh:
            for line in fh:
                m = re.match(r"cpu(\d+) (.*)", line)
                if m:
                    v = [int(x) for x in m.group(2).split()]
                    out[int(m.group(1))] = (sum(v), v[3] + v[4])
        return out
    a = sample()
    time.sleep(seconds)
    b = sample()
    fr = []
    for c in cpus:
        tot = b[c][0] - a[c][0]
        idle = b[c][1] - a[c][1]
        fr.append(1 - idle / tot if tot else 0.0)
    return sum(fr) / len(fr) if fr else 0.0


def git_state(pipeline_dir):
    r = subprocess.run(["git", "-C", pipeline_dir, "rev-parse", "HEAD"], capture_output=True, text=True)
    d = subprocess.run(["git", "-C", pipeline_dir, "status", "--porcelain", "--untracked-files=no",
                        "--", "scripts/g2f_anchorwave.py"], capture_output=True, text=True)
    t = subprocess.run(["git", "-C", pipeline_dir, "status", "--porcelain", "--untracked-files=no"],
                       capture_output=True, text=True)
    return {"commit": r.stdout.strip() or None,
            "script_modified_or_untracked": bool(d.stdout.strip()) or not _tracked(pipeline_dir),
            "tracked_files_modified": [l[3:] for l in t.stdout.splitlines()],
            "script_sha256": SCRIPT_SHA256,
            "script_sha256_on_disk_now": file_digest(os.path.abspath(__file__))}


def _tracked(pipeline_dir):
    r = subprocess.run(["git", "-C", pipeline_dir, "ls-files", "--error-unmatch", "scripts/g2f_anchorwave.py"],
                       capture_output=True, text=True)
    return r.returncode == 0


def fasta_index(path):
    """Name -> length from <path>.fai (required: it proves uniform line lengths)."""
    fai = path + ".fai"
    if not os.path.exists(fai):
        raise StageError(f"missing FASTA index {fai} (samtools faidx)")
    out = {}
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = int(f[1])
    return out


# ---------------------------------------------------------------- config

class Config:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(path) as fh:
            c = yaml.safe_load(fh)
        self.raw = c
        for key in ("outdir", "tmpdir", "science", "runtime", "resources", "reference", "samples"):
            if key not in c:
                raise SystemExit(f"config {path}: missing '{key}'")
        self.profile = c.get("profile", "anchorwave_v1")
        self.outdir = os.path.join(c["outdir"], self.profile)
        self.tmpdir = c["tmpdir"]
        sci = c["science"]
        self.science_file = sci if os.path.isabs(sci) else os.path.join(self.pipeline_dir, sci)
        with open(self.science_file) as fh:
            self.science = yaml.safe_load(fh)
        self.runtime = c["runtime"]
        self.res = c["resources"]
        self.ref = c["reference"]
        self.samples = c["samples"]
        self.memory_cap = None  # set by the run command after check_memory_cap()
        self.run_id = None      # set by the run command
        self.holding_lock = False
        self._validate()

    def _validate(self):
        errs = []
        for k in ("id", "fasta", "gff3"):
            if k not in self.ref:
                errs.append(f"reference.{k} missing")
        if "image" not in self.runtime:
            errs.append("runtime.image missing")
        for k in ("threads", "cpus", "max_rss_gb", "min_free_ram_gb", "reserve_disk_gb"):
            if k not in self.res:
                errs.append(f"resources.{k} missing")
        if self.res.get("memory_cap", "none") not in ("none", "systemd-scope"):
            errs.append("resources.memory_cap must be 'none' or 'systemd-scope'")
        variants = self.science.get("proali", {}).get("variants", {})
        for name, v in variants.items():
            for k in ("R", "Q", "e"):
                if not isinstance(v.get(k), int) or v[k] < 1:
                    errs.append(f"science proali.variants.{name}.{k} must be a positive integer")
        for name in self.science.get("proali", {}).get("align_variants", []):
            if name not in variants:
                errs.append(f"align variant {name} not defined in proali.variants")
        for s, sd in self.samples.items():
            for a, ad in sd.get("assemblies", {}).items():
                if "fasta" not in ad:
                    errs.append(f"samples.{s}.assemblies.{a}.fasta missing")
                if ad.get("fasta", "").endswith((".gz", ".bgz")):
                    errs.append(f"samples.{s}.assemblies.{a}.fasta: use an uncompressed, faidx-indexed FASTA")
        if errs:
            raise SystemExit("config errors:\n  " + "\n  ".join(errs))

    def mem_gib(self):
        return self.res.get("anchorwave_mem_gib")


# ---------------------------------------------------------------- stages

class Stage:
    def __init__(self, sid, kind, outdir, deps, inputs, params, estimate, sample=None, asm=None, variant=None):
        self.sid, self.kind, self.outdir, self.deps = sid, kind, outdir, deps
        self.inputs, self.params, self.estimate = inputs, params, estimate
        self.sample, self.asm, self.variant = sample, asm, variant

    @property
    def partial(self):
        return self.outdir + ".partial"

    def state(self, cfg):
        v = os.path.join(self.outdir, "VALIDATED.json")
        if os.path.exists(v):
            rec = read_json(v)
            same = (rec.get("params_hash") == short_hash(self.params)
                    and rec.get("inputs") == file_signature(self.inputs)
                    and rec.get("code_version") == CODE_VERSION)
            return "validated" if same else "stale"
        if os.path.isdir(self.outdir):
            return "invalid"  # final directory without a validation record: never touched
        if os.path.isdir(self.partial):
            if os.path.exists(os.path.join(self.partial, "FAILED.json")):
                return "failed"
            if not cfg.holding_lock and run_lock_held(cfg):  # flock is per open file: skip if we hold it
                return "running"
            pids, _ = orphan_processes(self.partial)
            # None = recorded on another host: cannot prove the command is gone, treat as orphaned
            return "orphaned" if pids is None or pids else "interrupted"
        return "planned"


def run_lock_held(cfg):
    """True if another process holds the run lock of this output profile."""
    p = os.path.join(cfg.outdir, ".runner.lock")
    if not os.path.exists(p):
        return False
    with open(p) as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(fh, fcntl.LOCK_UN)
        return False


def build_stages(cfg):
    sci = cfg.science
    od = cfg.outdir
    rid = cfg.ref["id"]
    refdir = os.path.join(od, "reference", rid)
    mm_args = sci["minimap2"]["args"]
    gm = int(sci["gff2seq"].get("min_exon", 20))
    stages = []
    cds_fa = os.path.join(refdir, "cds", "cds.fa")
    ref_sam = os.path.join(refdir, "ref_map", "ref.sam")
    stages.append(Stage("ref_cds", "ref_cds", os.path.join(refdir, "cds"), [],
                        [cfg.ref["fasta"], cfg.ref["gff3"]], {"min_exon": gm, "use_exon": False},
                        {"threads": 1, "ram_gb": 8, "disk_gb": 0.2, "hours": 0.1}))
    stages.append(Stage("ref_map", "ref_map", os.path.join(refdir, "ref_map"), ["ref_cds"],
                        [cfg.ref["fasta"], cds_fa], {"minimap2": mm_args},
                        {"threads": cfg.res["threads"], "ram_gb": 60, "disk_gb": 2, "hours": 0.5}))
    variants = sci["proali"]["variants"]
    common = sci["proali"].get("common_args", [])
    for s, sd in cfg.samples.items():
        for a, ad in sd["assemblies"].items():
            base = os.path.join(od, "samples", s, a)
            qsam = os.path.join(base, "query_map", "query.sam")
            gb = sum(fasta_index(ad["fasta"]).values()) / 1e9 if os.path.exists(ad["fasta"] + ".fai") else 3.0
            stages.append(Stage(f"query_map:{s}:{a}", "query_map", os.path.join(base, "query_map"),
                                ["ref_cds"], [ad["fasta"], cds_fa], {"minimap2": mm_args},
                                {"threads": cfg.res["threads"], "ram_gb": 30 * gb, "disk_gb": 2, "hours": 0.5},
                                s, a))
            pro_inputs = [cfg.ref["gff3"], cfg.ref["fasta"], cds_fa, ref_sam, ad["fasta"], qsam]
            deps = ["ref_cds", "ref_map", f"query_map:{s}:{a}"]
            for vname, v in variants.items():
                p = {"R": v["R"], "Q": v["Q"], "e": v["e"], "args": v.get("args", []), "min_exon": gm,
                     "common_args": common, "ns": True}
                stages.append(Stage(f"anchors:{s}:{a}:{vname}", "anchors",
                                    os.path.join(base, f"anchors_{vname}"), deps, pro_inputs, p,
                                    {"threads": cfg.res["threads"], "ram_gb": 40, "disk_gb": 0.5, "hours": 0.5},
                                    s, a, vname))
            for vname in sci["proali"].get("align_variants", []):
                v = variants[vname]
                p = {"R": v["R"], "Q": v["Q"], "e": v["e"], "args": v.get("args", []), "min_exon": gm,
                     "common_args": common, "ns": False, "align_args": sci["proali"].get("align_args", [])}
                stages.append(Stage(f"align:{s}:{a}:{vname}", "align", os.path.join(base, f"align_{vname}"),
                                    deps + [f"anchors:{s}:{a}:{vname}"], pro_inputs, p,
                                    {"threads": cfg.res["threads"], "ram_gb": cfg.res["max_rss_gb"],
                                     "disk_gb": 25 * gb * v["R"] / 2, "hours": 6},
                                    s, a, vname))
    return stages


def select(stages, kinds_or_ids, samples):
    out = []
    for st in stages:
        if samples and st.sample is not None and st.sample not in samples:
            continue
        if kinds_or_ids and st.kind not in kinds_or_ids and st.sid not in kinds_or_ids:
            continue
        out.append(st)
    return out


# ---------------------------------------------------------------- execution

class Runtime:
    def __init__(self, cfg):
        self.cfg = cfg
        rt = cfg.runtime
        self.apptainer = rt.get("apptainer", "apptainer")
        self.image = rt["image"]
        self.binds = rt.get("bind", [])

    def wrap(self, argv, tmp):
        cmd = [self.apptainer, "exec", "--cleanenv", "--no-home", "--env", f"TMPDIR={tmp}"]
        for b in self.binds:
            cmd += ["-B", f"{b}:{b}"]
        return cmd + [self.image] + argv

    def versions(self):
        out = {"image": self.image, "image_sha256": None}
        out["image_sha256"] = file_digest(self.image)
        for name, argv in (("anchorwave", ["anchorwave", "--version"]), ("minimap2", ["minimap2", "--version"])):
            r = subprocess.run(self.wrap(argv, self.cfg.tmpdir), capture_output=True, text=True)
            out[name] = (r.stdout + r.stderr).strip().splitlines()[0] if (r.stdout + r.stderr).strip() else None
            out[name + "_exit"] = r.returncode
        r = subprocess.run([self.apptainer, "--version"], capture_output=True, text=True)
        out["apptainer"] = r.stdout.strip()
        return out


class Watchdog(threading.Thread):
    """Stops a stage's process group if it exceeds the RSS limit, or if free RAM
    or free disk fall below their floors. Records why in .reason."""

    def __init__(self, pgid, cfg, disk_path, interval=10):
        super().__init__(daemon=True)
        self.pgid, self.cfg, self.disk_path, self.interval = pgid, cfg, disk_path, interval
        self.reason = None
        self.peak_rss_gb = 0.0
        self._halt = threading.Event()

    def group_rss_gb(self):
        tot = 0
        for d in os.listdir("/proc"):
            if not d.isdigit():
                continue
            try:
                with open(f"/proc/{d}/stat") as fh:
                    f = fh.read().rsplit(")", 1)[1].split()
                if int(f[2]) == self.pgid:  # field 5 (pgrp)
                    tot += int(f[21]) * os.sysconf("SC_PAGE_SIZE")  # field 24 (rss pages)
            except (OSError, IndexError, ValueError):
                pass
        return tot / 1e9

    def run(self):
        r = self.cfg.res
        while not self._halt.wait(self.interval):
            rss = self.group_rss_gb()
            self.peak_rss_gb = max(self.peak_rss_gb, rss)
            why = None
            if rss > r["max_rss_gb"]:
                why = f"process-group RSS {rss:.1f} GB > max_rss_gb {r['max_rss_gb']}"
            elif available_ram_gb() < r["min_free_ram_gb"]:
                why = f"system MemAvailable {available_ram_gb():.1f} GB < min_free_ram_gb {r['min_free_ram_gb']}"
            elif free_disk_gb(self.disk_path) < r["reserve_disk_gb"] / 2:
                why = f"free disk {free_disk_gb(self.disk_path):.1f} GB < reserve_disk_gb/2"
            if why:
                self.reason = why
                log(f"watchdog: {why}; stopping process group {self.pgid}")
                try:
                    os.killpg(self.pgid, signal.SIGTERM)
                    time.sleep(30)
                    os.killpg(self.pgid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                return

    def stop(self):
        self._halt.set()


def memory_cap_argv(res, unit=None):
    """Kernel-enforced memory limit for one command: a transient systemd user scope
    with MemoryMax = max_rss_gb and no swap. The limit applies to the whole cgroup,
    i.e. the command and every process it starts. Exceeding it kills the command
    inside its own cgroup (other users' processes are never affected). The watchdog
    stays as a second layer for the system-wide floors. `unit` names the scope so
    that an orphaned command can be identified later (see orphan_processes)."""
    if res.get("memory_cap", "none") == "none":
        return []
    return (["systemd-run", "--user", "--scope", "--quiet"] + (["--unit", unit] if unit else [])
            + ["-p", f"MemoryMax={int(res['max_rss_gb'] * 1e9)}", "-p", "MemorySwapMax=0"])


CGROUP_PROBE = 'p=$(cut -d: -f3 /proc/self/cgroup); cat "/sys/fs/cgroup$p/memory.max" "/sys/fs/cgroup$p/memory.swap.max"'


def check_memory_cap(res):
    """Fail closed: a configured cap that cannot be created, or whose cgroup does not
    carry the requested limits, must not silently fall back to no cap."""
    cap = memory_cap_argv(res)
    if not cap:
        return None
    want = int(res["max_rss_gb"] * 1e9)
    try:
        r = subprocess.run(cap + ["sh", "-c", CGROUP_PROBE], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise SystemExit(f"resources.memory_cap: systemd-scope unavailable ({e}); "
                         "fix the user systemd session or set memory_cap: none explicitly")
    lines = r.stdout.split()
    if r.returncode != 0 or len(lines) != 2:
        raise SystemExit(f"resources.memory_cap: systemd-scope unavailable ({(r.stderr or r.stdout).strip()}); "
                         "fix the user systemd session or set memory_cap: none explicitly")
    mem, swap = lines
    if not mem.isdigit() or abs(int(mem) - want) >= 4096 or swap != "0":  # memory.max is page-rounded
        raise SystemExit(f"resources.memory_cap: scope cgroup has memory.max={mem}, memory.swap.max={swap}; "
                         f"expected ~{want} and 0")
    return f"systemd-scope MemoryMax={want} bytes (cgroup memory.max={mem}), MemorySwapMax=0 (verified)"


def proc_stat(pid):
    """(pgrp, starttime in clock ticks since boot) of a process, or None."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            f = fh.read().rsplit(")", 1)[1].split()
        return int(f[2]), int(f[19])
    except (OSError, IndexError, ValueError):
        return None


def orphan_processes(workdir):
    """Live processes of a command recorded in <workdir>/process.json, identified by
    their systemd scope (cgroup path contains the recorded unit) or, without a scope,
    by process group plus a start time not earlier than the recorded one. Returns
    (pids, note). Nothing is ever signalled here."""
    p = os.path.join(workdir, "process.json")
    if not os.path.exists(p):
        return [], "no process record"
    rec = read_json(p)
    if rec.get("host") != socket.gethostname():
        return None, f"recorded on host {rec.get('host')}; cannot check from here"
    unit, pgid, t0 = rec.get("scope_unit"), rec.get("pgid"), rec.get("pgid_starttime", 0)
    pids = []
    for d in os.listdir("/proc"):
        if not d.isdigit() or int(d) == os.getpid():
            continue
        st = proc_stat(d)
        if st is None:
            continue
        hit = False
        if unit:
            try:
                with open(f"/proc/{d}/cgroup") as fh:
                    hit = f"/{unit}.scope" in fh.read()
            except OSError:
                pass
        else:
            hit = st[0] == pgid and st[1] >= t0
        if hit:
            pids.append(int(d))
    return sorted(pids), (f"scope {unit}.scope" if unit else f"process group {pgid}")


def stop_group(p, grace=30):
    """SIGTERM the command's process group, SIGKILL it if still alive after `grace` seconds."""
    try:
        os.killpg(p.pid, signal.SIGTERM)
        p.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()
    except ProcessLookupError:
        pass


def run_cmd(cfg, rt, argv, workdir, name, stdout_file=None):
    """Run one containerised command under taskset/nice and /usr/bin/time, with
    the optional memory cap and the watchdog. Returns a dict with exit code,
    timings and peak RSS."""
    res = cfg.res
    tmp = os.path.join(cfg.tmpdir, f"{os.path.basename(workdir.rstrip('/'))}.{name}.{stamp()}")
    os.makedirs(tmp, exist_ok=False)
    timefile = os.path.join(workdir, f"{name}.time.txt")
    unit = None
    if memory_cap_argv(res):
        tag = re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.basename(workdir.rstrip("/")))
        unit = f"g2f-aw-{tag}-{name}-{stamp()}-{os.getpid()}"
    full = (["/usr/bin/time", "-v", "-o", timefile] + memory_cap_argv(res, unit)
            + ["nice", "-n", str(res.get("nice", 10)), "taskset", "-c", str(res["cpus"])] + rt.wrap(argv, tmp))
    with open(os.path.join(workdir, "commands.sh"), "a") as fh:
        fh.write(" ".join(shlex_quote(x) for x in full)
                 + (f" > {shlex_quote(stdout_file)}" if stdout_file else "") + "\n")
    errpath = os.path.join(workdir, f"{name}.log")
    t0 = time.time()
    with open(errpath, "w") as err, (open(stdout_file, "w") if stdout_file else open(os.devnull, "w")) as out:
        p = subprocess.Popen(full, stdout=out if stdout_file else err, stderr=err, cwd=workdir,
                             preexec_fn=os.setpgrp)
        # identity of the running command, so that a later run can detect it if this runner is
        # killed without cleanup (SIGKILL, OOM, reboot of the session) -- see orphan_processes()
        st0 = proc_stat(p.pid)
        write_json(os.path.join(workdir, "process.json"),
                   {"host": socket.gethostname(), "runner_pid": os.getpid(), "pgid": p.pid,
                    "pgid_starttime": st0[1] if st0 else 0, "scope_unit": unit, "command": name, "start": now()})
        wd = Watchdog(p.pid, cfg, workdir)
        wd.start()
        try:
            rc = p.wait()
        except BaseException:  # Ctrl-C, or SIGTERM/SIGHUP to the runner: never leave the command running
            stop_group(p)
            raise
        finally:
            wd.stop()
    info = {"name": name, "exit_code": rc, "elapsed_s": round(time.time() - t0, 1),
            "watchdog_peak_rss_gb": round(wd.peak_rss_gb, 2), "watchdog_reason": wd.reason,
            "scope_unit": unit, "tmpdir": tmp}
    if os.path.exists(timefile):
        with open(timefile) as fh:
            for line in fh:
                if "Maximum resident set size" in line:
                    info["max_rss_kb"] = int(line.split(":")[1])
                elif "Elapsed (wall clock)" in line:
                    info["wall"] = line.split("):", 1)[1].strip()
                elif "User time" in line:
                    info["user_s"] = float(line.split(":")[1])
    # the per-command tmp dir belongs to this workflow; it is removed only if empty
    try:
        os.rmdir(tmp)
        info["tmpdir"] = None
    except OSError:
        pass
    if rc != 0 or wd.reason:
        why = wd.reason or ""
        if not why and rc in (137, -9) and memory_cap_argv(res):
            why = f"killed (SIGKILL): probably the memory cap, MemoryMax {res['max_rss_gb']} GB"
        raise StageError(f"{name} failed (exit {rc}{'; ' + why if why else ''}); see {errpath}", info)
    return info


def shlex_quote(s):
    import shlex
    return shlex.quote(str(s))


# ---------------------------------------------------------------- validation helpers

def read_fasta_ids(path):
    ids = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def gff_index(gff):
    """seqids, genes with CDS, and ID -> (seqid, start, end) for every feature."""
    seqids, loc, cds_parents, region_len = set(), {}, set(), {}
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqids.add(f[0])
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] in ("chromosome", "region", "scaffold", "contig"):
                region_len[f[0]] = int(f[4])
            if "ID" in attrs:
                loc[attrs["ID"]] = (f[0], int(f[3]), int(f[4]))
            if f[2] == "CDS":
                cds_parents.update(attrs.get("Parent", "").split(","))
    return seqids, loc, cds_parents, region_len


def sam_summary(path, nrec_expected=None):
    """Per-query hit counts, primary placements and an index-split check."""
    hits, primary, sq = {}, {}, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                sq += line.startswith("@SQ")
                continue
            f = line.split("\t", 10)
            flag = int(f[1])
            if flag & 4:
                hits.setdefault(f[0], 0)
                continue
            if flag & 2048:
                continue
            hits[f[0]] = hits.get(f[0], 0) + 1
            if not flag & 256:
                primary[f[0]] = (f[2], int(f[3]))
    dist = {}
    for n in hits.values():
        k = str(n) if n < 5 else "5+"
        dist[k] = dist.get(k, 0) + 1
    return {"sq_lines": sq, "queries": len(hits), "mapped_queries": len(primary),
            "hits_per_query": dict(sorted(dist.items())), "_primary": primary, "_hits": hits}


def check_single_index_part(logpath):
    with open(logpath) as fh:
        txt = fh.read()
    parts = len(re.findall(r"loaded/built the index for \d+ target sequence", txt))
    if parts != 1:
        raise StageError(f"minimap2 index was built in {parts} parts ({logpath}); SAM would be incomplete")


def parse_anchors(path):
    """AnchorWave anchors file -> list of dicts (header line names the columns)."""
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                missing = [c for c in ANCHOR_COLUMNS if c not in header]
                if missing:
                    raise StageError(f"{path}: anchors header lacks {missing}")
                continue
            rows.append(dict(zip(header, f)))
    if header is None:
        raise StageError(f"{path}: no header")
    return rows


def anchor_summary(rows):
    genes = [r for r in rows if r["gene"] not in ("interanchor",) and not r["gene"].startswith("localAlignment")]
    per_gene = {}
    for r in genes:
        per_gene[r["gene"]] = per_gene.get(r["gene"], 0) + 1
    dist = {}
    for n in per_gene.values():
        k = str(n) if n < 4 else "4+"
        dist[k] = dist.get(k, 0) + 1
    blocks = {(r["blockIndex"], r["refChr"], r["queryChr"], r["strand"]) for r in genes}
    return {"rows": len(rows), "gene_anchors": len(genes), "distinct_genes": len(per_gene),
            "query_copies_per_gene": dict(sorted(dist.items())), "blocks": len(blocks),
            "query_sequences": len({r["queryChr"] for r in genes}),
            "strand": {s: sum(1 for r in genes if r["strand"] == s) for s in ("+", "-")}}


def maf_summary(path):
    blocks, ref_bp, qry_bp, names = 0, 0, 0, set()
    with open(path) as fh:
        in_block, row = False, 0
        for line in fh:
            if line.startswith("a"):
                blocks += 1
                in_block, row = True, 0
            elif line.startswith("s") and in_block:
                f = line.split()
                if row == 0:
                    ref_bp += int(f[3])
                else:
                    qry_bp += int(f[3])
                    names.add(f[1])
                row += 1
            elif not line.strip():
                in_block = False
    if blocks == 0:
        raise StageError(f"{path}: no alignment blocks")
    return {"blocks": blocks, "ref_aligned_bp": ref_bp, "query_aligned_bp": qry_bp,
            "query_sequences": len(names)}


# ---------------------------------------------------------------- stage bodies

def body_ref_cds(cfg, rt, st, wd):
    ref, gff = cfg.ref["fasta"], cfg.ref["gff3"]
    checks = {}
    for key, path in (("fasta_md5", ref), ("gff3_md5", gff)):
        exp = cfg.ref.get(key)
        if exp:
            got = file_digest(path, "md5")
            checks[key] = {"expected": exp, "observed": got}
            if got != exp:
                raise StageError(f"{path}: md5 {got} != expected {exp}")
    fa = fasta_index(ref)
    seqids, loc, cds_parents, region_len = gff_index(gff)
    absent = sorted(seqids - set(fa))
    if absent:
        raise StageError(f"GFF3 sequence ids not in FASTA: {absent[:10]}")
    bad_len = sorted(s for s, n in region_len.items() if fa.get(s) != n)
    if bad_len:
        raise StageError(f"GFF3 region lengths differ from FASTA for {bad_len[:10]}")
    out = os.path.join(wd, "cds.fa")
    info = run_cmd(cfg, rt, ["anchorwave", "gff2seq", "-i", gff, "-r", ref, "-o", out,
                             "-m", str(st.params["min_exon"])], wd, "gff2seq")
    ids = read_fasta_ids(out)
    if not ids:
        raise StageError("gff2seq wrote no sequences")
    if len(ids) != len(set(ids)):
        raise StageError("duplicate ids in cds.fa")
    unlocated = [i for i in ids if i not in loc]
    return {"commands": [info], "outputs": ["cds.fa"],
            "validation": {**checks, "gff_seqids": len(seqids), "gff_seqids_in_fasta": True,
                           "gff_region_lengths_match": len(region_len),
                           "transcripts_with_cds": len(cds_parents), "cds_records": len(ids),
                           "cds_ids_not_in_gff": len(unlocated)}}


def _map(cfg, rt, st, wd, target, out_name):
    cds = os.path.join(cfg.outdir, "reference", cfg.ref["id"], "cds", "cds.fa")
    argv = ["minimap2"] + st.params["minimap2"].split() + ["-t", str(cfg.res["threads"]), target, cds]
    out = os.path.join(wd, out_name)
    info = run_cmd(cfg, rt, argv, wd, "minimap2", stdout_file=out)
    check_single_index_part(os.path.join(wd, "minimap2.log"))
    s = sam_summary(out)
    n = len(read_fasta_ids(cds))
    if s["queries"] != n:
        raise StageError(f"{out_name}: {s['queries']} CDS in SAM, {n} in cds.fa")
    if s["sq_lines"] != len(fasta_index(target)):
        raise StageError(f"{out_name}: @SQ lines {s['sq_lines']} != target sequences")
    return info, s, n


def body_ref_map(cfg, rt, st, wd):
    info, s, n = _map(cfg, rt, st, wd, cfg.ref["fasta"], "ref.sam")
    _, loc, _, _ = gff_index(cfg.ref["gff3"])
    own = sum(1 for q, (chrom, pos) in s["_primary"].items()
              if q in loc and loc[q][0] == chrom and loc[q][1] - 1000 <= pos <= loc[q][2] + 1000)
    frac = own / n
    minf = cfg.science["validation"]["min_ref_self_placement"]
    if frac < minf:
        raise StageError(f"only {frac:.3f} of CDS map back to their own locus (< {minf})")
    return {"commands": [info], "outputs": ["ref.sam"],
            "validation": {"cds": n, "mapped": s["mapped_queries"], "primary_at_own_locus": own,
                           "primary_at_own_locus_frac": round(frac, 4), "hits_per_cds": s["hits_per_query"]}}


def body_query_map(cfg, rt, st, wd):
    ad = cfg.samples[st.sample]["assemblies"][st.asm]
    checks = {}
    if ad.get("expected_sha256"):
        got = file_digest(ad["fasta"])
        checks["fasta_sha256"] = {"expected": ad["expected_sha256"], "observed": got}
        if got != ad["expected_sha256"]:
            raise StageError(f"{ad['fasta']}: sha256 {got} != expected")
    idx = fasta_index(ad["fasta"])
    for k, obs in (("expected_contigs", len(idx)), ("expected_bp", sum(idx.values()))):
        if ad.get(k) is not None and ad[k] != obs:
            raise StageError(f"{ad['fasta']}: {k} {ad[k]} != observed {obs}")
        checks[k] = obs
    info, s, n = _map(cfg, rt, st, wd, ad["fasta"], "query.sam")
    frac = s["mapped_queries"] / n
    minf = cfg.science["validation"]["min_query_mapped"]
    if frac < minf:
        raise StageError(f"only {frac:.3f} of CDS map to the assembly (< {minf})")
    return {"commands": [info], "outputs": ["query.sam"],
            "validation": {**checks, "cds": n, "mapped": s["mapped_queries"], "mapped_frac": round(frac, 4),
                           "hits_per_cds": s["hits_per_query"]}}


def _proali_argv(cfg, st, wd, align):
    sd = os.path.join(cfg.outdir, "samples", st.sample, st.asm)
    rd = os.path.join(cfg.outdir, "reference", cfg.ref["id"])
    p = st.params
    argv = ["anchorwave", "proali", "-i", cfg.ref["gff3"], "-as", os.path.join(rd, "cds", "cds.fa"),
            "-r", cfg.ref["fasta"], "-a", os.path.join(sd, "query_map", "query.sam"),
            "-ar", os.path.join(rd, "ref_map", "ref.sam"), "-s", cfg.samples[st.sample]["assemblies"][st.asm]["fasta"],
            "-n", os.path.join(wd, "anchors"), "-R", str(p["R"]), "-Q", str(p["Q"]), "-e", str(p["e"]),
            "-m", str(p["min_exon"]), "-t", str(cfg.res["threads"])] + [str(x) for x in p["common_args"] + p["args"]]
    if align:
        argv += ["-o", os.path.join(wd, "alignment.maf"), "-f", os.path.join(wd, "alignment.f.maf"),
                 "-b", os.path.join(wd, "methods.bed")] + [str(x) for x in p["align_args"]]
        if cfg.mem_gib():
            argv += ["-M", str(cfg.mem_gib())]
    else:
        argv += ["-ns"]
    return argv


def body_anchors(cfg, rt, st, wd):
    info = run_cmd(cfg, rt, _proali_argv(cfg, st, wd, False), wd, "proali")
    summ = anchor_summary(parse_anchors(os.path.join(wd, "anchors")))
    if summ["gene_anchors"] == 0:
        raise StageError("no gene anchors")
    return {"commands": [info], "outputs": ["anchors"], "validation": summ}


def body_align(cfg, rt, st, wd):
    info = run_cmd(cfg, rt, _proali_argv(cfg, st, wd, True), wd, "proali")
    summ = anchor_summary(parse_anchors(os.path.join(wd, "anchors")))
    m = maf_summary(os.path.join(wd, "alignment.maf"))
    return {"commands": [info], "outputs": ["anchors", "alignment.maf", "alignment.f.maf", "methods.bed"],
            "validation": {"anchors": summ, "maf": m}}


BODIES = {"ref_cds": body_ref_cds, "ref_map": body_ref_map, "query_map": body_query_map,
          "anchors": body_anchors, "align": body_align}


# ---------------------------------------------------------------- runner

def preflight(cfg, st):
    res, est = cfg.res, st.estimate
    problems = []
    disk = free_disk_gb(cfg.outdir)
    if disk < res["reserve_disk_gb"] + est["disk_gb"]:
        problems.append(f"free disk {disk:.0f} GB < reserve {res['reserve_disk_gb']} + estimate {est['disk_gb']:.0f}")
    ram = available_ram_gb()
    if ram < res["min_free_ram_gb"] + est["ram_gb"]:
        problems.append(f"MemAvailable {ram:.0f} GB < floor {res['min_free_ram_gb']} + estimate {est['ram_gb']:.0f}")
    busy = cpu_busy_fraction(parse_cpus(res["cpus"]))
    if busy > res.get("max_cpu_busy_frac", 0.2):
        problems.append(f"CPUs {res['cpus']} are {busy:.0%} busy (> {res.get('max_cpu_busy_frac', 0.2):.0%})")
    return problems, {"free_disk_gb": round(disk, 1), "mem_available_gb": round(ram, 1),
                      "cpus_busy_frac": round(busy, 3)}


def run_stage(cfg, rt, st, versions, rerun_stale):
    state = st.state(cfg)
    if state == "validated":
        log(f"{st.sid}: validated, skipped")
        return True
    if state == "orphaned":
        raise SystemExit(f"{st.sid}: processes of an earlier attempt are still alive ({st.partial}); see status")
    if state == "invalid":
        raise SystemExit(f"{st.sid}: {st.outdir} exists without VALIDATED.json; inspect it, refusing to touch it")
    if state == "stale":
        if not rerun_stale:
            log(f"{st.sid}: STALE (inputs/parameters changed); not rerun without --rerun-stale")
            return False
        dest = f"{st.outdir}.previous-{stamp()}"
        os.rename(st.outdir, dest)
        log(f"{st.sid}: stale result kept as {dest}")
    if os.path.isdir(st.partial):
        dest = f"{st.outdir}.attempt-{stamp()}"
        os.rename(st.partial, dest)
        log(f"{st.sid}: earlier attempt kept as {dest}")
    problems, snapshot = preflight(cfg, st)
    if problems:
        log(f"{st.sid}: preflight refused: " + "; ".join(problems))
        return False
    os.makedirs(st.partial)
    rec = {"schema": SCHEMA, "stage": st.sid, "start": now(), "host": socket.gethostname(), "user": getpass.getuser(),
           "pid": os.getpid(), "run_id": cfg.run_id, "config": cfg.path, "profile": cfg.profile,
           "science_file": cfg.science_file,
           "params": st.params, "params_hash": short_hash(st.params), "code_version": CODE_VERSION,
           "pipeline_git": git_state(cfg.pipeline_dir), "inputs": file_signature(st.inputs),
           "resources": cfg.res, "memory_cap": cfg.memory_cap, "estimate": st.estimate, "preflight": snapshot,
           "runtime": versions}
    log(f"{st.sid}: running in {st.partial}")
    t0 = time.time()
    try:
        result = BODIES[st.kind](cfg, rt, st, st.partial)
    except StageError as e:
        rec.update(status="failed", end=now(), elapsed_s=round(time.time() - t0, 1), error=str(e.args[0]),
                   command=e.args[1] if len(e.args) > 1 else None)
        write_json(os.path.join(st.partial, "FAILED.json"), rec)
        log(f"{st.sid}: FAILED: {e.args[0]}")
        return False
    except (Interrupted, KeyboardInterrupt) as e:
        why = e.args[0] if e.args else "SIGINT"
        rec.update(status="interrupted", end=now(), elapsed_s=round(time.time() - t0, 1),
                   error=f"runner received {why}; stage process group stopped")
        write_json(os.path.join(st.partial, "FAILED.json"), rec)
        log(f"{st.sid}: INTERRUPTED ({why}); process group stopped, attempt kept in {st.partial}")
        raise
    rec.update(status="validated", end=now(), elapsed_s=round(time.time() - t0, 1),
               commands=result["commands"], validation=result["validation"],
               outputs=[[f, os.path.getsize(os.path.join(st.partial, f))] for f in result["outputs"]])
    write_json(os.path.join(st.partial, "VALIDATED.json"), rec)
    os.rename(st.partial, st.outdir)
    log(f"{st.sid}: validated in {rec['elapsed_s'] / 60:.1f} min: {json.dumps(result['validation'])[:400]}")
    return True


def cmd_check_env(cfg, args):
    rt = Runtime(cfg)
    ok = True
    for label, p in [("reference fasta", cfg.ref["fasta"]), ("reference gff3", cfg.ref["gff3"]),
                     ("image", rt.image)] + [(f"{s}:{a}", ad["fasta"]) for s, sd in cfg.samples.items()
                                             for a, ad in sd["assemblies"].items()]:
        e = os.path.exists(p) and (label in ("image", "reference gff3") or os.path.exists(p + ".fai"))
        print(f"{'ok ' if e else 'MISSING'} {label}: {p}")
        ok &= e
    v = rt.versions()
    print(json.dumps(v, indent=2))
    print(f"memory cap: {check_memory_cap(cfg.res) or 'none (watchdog only)'}")
    ok &= v["anchorwave_exit"] == 0 and v["minimap2_exit"] == 0
    exp = cfg.runtime.get("expect", {})
    for k in ("anchorwave", "minimap2", "image_sha256"):
        if exp.get(k) and exp[k] != v[k]:
            print(f"MISMATCH {k}: expected {exp[k]}, got {v[k]}")
            ok = False
    return 0 if ok else 1


def cmd_plan(cfg, args):
    stages = select(build_stages(cfg), args.stages, args.samples)
    for st in stages:
        e = st.estimate
        print(f"{st.state(cfg):11s} {st.sid:45s} ~{e['hours']:.1f} h  {e['ram_gb']:.0f} GB RAM  "
              f"{e['disk_gb']:.1f} GB disk  -> {st.outdir}")
    return 0


def orphan_report(stages, cfg):
    """Stages whose earlier attempt still has live processes, with what to do about them."""
    out = []
    for st in stages:
        if os.path.isdir(st.partial) and st.state(cfg) == "orphaned":
            pids, how = orphan_processes(st.partial)
            rec = read_json(os.path.join(st.partial, "process.json"))
            stop = (f"systemctl --user stop {rec['scope_unit']}.scope" if rec.get("scope_unit")
                    else f"kill -TERM -- -{rec['pgid']}")
            out.append(f"{st.sid}: live processes {pids} ({how}); started {rec.get('start')} by runner pid "
                       f"{rec.get('runner_pid')}. Inspect them (ps -o pid,etime,rss,args -p "
                       f"{','.join(map(str, pids or []))}); if they are this stage's command, stop them with: {stop}")
    return out


def cmd_status(cfg, args):
    for st in build_stages(cfg):
        line = f"{st.state(cfg):11s} {st.sid}"
        v = os.path.join(st.outdir, "VALIDATED.json")
        if os.path.exists(v):
            r = read_json(v)
            line += f"  ({r['elapsed_s'] / 60:.1f} min, end {r['end']})"
        print(line)
    for line in orphan_report(build_stages(cfg), cfg):
        print("ORPHANED " + line)
    return 0


def cmd_run(cfg, args):
    os.makedirs(cfg.outdir, exist_ok=True)
    os.makedirs(cfg.tmpdir, exist_ok=True)
    with open(os.path.join(cfg.outdir, ".runner.lock"), "w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("another g2f_anchorwave run holds the lock")
        cfg.holding_lock = True
        for s in (signal.SIGTERM, signal.SIGHUP):
            signal.signal(s, _on_signal)
        try:
            return _run_locked(cfg, args)
        except Interrupted as e:
            log(f"run stopped: {e.args[0]}")
            return 143


def write_run_record(cfg, args, versions, todo):
    """<outdir>/runs/<run_id>/: RUN.json (identity and provenance of this invocation) plus
    copies of the run configuration and the science file as they were at launch."""
    d = os.path.join(cfg.outdir, "runs", cfg.run_id)
    os.makedirs(d)
    shutil.copy2(cfg.path, os.path.join(d, "config.yaml"))
    shutil.copy2(cfg.science_file, os.path.join(d, "science.yaml"))
    rec = {"schema": "g2f-anchorwave-run/1", "run_id": cfg.run_id, "start": now(), "host": socket.gethostname(),
           "user": getpass.getuser(), "argv": sys.argv, "cwd": os.getcwd(),
           "config": cfg.path, "config_sha256": file_digest(cfg.path),
           "science_file": cfg.science_file, "science_sha256": file_digest(cfg.science_file),
           "pipeline_git": git_state(cfg.pipeline_dir), "code_version": CODE_VERSION,
           "runtime": versions, "memory_cap": cfg.memory_cap, "resources": cfg.res,
           "reference": {k: cfg.ref.get(k) for k in ("id", "fasta", "gff3", "fasta_md5", "gff3_md5")},
           "assemblies": {f"{s}:{a}": {k: ad.get(k) for k in ("fasta", "method_id", "expected_sha256",
                                                              "expected_contigs", "expected_bp")}
                          for s, sd in cfg.samples.items() for a, ad in sd["assemblies"].items()},
           "selected_stages": [s.sid for s in todo], "rerun_stale": args.rerun_stale}
    path = os.path.join(d, "RUN.json")
    write_json(path, rec)
    return path, rec


def _run_locked(cfg, args):
    every = {s.sid: s for s in build_stages(cfg)}
    orphans = orphan_report(every.values(), cfg)
    if orphans:  # never start next to (or move away the directory of) a command that is still running
        raise SystemExit("refusing to run: an earlier runner was killed and left processes running:\n  "
                         + "\n  ".join(orphans))
    cfg.memory_cap = check_memory_cap(cfg.res)
    log(f"memory cap: {cfg.memory_cap or 'none (watchdog only)'}")
    rt = Runtime(cfg)
    versions = rt.versions()
    exp = cfg.runtime.get("expect", {})
    for k in ("anchorwave", "minimap2", "image_sha256"):
        if exp.get(k) and exp[k] != versions[k]:
            raise SystemExit(f"runtime {k}: expected {exp[k]}, got {versions[k]}")
    todo = select(list(every.values()), args.stages, args.samples)
    cfg.run_id = (f"{datetime.datetime.now():%Y%m%dT%H%M%S.%f}_{socket.gethostname().split('.')[0]}"
                  f"_{os.getpid()}")  # never reused: write_run_record refuses an existing directory
    run_json, run_rec = write_run_record(cfg, args, versions, todo)
    git = run_rec["pipeline_git"]
    log(f"run {cfg.run_id}: commit {git['commit']}"
        + (f" WITH UNCOMMITTED CHANGES {git['tracked_files_modified']}" if git["tracked_files_modified"]
           or git["script_modified_or_untracked"] else " (clean)") + f"; record {run_json}")
    try:
        rc = _run_stages(cfg, args, rt, versions, every, todo)
    except BaseException as e:
        run_rec.update(end=now(), status=f"aborted: {type(e).__name__} {e}")
        write_json(run_json, run_rec)
        raise
    run_rec.update(end=now(), status="all selected stages validated" if rc == 0 else "some stages not validated",
                   exit_code=rc, stage_states={s.sid: s.state(cfg) for s in todo})
    write_json(run_json, run_rec)
    return rc


def _run_stages(cfg, args, rt, versions, every, todo):
    chosen = {s.sid for s in todo}
    ok_all = True
    for st in todo:
        unmet = [d for d in st.deps if d not in chosen and every[d].state(cfg) != "validated"]
        failed_dep = [d for d in st.deps if d in chosen and every[d].state(cfg) != "validated"]
        if unmet or failed_dep:
            log(f"{st.sid}: skipped, dependencies not validated: {unmet + failed_dep}")
            ok_all = False
            continue
        ok_all &= run_stage(cfg, rt, st, versions, args.rerun_stale)
    log(f"run finished: {'all selected stages validated' if ok_all else 'some stages not validated'}")
    return 0 if ok_all else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("check-env", "plan", "run", "status"):
        p = sub.add_parser(name)
        if name in ("plan", "run"):
            p.add_argument("--stages", type=lambda s: s.split(","), default=None)
            p.add_argument("--samples", type=lambda s: s.split(","), default=None)
        if name == "run":
            p.add_argument("--rerun-stale", action="store_true")
    args = ap.parse_args(argv)
    cfg = Config(args.config)
    return {"check-env": cmd_check_env, "plan": cmd_plan, "run": cmd_run, "status": cmd_status}[args.cmd](cfg, args)


if __name__ == "__main__":
    sys.exit(main())
