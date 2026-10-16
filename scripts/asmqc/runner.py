"""Stage execution with explicit states, provenance and safe resumption.

States of a stage (derived from files in its output directory):
  planned    no output directory
  running    <dir>.partial exists and its lock is held by a live process
  failed     <dir>.partial holds a FAILED.json (kept for diagnosis)
  completed  never persisted: a stage that ran but fails validation is failed
  validated  <dir>/VALIDATED.json exists and its input signature, science
             hash and code version still match
  stale      VALIDATED.json exists but inputs or parameters changed

A stage runs in <dir>.partial and is renamed to <dir> only after its outputs
validate, so a final directory always holds a validated result. Nothing is
deleted: an interrupted or failed attempt is moved to <dir>.attempt-<stamp>
before a new one starts, and --rerun-stale moves an old result to
<dir>.previous-<stamp>.
"""
import datetime
import fcntl
import getpass
import json
import os
import shutil
import socket
import subprocess
import time

SCHEMA = "g2f-asmqc-stage/1"
CODE_VERSION = 1  # bump when stage outputs change meaning


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def stamp():
    return datetime.datetime.now().strftime("%Y%m%dT%H%M%S")


def write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def file_signature(paths):
    """Size + mtime of every input (directories are not followed)."""
    sig = []
    for p in paths:
        try:
            st = os.stat(p)
            sig.append([p, st.st_size, int(st.st_mtime)])
        except OSError:
            sig.append([p, None, None])
    return sig


def code_digest():
    """sha256 over this package's source, so uncommitted code is identifiable."""
    import hashlib
    d = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    for f in sorted(os.listdir(d)):
        if f.endswith(".py"):
            h.update(f.encode())
            with open(os.path.join(d, f), "rb") as fh:
                h.update(fh.read())
    return h.hexdigest()[:16]


def git_commit(pipeline_dir):
    r = subprocess.run(["git", "-C", pipeline_dir, "rev-parse", "--short=12", "HEAD"],
                       capture_output=True, text=True)
    d = subprocess.run(["git", "-C", pipeline_dir, "status", "--porcelain", "--untracked-files=no"],
                       capture_output=True, text=True)
    return {"commit": r.stdout.strip() or None, "dirty": bool(d.stdout.strip()), "asmqc_code_digest": code_digest()}


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


class StageError(Exception):
    pass


class Stage:
    """One unit of work. Subclasses (or instances) provide:
    sid      unique id, e.g. "map_hifi:SAMPLE_A:unpurged"
    outdir   final output directory
    deps     list of stage ids
    inputs() list of input file paths (signature)
    params() dict of science parameters this stage uses (hashed)
    run(ctx, workdir) -> dict with "outputs" (list of file names) and
             "validation" (dict); raise StageError on failure
    estimate(ctx) -> dict(threads, ram_gb, disk_final_gb, disk_tmp_gb, hours)
    heavy    True for whole-genome alignment stages
    """
    heavy = False

    def __init__(self, sid, outdir, deps=(), run=None, inputs=None, params=None, estimate=None,
                 heavy=False, describe=""):
        self.sid, self.outdir, self.deps = sid, outdir, list(deps)
        self._run, self._inputs, self._params, self._estimate = run, inputs, params, estimate
        self.heavy, self.describe = heavy, describe

    def inputs(self):
        return self._inputs() if self._inputs else []

    def params(self):
        return self._params() if self._params else {}

    def estimate(self, ctx):
        return self._estimate(ctx) if self._estimate else {}

    def run(self, ctx, workdir):
        return self._run(ctx, workdir)

    @property
    def marker(self):
        return os.path.join(self.outdir, "VALIDATED.json")

    @property
    def partial(self):
        return self.outdir + ".partial"


def _params_hash(stage):
    from .config import _hash
    return _hash({"params": stage.params(), "code_version": CODE_VERSION})


def stage_state(stage):
    if os.path.exists(stage.marker):
        rec = read_json(stage.marker)
        cur_sig = file_signature(stage.inputs())
        if rec.get("params_hash") != _params_hash(stage):
            return "stale", "parameters or code version changed"
        if rec.get("inputs") != cur_sig:
            changed = [s[0] for s, c in zip(rec.get("inputs", []), cur_sig) if s != c]
            if len(rec.get("inputs", [])) != len(cur_sig):
                changed = ["input list changed"]
            return "stale", "inputs changed: " + ", ".join(changed[:3])
        return "validated", rec.get("end", "")
    if os.path.exists(stage.outdir):
        return "unknown", "output directory exists without VALIDATED.json"
    if os.path.exists(stage.partial):
        lock = os.path.join(stage.partial, "RUNNING.lock")
        if os.path.exists(lock):
            try:
                with open(lock) as fh:
                    fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fcntl.flock(fh, fcntl.LOCK_UN)
            except OSError:
                return "running", "lock held"
        if os.path.exists(os.path.join(stage.partial, "FAILED.json")):
            return "failed", read_json(os.path.join(stage.partial, "FAILED.json")).get("error", "")[:200]
        return "failed", "interrupted attempt (no FAILED.json)"
    return "planned", ""


class Context:
    """Everything a stage needs at run time."""

    def __init__(self, cfg, pipeline_dir, log=print, dry_run=False):
        self.cfg, self.pipeline_dir, self.log, self.dry_run = cfg, pipeline_dir, log, dry_run
        self._versions = None

    @property
    def res(self):
        return self.cfg["resources"]

    @property
    def sci(self):
        return self.cfg["science"]

    def prefix(self):
        """CPU pinning and niceness for external commands."""
        p = []
        if self.res.get("cpus"):
            p += ["taskset", "-c", str(self.res["cpus"])]
        if self.res.get("nice"):
            p += ["nice", "-n", str(int(self.res["nice"]))]
        return " ".join(p)

    def tool_versions(self):
        if self._versions is None:
            out = {}
            for name, t in self.cfg["tools"].items():
                r = subprocess.run(["bash", "-c", t["version_cmd"]], capture_output=True, text=True)
                text = (r.stdout.strip() or r.stderr.strip()).splitlines()
                out[name] = {"cmd": t["cmd"], "exit_code": r.returncode,
                             "reported": text[0].strip() if text else ""}
            self._versions = out
        return self._versions

    def sh(self, command, workdir, logname="commands.log", pin=True):
        """Run one shell command (bash, pipefail) and append it and its output to the stage log."""
        if pin and self.prefix():
            command = f"{self.prefix()} bash -o pipefail -c {_quote(command)}"
        log = os.path.join(workdir, logname)
        with open(log, "a") as fh:
            fh.write(f"\n### {now()}\n$ {command}\n")
            fh.flush()
            env = dict(os.environ, TMPDIR=self.cfg["tmpdir"])
            r = subprocess.run(["bash", "-o", "pipefail", "-c", command], stdout=fh, stderr=subprocess.STDOUT,
                               cwd=workdir, env=env)
            fh.write(f"### exit {r.returncode} at {now()}\n")
        cmds = os.path.join(workdir, "commands.sh")
        with open(cmds, "a") as fh:
            fh.write(command + "\n")
        if r.returncode != 0:
            raise StageError(f"command failed (exit {r.returncode}); see {log}")


def _quote(s):
    return "'" + s.replace("'", "'\"'\"'") + "'"


def check_resources(ctx, stage):
    """Refuse to start a stage whose estimate does not fit the free disk / RAM."""
    est = stage.estimate(ctx) or {}
    need_disk = est.get("disk_final_gb", 0) + est.get("disk_tmp_gb", 0) + ctx.res["reserve_disk_gb"]
    free = free_disk_gb(stage.outdir)
    if free < need_disk:
        return f"free disk {free:.0f} GB < estimate {need_disk:.0f} GB (incl. {ctx.res['reserve_disk_gb']} GB reserve)"
    tmp_free = free_disk_gb(ctx.cfg["tmpdir"])
    if tmp_free < est.get("disk_tmp_gb", 0) + ctx.res["reserve_disk_gb"]:
        return f"free space in tmpdir {tmp_free:.0f} GB too small"
    need_ram = max(est.get("ram_gb", 0), ctx.res["min_free_ram_gb"] if stage.heavy else 0)
    ram = available_ram_gb()
    if ram < need_ram:
        return f"available RAM {ram:.0f} GB < {need_ram:.0f} GB"
    return None


def execute(ctx, stage):
    """Run one stage. Returns the state after the attempt."""
    state, why = stage_state(stage)
    if state == "validated":
        return state
    if state in ("stale", "unknown", "running"):
        raise StageError(f"{stage.sid}: {state} ({why}); not touching it")
    problem = check_resources(ctx, stage)
    if problem:
        raise StageError(f"{stage.sid}: resource gate: {problem}")
    if os.path.exists(stage.partial):  # preserve a failed/interrupted attempt
        dest = f"{stage.outdir}.attempt-{stamp()}"
        os.rename(stage.partial, dest)
        ctx.log(f"  previous attempt kept in {dest}")
    os.makedirs(stage.partial)
    os.makedirs(ctx.cfg["tmpdir"], exist_ok=True)
    lock_path = os.path.join(stage.partial, "RUNNING.lock")
    lockfh = open(lock_path, "w")
    fcntl.flock(lockfh, fcntl.LOCK_EX)
    lockfh.write(f"{socket.gethostname()} pid {os.getpid()} {now()}\n")
    lockfh.flush()
    rec = {
        "schema": SCHEMA, "stage": stage.sid, "start": now(), "host": socket.getfqdn(),
        "user": getpass.getuser(), "pid": os.getpid(),
        "config": ctx.cfg["config_path"], "profile": ctx.cfg["profile"],
        "science_file": ctx.cfg["science_path"], "params": stage.params(),
        "params_hash": _params_hash(stage), "code_version": CODE_VERSION,
        "pipeline_git": git_commit(ctx.pipeline_dir),
        "inputs": file_signature(stage.inputs()),
        "resources": ctx.res, "estimate": stage.estimate(ctx),
    }
    t0 = time.time()
    try:
        result = stage.run(ctx, stage.partial) or {}
        rec["tools"] = ctx.tool_versions()
        for name in result.get("outputs", []):
            if not os.path.exists(os.path.join(stage.partial, name)):
                raise StageError(f"declared output missing: {name}")
        rec.update({"status": "validated", "end": now(), "elapsed_s": round(time.time() - t0, 1),
                    "outputs": [[n, os.path.getsize(os.path.join(stage.partial, n))]
                                for n in result.get("outputs", [])],
                    "validation": result.get("validation", {})})
        write_json(os.path.join(stage.partial, "VALIDATED.json"), rec)
        fcntl.flock(lockfh, fcntl.LOCK_UN)
        lockfh.close()
        os.remove(lock_path)
        os.rename(stage.partial, stage.outdir)
        return "validated"
    except BaseException as e:  # record every failure, including Ctrl-C
        rec.update({"status": "failed", "end": now(), "elapsed_s": round(time.time() - t0, 1),
                    "error": f"{type(e).__name__}: {e}"})
        write_json(os.path.join(stage.partial, "FAILED.json"), rec)
        fcntl.flock(lockfh, fcntl.LOCK_UN)
        lockfh.close()
        if isinstance(e, (StageError, KeyboardInterrupt)):
            raise
        raise StageError(f"{stage.sid}: {type(e).__name__}: {e}") from e


def move_stale(ctx, stage):
    """--rerun-stale: keep the old result under a new name, never delete it."""
    dest = f"{stage.outdir}.previous-{stamp()}"
    os.rename(stage.outdir, dest)
    ctx.log(f"  stale result of {stage.sid} kept in {dest}")


def readtext(path, errors=None):
    with open(path, errors=errors) as fh:
        return fh.read()


def copy_if_exists(src, dst):
    if src and os.path.exists(src):
        shutil.copy2(src, dst)
        return True
    return False
