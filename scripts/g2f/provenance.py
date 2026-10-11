"""Execution-time provenance: run one workflow step and record what happened.

Every heavy rule runs its command through `g2f_tool.py run`, which
  1. runs each tool's version command and (optionally) checks it against the
     version required by the scientific method, refusing to run on mismatch;
  2. executes the exact command string with `bash -euo pipefail -c`;
  3. writes a JSON record (command, versions, inputs/outputs with sizes, host,
     timing, exit status, pipeline Git commit) whether the step succeeds or
     fails. The record is a Snakemake log file, so it survives failures.
"""
import datetime
import getpass
import json
import os
import socket
import subprocess
import sys
import time

SCHEMA = "g2f-provenance/1"


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _file_info(path):
    info = {"path": path, "realpath": os.path.realpath(path)}
    try:
        st = os.stat(path)
        info["size"] = st.st_size if not os.path.isdir(path) else None
        info["is_dir"] = os.path.isdir(path)
        info["mtime"] = datetime.datetime.fromtimestamp(st.st_mtime).astimezone().isoformat(
            timespec="seconds")
    except OSError:
        info["exists"] = False
    return info


def _git_state(pipeline_dir):
    def git(*args):
        r = subprocess.run(["git", "-C", pipeline_dir, *args], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    commit = git("rev-parse", "HEAD")
    dirty = git("status", "--porcelain", "--untracked-files=no")
    return {"commit": commit, "dirty": bool(dirty) if dirty is not None else None}


def tool_versions(version_cmds):
    out = {}
    for name, cmd in version_cmds.items():
        r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        text = (r.stdout.strip() or r.stderr.strip()).splitlines()
        out[name] = {"command": cmd, "exit_code": r.returncode,
                     "reported": text[0].strip() if text else ""}
    return out


def write_json(path, data):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, path)


def run_step(prov, step, command, sample=None, method=None, version_cmds=None, expect=None,
             inputs=(), outputs=(), pipeline_dir=None, extra=None):
    rec = {
        "schema": SCHEMA,
        "step": step,
        "sample": sample,
        "method_id": method,
        "command": command,
        "cwd": os.getcwd(),
        "host": socket.getfqdn(),
        "user": getpass.getuser(),
        "start": _now(),
        "environment": {k: os.environ.get(k) for k in
                        ("PATH", "LOADEDMODULES", "TMPDIR", "OMP_NUM_THREADS", "SLURM_JOB_ID")},
        "python": sys.version.split()[0],
        "pipeline_git": _git_state(pipeline_dir) if pipeline_dir else None,
        "inputs": [_file_info(p) for p in inputs],
        "extra": extra or {},
    }
    rec["tools"] = tool_versions(version_cmds or {})
    mismatches = []
    for name, want in (expect or {}).items():
        got = rec["tools"].get(name, {}).get("reported", "")
        if want not in got:
            mismatches.append(f"{name}: expected version '{want}', tool reported '{got}'")
    t0 = time.time()
    if mismatches:
        rec.update(status="failed", failure="version_mismatch", errors=mismatches, exit_code=3)
        write_json(prov, rec)
        sys.stderr.write("Refusing to run: " + "; ".join(mismatches) + "\n")
        return 3
    write_json(prov, {**rec, "status": "running"})
    rc = subprocess.run(["bash", "-euo", "pipefail", "-c", command]).returncode
    rec.update(end=_now(), wall_seconds=round(time.time() - t0, 1), exit_code=rc,
               status="completed" if rc == 0 else "failed",
               outputs=[_file_info(p) for p in outputs])
    if rc != 0:
        rec["failure"] = "nonzero_exit"
    write_json(prov, rec)
    return rc
