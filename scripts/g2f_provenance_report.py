#!/usr/bin/env python3
"""Tabulate pipeline provenance records (schema g2f-provenance/1) found under one or more directories.

One row per record: sample, step, method_id, status, exit code, start/end, pipeline Git commit and whether the
working tree was dirty, the tool versions the record itself reports, and the record path. Records are read only;
nothing is re-executed, so the table states what was logged, not what was reproduced.

  g2f_provenance_report.py DIR [DIR ...] [--sample SAMPLE_A] [--step hifiasm] [--out table.tsv]
"""
import argparse
import csv
import json
import os
import sys

FIELDS = ["sample", "step", "method_id", "status", "exit_code", "start", "end", "pipeline_commit", "pipeline_dirty",
          "tools_reported", "record"]


def records(dirs):
    for top in dirs:
        for root, _, files in os.walk(top):
            for f in files:
                if not f.endswith(".json"):
                    continue
                p = os.path.join(root, f)
                try:
                    with open(p) as fh:
                        d = json.load(fh)
                except (OSError, ValueError):
                    continue
                if isinstance(d, dict) and d.get("schema") == "g2f-provenance/1":
                    yield p, d


def row(p, d):
    git = d.get("pipeline_git") or {}
    tools = d.get("tools") or {}
    rep = ";".join(f"{k}={(v or {}).get('reported', '')}" for k, v in sorted(tools.items()) if isinstance(v, dict))
    val = lambda x: "" if x is None else x  # noqa: E731
    return {"sample": val(d.get("sample")), "step": val(d.get("step")), "method_id": val(d.get("method_id")),
            "status": val(d.get("status")), "exit_code": val(d.get("exit_code")), "start": val(d.get("start")),
            "end": val(d.get("end")), "pipeline_commit": val(git.get("commit")), "pipeline_dirty": val(git.get("dirty")),
            "tools_reported": rep, "record": os.path.abspath(p)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--sample")
    ap.add_argument("--step")
    ap.add_argument("--out")
    a = ap.parse_args()
    rows = [row(p, d) for p, d in records(a.dirs)
            if (not a.sample or d.get("sample") == a.sample) and (not a.step or d.get("step") == a.step)]
    rows.sort(key=lambda r: (r["sample"], r["method_id"], r["start"], r["step"]))
    fh = open(a.out, "w") if a.out else sys.stdout
    w = csv.DictWriter(fh, FIELDS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    if a.out:
        fh.close()
        print(f"{len(rows)} records -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
