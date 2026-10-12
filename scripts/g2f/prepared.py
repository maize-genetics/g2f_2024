"""Verification of prepared hifiasm inputs, and raw-file release accounting.

Prepared inputs are the files hifiasm actually reads:
  * HiFi: FASTQ converted from each BAM (results/<LINE>/reads/hifi/converted/),
    or a FASTQ/FASTA input used in place;
  * ONT: reads >= the method's length threshold per run
    (results/<LINE>/reads/ont/ul<L>/), or an acknowledged pre-filtered file
    used in place.

`verify_prepared` (a workflow gate before hifiasm) checks each prepared file
in one full pass (gzip integrity, MD5, record count) against the read QC
that produced it, and reconciles the runs with the raw-data manifest.

`release_report` (run any time, cheap) decides which LOCAL raw copies
(BAM/tar) a person may delete: only when the prepared file derived from
them was verified, is still present and unchanged, and the original
remote copy was confirmed by `scripts/check_remote_copies.sh`. Nothing here
deletes anything.
"""
import csv
import datetime
import hashlib
import json
import os
import zlib

RAW_FORMATS = ("bam", "tar")


# --------------------------------------------------------------------------
# One-pass check of a gzip file
# --------------------------------------------------------------------------

def check_gz(path, chunk=1 << 24):
    """MD5 of the compressed file; full decompression (all gzip members);
    line and record (FASTQ: lines/4) counts. A truncated last member is an error."""
    md5, lines, raw = hashlib.md5(), 0, 0
    d, in_member = zlib.decompressobj(31), False
    ok, err = True, ""
    with open(path, "rb") as fh:
        while ok:
            buf = fh.read(chunk)
            if not buf:
                break
            md5.update(buf)
            data = buf
            while data:
                try:
                    out = d.decompress(data)
                except zlib.error as exc:
                    ok, err = False, f"gzip error: {exc}"
                    break
                in_member = True
                lines += out.count(b"\n")
                raw += len(out)
                if d.eof:
                    data, in_member = d.unused_data, False
                    d = zlib.decompressobj(31)
                else:
                    data = b""
    if ok and in_member:
        ok, err = False, "truncated gzip member at end of file"
    if ok and raw == 0:
        ok, err = False, "empty file"
    st = os.stat(path)
    return {"path": os.path.abspath(path), "size": st.st_size, "mtime": int(st.st_mtime),
            "md5": md5.hexdigest(), "gzip_ok": ok, "error": err,
            "decompressed_bytes": raw, "lines": lines,
            "records": lines // 4 if lines % 4 == 0 else None}


# --------------------------------------------------------------------------
# Manifest (optional)
# --------------------------------------------------------------------------

def load_manifest(path):
    """Rows of the raw-data manifest (TSV: genotype, technology, source_path, size_bytes, ...)."""
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [{k.strip(): (v or "").strip() for k, v in r.items()} for r in rows]


def _tech(t):
    return {"hifi": "hifi", "ont": "ont"}.get(t.lower(), t.lower())


# --------------------------------------------------------------------------
# Gate before hifiasm
# --------------------------------------------------------------------------

def verify_prepared(sample, inventory, hifi_val, ont_val, checks, raw_md5, hifiasm_inputs,
                    min_len=None, manifest_rows=None):
    """checks: {(tech, input_id): check_gz result}; raw_md5: {(tech, input_id): md5}."""
    errors, warnings, runs = [], [], []
    hifi_per = (hifi_val or {}).get("per_input", {})
    ont_per = (ont_val or {}).get("per_input", {})
    for rec in inventory["inputs"]:
        tech, iid, status = rec["tech"], rec["input_id"], rec["status"]
        if status == "ignored":
            continue
        run = {"tech": tech, "input_id": iid, "raw_path": rec["path"], "raw_format": rec["format"],
               "raw_size": rec["size"], "raw_fingerprint": rec["fingerprint"], "status": status,
               "notes": rec.get("notes", [])}
        if status == "excluded":
            run["role"] = "excluded by human decision (samples.yaml); not used"
            runs.append(run)
            continue
        run["raw_md5"] = raw_md5.get((tech, iid))
        c = checks.get((tech, iid))
        if c is None:
            run["role"] = "used in place: the raw file is itself the hifiasm input"
            run["prepared_path"] = rec["path"]
        else:
            run["role"] = "converted/filtered: hifiasm reads the prepared file"
            run.update(prepared_path=c["path"], prepared_size=c["size"], prepared_md5=c["md5"],
                       prepared_mtime=c["mtime"], prepared_records=c["records"],
                       gzip_ok=c["gzip_ok"])
            if not c["gzip_ok"]:
                errors.append(f"{tech}:{iid}: prepared file failed gzip check: {c['error']}")
            if c["records"] is None:
                errors.append(f"{tech}:{iid}: prepared file is not whole FASTQ records")
            if tech == "hifi":
                exp = (hifi_per.get(iid) or {}).get("n_reads")
                bam = (hifi_per.get(iid) or {}).get("bam_primary_records")
                run.update(expected_records=exp, bam_primary_records=bam,
                           bases=(hifi_per.get(iid) or {}).get("total_bases"))
                if exp is None or c["records"] != exp:
                    errors.append(f"hifi:{iid}: {c['records']} FASTQ records, read QC counted {exp}")
                if bam is not None and c["records"] != bam:
                    errors.append(f"hifi:{iid}: {c['records']} FASTQ records but BAM has {bam}")
            else:
                ret = (ont_per.get(iid) or {}).get("retained") or {}
                sup = (ont_per.get(iid) or {}).get("supplied") or {}
                run.update(expected_records=ret.get("n_reads"), bases=ret.get("total_bases"),
                           supplied_reads=sup.get("n_reads"), supplied_bases=sup.get("total_bases"),
                           min_read_len=min_len)
                if c["records"] != ret.get("n_reads"):
                    errors.append(f"ont:{iid}: {c['records']} FASTQ records, filtered-read QC "
                                  f"counted {ret.get('n_reads')}")
                if not ret.get("n_reads"):
                    errors.append(f"ont:{iid}: no reads >= {min_len} bp retained")
            if not run.get("bases"):
                errors.append(f"{tech}:{iid}: zero bases")
        runs.append(run)
    # every hifiasm input must be one of the verified prepared files
    prepared = {r.get("prepared_path") for r in runs if r.get("prepared_path")}
    missing = [p for p in hifiasm_inputs if os.path.abspath(p) not in
               {os.path.abspath(x) for x in prepared}]
    if missing:
        errors.append("hifiasm inputs not covered by verification: " + ", ".join(missing))
    # reconcile runs with the manifest (if one is configured)
    manifest = []
    if manifest_rows is not None:
        rows = [r for r in manifest_rows if r.get("genotype") == sample]
        names = {os.path.basename(r["source_path"]): r for r in rows}
        raw_runs = [r for r in runs if r["raw_format"] in RAW_FORMATS or r["status"] == "excluded"]
        local = {os.path.basename(r["raw_path"]): r for r in raw_runs}
        derived = [r for r in runs if r["raw_format"] not in RAW_FORMATS and r["status"] == "use"]
        for n, r in sorted(local.items()):
            if n not in names:
                errors.append(f"{r['tech']}:{r['input_id']}: raw file {n} is not in the manifest "
                              "(unexpected run)")
            else:
                m = names[n]
                r["manifest_source_path"] = m["source_path"]
                r["manifest_size"] = int(m.get("size_bytes") or 0)
                if _tech(m.get("technology", "")) != r["tech"]:
                    errors.append(f"{n}: manifest technology {m.get('technology')} != {r['tech']}")
                if r["manifest_size"] and r["raw_size"] and r["manifest_size"] != r["raw_size"]:
                    errors.append(f"{n}: local size {r['raw_size']} != manifest {r['manifest_size']}")
        absent = sorted(set(names) - set(local))
        if absent and derived:
            warnings.append(f"{len(absent)} manifest run(s) not present locally; this sample uses "
                            "prepared derivatives in place (" + ", ".join(sorted(
                                os.path.basename(r['raw_path']) for r in derived)) + ")")
        elif absent:
            errors.append("manifest run(s) missing locally (expected runs not included): "
                          + ", ".join(absent))
        manifest = [{"source_path": r["source_path"], "technology": r.get("technology"),
                     "size_bytes": r.get("size_bytes"),
                     "local": os.path.basename(r["source_path"]) in local} for r in rows]
    return {"ok": not errors, "errors": errors, "warnings": warnings, "sample": sample,
            "runs": runs, "hifiasm_inputs": [os.path.abspath(p) for p in hifiasm_inputs],
            "manifest_runs": manifest,
            "verified_at": datetime.datetime.now().isoformat(timespec="seconds")}


# --------------------------------------------------------------------------
# Raw-copy release report (no deletion)
# --------------------------------------------------------------------------

def load_remote_table(path):
    """TSV from scripts/check_remote_copies.sh: source_path, exists, size, mtime, checked_at."""
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="") as fh:
        return {r["source_path"]: r for r in csv.DictReader(fh, delimiter="\t")}


def release_report(verified, remote, hifi_temporary=False):
    """Per raw file: eligible for manual deletion, with reasons. `verified` is
    the output of verify_prepared (or None if not yet verified). With
    `hifi_temporary`, the converted HiFi FASTQ is deleted after hifiasm in the
    latest run, so the BAMs stay the local source even while the FASTQ exists."""
    out = []
    for r in (verified or {}).get("runs", []):
        if r["raw_format"] not in RAW_FORMATS:
            continue
        row = {"tech": r["tech"], "input_id": r["input_id"], "local_path": r["raw_path"],
               "size": r["raw_size"], "raw_md5": r.get("raw_md5"),
               "fingerprint": r["raw_fingerprint"], "status": r["status"],
               "source_path": r.get("manifest_source_path"),
               "prepared_path": r.get("prepared_path"), "reasons": []}
        if not verified.get("ok"):
            row["reasons"].append("prepared inputs not verified")
        if r["status"] == "use":
            p = r.get("prepared_path")
            if not p or p == r["raw_path"]:
                row["reasons"].append("no separate prepared file")
            elif r["tech"] == "hifi" and hifi_temporary:
                row["reasons"].append("prepared HiFi FASTQ is temporary in the latest run "
                                      "(deleted after hifiasm; raw copy is the local source)")
            elif not os.path.exists(p):
                row["reasons"].append("prepared file not retained (raw copy is the local source)")
            else:
                st = os.stat(p)
                if st.st_size != r.get("prepared_size") or int(st.st_mtime) != r.get("prepared_mtime"):
                    row["reasons"].append("prepared file changed since verification")
            if not r.get("raw_md5"):
                row["reasons"].append("no checksum of the raw file")
        if not row["source_path"]:
            row["reasons"].append("not in the manifest (no known remote copy)")
        else:
            rem = remote.get(row["source_path"])
            if not rem:
                row["reasons"].append("remote copy not checked (run scripts/check_remote_copies.sh)")
            elif rem.get("exists") != "yes" or str(rem.get("size")) != str(r.get("manifest_size")):
                row["reasons"].append(f"remote copy missing or size differs ({rem.get('size')})")
            else:
                row["remote_checked_at"] = rem.get("checked_at")
        row["local_present"] = os.path.exists(r["raw_path"])
        row["eligible"] = not row["reasons"]
        out.append(row)
    return out
