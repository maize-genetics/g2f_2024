#!/usr/bin/env python3
"""Reference-guided pseudomolecules from sequence-anchor markers with ALLMAPS
(HISAT2 marker mapping, ALLMAPS merge + path), as resumable stages with explicit
states. The contigs are never split, trimmed or edited: a pseudomolecule is the
source contigs in ALLMAPS order and orientation, joined by fixed-length N gaps.

Stages (in order):
  markers        validate the marker FASTA + metadata; one unique marker per distinct
                 sequence; full marker-ID -> unique-marker -> map-position table
  index:S        hisat2-build of assembly S, inside its own stage directory
  map:S          HISAT2 mapping of the unique markers; unique-alignment filter
  scaffold:S     per-contig marker support and LG conflicts, ALLMAPS input map,
                 ALLMAPS merge + path, independent AGP reconstruction check
  scaffold_refpos:S   optional variant (config reference_positions): the same, but each
                 marker's map position is replaced by its own unique position on its LG's
                 sequence in a reference sample of this config (validated map:REF), in Mb.
                 Use when the supplied positions are coordinates of an older reference
                 version whose local orientation differs from the current one.

States, directories and the no-overwrite rules are those of g2f_anchorwave.py:
  planned / running / failed / interrupted / orphaned / validated / stale / invalid.
A stage runs in <dir>.partial and becomes <dir> only after its outputs validate.

  g2f_allmaps.py -c CONFIG check-env
  g2f_allmaps.py -c CONFIG plan   [--stages ...] [--samples ...]
  g2f_allmaps.py -c CONFIG run    [--stages ...] [--samples ...] [--rerun-stale]
  g2f_allmaps.py -c CONFIG status

Reference-guided ordering is not independent validation: a contig placed on a
chromosome is not thereby structurally correct, and gaps between contigs have
no measured size.
"""
import argparse
import bisect
import collections
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
import time

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2f_anchorwave as aw  # noqa: E402  shared runner, watchdog, locks and orphan detection

SCHEMA = "g2f-allmaps-stage/1"
# per stage kind; bump a kind's version when its outputs change meaning (only that kind becomes stale)
#   scaffold 2: fixed map file name; merged-marker count and pseudomolecule presence checked (defensive:
#               jcvi 1.6.6 already replaces '-' in map names, the version-1 outputs were not affected)
#   scaffold 3 / scaffold_refpos 2: robust map intervals (5th-95th percentile) and the placement flags
#               OVERLAPS_LARGER (redundant or allelic copy candidate) and ORDER_AMBIGUOUS
CODE_VERSIONS = {"markers": 1, "index": 1, "map": 1, "scaffold": 3, "scaffold_refpos": 2}
KINDS = ("markers", "index", "map", "scaffold", "scaffold_refpos")
with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()
MAP_HEADER = "Scaffold ID,scaffold position,LG,genetic position"
COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")

StageError = aw.StageError
log, now, stamp, write_json, read_json = aw.log, aw.now, aw.stamp, aw.write_json, aw.read_json
short_hash, file_digest, file_signature, open_text = aw.short_hash, aw.file_digest, aw.file_signature, aw.open_text


def revcomp(s):
    return s.translate(COMP)[::-1]


def lg_sort_key(lg):
    return (0, int(lg), "") if str(lg).isdigit() else (1, 0, str(lg))


def git_state(pipeline_dir):
    rel = "scripts/g2f_allmaps.py"
    r = subprocess.run(["git", "-C", pipeline_dir, "rev-parse", "HEAD"], capture_output=True, text=True)
    d = subprocess.run(["git", "-C", pipeline_dir, "status", "--porcelain", "--untracked-files=no", "--", rel],
                       capture_output=True, text=True)
    t = subprocess.run(["git", "-C", pipeline_dir, "status", "--porcelain", "--untracked-files=no"],
                       capture_output=True, text=True)
    tracked = subprocess.run(["git", "-C", pipeline_dir, "ls-files", "--error-unmatch", rel],
                             capture_output=True, text=True).returncode == 0
    return {"commit": r.stdout.strip() or None,
            "script_modified_or_untracked": bool(d.stdout.strip()) or not tracked,
            "tracked_files_modified": [l[3:] for l in t.stdout.splitlines()],
            "script_sha256": SCRIPT_SHA256,
            "script_sha256_on_disk_now": file_digest(os.path.abspath(__file__))}


# ---------------------------------------------------------------- FASTA access

class Fai:
    """Random access to an uncompressed FASTA through its .fai (required: proves uniform line lengths)."""

    def __init__(self, path):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            raise StageError(f"missing FASTA index {fai} (samtools faidx)")
        self.rec = collections.OrderedDict()
        with open(fai) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                self.rec[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self.fh = open(path, "rb")

    def lengths(self):
        return collections.OrderedDict((k, v[0]) for k, v in self.rec.items())

    def fetch(self, name):
        n, off, lb, lw = self.rec[name]
        if n == 0:
            return ""
        nlines = (n + lb - 1) // lb
        self.fh.seek(off)
        raw = self.fh.read((nlines - 1) * lw + (n - (nlines - 1) * lb))
        s = raw.replace(b"\n", b"").replace(b"\r", b"").decode()
        if len(s) != n:
            raise StageError(f"{self.path}: record {name} read {len(s)} bp, .fai says {n}")
        return s

    def close(self):
        self.fh.close()


def fasta_records(path):
    """(name, sequence) for every record, streamed; sequence upper-case preserved as is."""
    name, buf = None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


# ---------------------------------------------------------------- config

class Config:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(path) as fh:
            c = yaml.safe_load(fh)
        self.raw = c
        for key in ("outdir", "tmpdir", "science", "runtime", "resources", "markers", "samples"):
            if key not in c:
                raise SystemExit(f"config {path}: missing '{key}'")
        self.profile = c.get("profile", "allmaps_v1")
        self.outdir = os.path.join(c["outdir"], self.profile)
        self.tmpdir = c["tmpdir"]
        sci = c["science"]
        self.science_file = sci if os.path.isabs(sci) else os.path.join(self.pipeline_dir, sci)
        with open(self.science_file) as fh:
            self.science = yaml.safe_load(fh)
        self.runtime = c["runtime"]
        self.res = c["resources"]
        self.markers = c["markers"]
        self.samples = c["samples"]
        self.memory_cap = None
        self.run_id = None
        self.holding_lock = False
        self._validate()

    def _validate(self):
        errs = []
        for k in ("fasta", "metadata"):
            if k not in self.markers:
                errs.append(f"markers.{k} missing")
        for k in ("hisat2_image", "jcvi_image"):
            if k not in self.runtime:
                errs.append(f"runtime.{k} missing")
        for k in ("threads", "cpus", "max_rss_gb", "min_free_ram_gb", "reserve_disk_gb"):
            if k not in self.res:
                errs.append(f"resources.{k} missing")
        if self.res.get("memory_cap", "none") not in ("none", "systemd-scope"):
            errs.append("resources.memory_cap must be 'none' or 'systemd-scope'")
        f = self.science.get("filter", {})
        for k in ("min_mapq", "max_nh", "max_nm"):
            if not isinstance(f.get(k), int):
                errs.append(f"science filter.{k} must be an integer")
        for s, sd in self.samples.items():
            if not re.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", s):
                errs.append(f"sample id {s!r}: letters, digits, '_', '.', '-' only")
            if "fasta" not in sd:
                errs.append(f"samples.{s}.fasta missing")
            elif sd["fasta"].endswith((".gz", ".bgz")):
                errs.append(f"samples.{s}.fasta: use an uncompressed, faidx-indexed FASTA")
        if errs:
            raise SystemExit("config errors:\n  " + "\n  ".join(errs))


# ---------------------------------------------------------------- stages

class Stage:
    def __init__(self, sid, kind, outdir, deps, inputs, params, estimate, sample=None):
        self.sid, self.kind, self.outdir, self.deps = sid, kind, outdir, deps
        self.inputs, self.params, self.estimate, self.sample = inputs, params, estimate, sample

    @property
    def partial(self):
        return self.outdir + ".partial"

    def state(self, cfg):
        v = os.path.join(self.outdir, "VALIDATED.json")
        if os.path.exists(v):
            rec = read_json(v)
            same = (rec.get("params_hash") == short_hash(self.params)
                    and rec.get("inputs") == file_signature(self.inputs)
                    and rec.get("code_version") == CODE_VERSIONS[self.kind])
            return "validated" if same else "stale"
        if os.path.isdir(self.outdir):
            return "invalid"
        if os.path.isdir(self.partial):
            if os.path.exists(os.path.join(self.partial, "FAILED.json")):
                return "failed"
            if not cfg.holding_lock and aw.run_lock_held(cfg):
                return "running"
            pids, _ = aw.orphan_processes(self.partial)
            return "orphaned" if pids is None or pids else "interrupted"
        return "planned"


def build_stages(cfg):
    sci, od = cfg.science, cfg.outdir
    mdir = os.path.join(od, "markers")
    mk_inputs = [cfg.markers["fasta"], cfg.markers["metadata"]]
    stages = [Stage("markers", "markers", mdir, [], mk_inputs, {"markers": sci["markers"]},
                    {"threads": 1, "ram_gb": 6, "disk_gb": 0.5, "hours": 0.1})]
    ufa = os.path.join(mdir, "unique_markers.fa")
    utsv = os.path.join(mdir, "unique_markers.tsv")
    for s, sd in cfg.samples.items():
        base = os.path.join(od, "samples", s)
        fa = sd["fasta"]
        gb = (sum(v[0] for v in _fai_rows(fa)) / 1e9) if os.path.exists(fa + ".fai") else 3.0
        idx = os.path.join(base, "index")
        stages.append(Stage(f"index:{s}", "index", idx, [], [fa, fa + ".fai"],
                            {"tool": "hisat2-build", "expected_sha256": sd.get("expected_sha256")},
                            {"threads": cfg.res["threads"], "ram_gb": 12 * gb, "disk_gb": 1.6 * gb, "hours": 1.5 * gb},
                            s))
        mp = os.path.join(base, "map")
        stages.append(Stage(f"map:{s}", "map", mp, ["markers", f"index:{s}"],
                            [fa, ufa, utsv, os.path.join(idx, "VALIDATED.json")],
                            {"hisat2": sci["hisat2"], "filter": sci["filter"]},
                            {"threads": cfg.res["threads"], "ram_gb": 8 * gb, "disk_gb": 0.5, "hours": 0.2}, s))
        sc = os.path.join(base, "scaffold")
        stages.append(Stage(f"scaffold:{s}", "scaffold", sc, ["markers", f"map:{s}"],
                            [fa, fa + ".fai", os.path.join(mp, "hits.tsv"), utsv],
                            {"scaffold": sci["scaffold"], "allmaps": sci["allmaps"], "lg_labels": sci["markers"]["lg_labels"],
                             "label": sd.get("label", s)},
                            {"threads": min(cfg.res["threads"], 16), "ram_gb": 4 * gb + 8, "disk_gb": 3 * gb,
                             "hours": 0.5}, s))
        rp = cfg.raw.get("reference_positions")
        if rp:
            ref = rp["sample"]
            if ref not in cfg.samples:
                raise SystemExit(f"reference_positions.sample {ref} is not a sample of this config")
            stages.append(Stage(f"scaffold_refpos:{s}", "scaffold_refpos", os.path.join(base, "scaffold_refpos"),
                                ["markers", f"map:{s}", f"map:{ref}"],
                                [fa, fa + ".fai", os.path.join(mp, "hits.tsv"), utsv,
                                 os.path.join(od, "samples", ref, "map", "hits.tsv")],
                                {"scaffold": sci["scaffold"], "allmaps": sci["allmaps"],
                                 "lg_labels": sci["markers"]["lg_labels"], "label": sd.get("label", s),
                                 "position_source": {"sample": ref, "lg_to_seq": rp["lg_to_seq"], "units": "Mb"}},
                                {"threads": min(cfg.res["threads"], 16), "ram_gb": 4 * gb + 8, "disk_gb": 3 * gb,
                                 "hours": 0.5}, s))
    return stages


def _fai_rows(fa):
    with open(fa + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            yield int(f[1]), f[0]


def select(stages, kinds_or_ids, samples):
    out = []
    for st in stages:
        if samples and st.sample is not None and st.sample not in samples:
            continue
        if kinds_or_ids and st.kind not in kinds_or_ids and st.sid not in kinds_or_ids:
            continue
        out.append(st)
    return out


# ---------------------------------------------------------------- runtime

class Image:
    """apptainer exec wrapper for one image (duck-typed for g2f_anchorwave.run_cmd)."""

    def __init__(self, cfg, image):
        self.cfg, self.image = cfg, image
        self.apptainer = cfg.runtime.get("apptainer", "apptainer")
        self.binds = cfg.runtime.get("bind", [])

    def wrap(self, argv, tmp):
        cmd = [self.apptainer, "exec", "--cleanenv", "--no-home", "--env", f"TMPDIR={tmp}",
               "--env", f"MPLCONFIGDIR={tmp}"]
        for b in self.binds:
            cmd += ["-B", f"{b}:{b}"]
        return cmd + [self.image] + argv


def versions(cfg):
    out = {}
    for key, argv in (("hisat2_image", ["hisat2", "--version"]),
                      ("jcvi_image", ["python", "-c", "import jcvi; print('jcvi', jcvi.__version__)"])):
        img = Image(cfg, cfg.runtime[key])
        out[key] = img.image
        out[key + "_sha256"] = file_digest(img.image) if os.path.exists(img.image) else None
        r = subprocess.run(img.wrap(argv, cfg.tmpdir), capture_output=True, text=True)
        lines = [l for l in (r.stdout + r.stderr).splitlines() if l.strip() and "WARNING" not in l]
        name = "hisat2" if key.startswith("hisat2") else "jcvi"
        out[name] = (lines[0].strip() if lines else None)
        out[name + "_exit"] = r.returncode
    r = subprocess.run([cfg.runtime.get("apptainer", "apptainer"), "--version"], capture_output=True, text=True)
    out["apptainer"] = r.stdout.strip()
    return out


# ---------------------------------------------------------------- stage bodies

def body_markers(cfg, rt, st, wd):
    sci = st.params["markers"]
    meta = {}
    n_meta = dup_meta = bad_meta = 0
    with open_text(cfg.markers["metadata"]) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 3:
                bad_meta += 1
                continue
            lg, mid, pos = f
            try:
                float(pos)
            except ValueError:
                bad_meta += 1
                continue
            n_meta += 1
            if mid in meta:
                dup_meta += 1
            meta[mid] = (lg, pos)
    if bad_meta or dup_meta:
        raise StageError(f"metadata: {bad_meta} malformed lines, {dup_meta} duplicate marker IDs")
    by_seq = collections.defaultdict(list)
    n_fa = 0
    lengths = collections.Counter()
    fa_ids = set()
    for name, seq in fasta_records(cfg.markers["fasta"]):
        n_fa += 1
        if name in fa_ids:
            raise StageError(f"marker FASTA: duplicate ID {name}")
        fa_ids.add(name)
        s = seq.upper()
        lengths[len(s)] += 1
        by_seq[s].append(name)
    only_fa, only_meta = fa_ids - set(meta), set(meta) - fa_ids
    if only_fa or only_meta:
        raise StageError(f"marker IDs differ: {len(only_fa)} only in FASTA, {len(only_meta)} only in metadata")
    unknown_lg = sorted({meta[m][0] for m in meta} - set(sci["lg_labels"]), key=lg_sort_key)
    if unknown_lg:
        raise StageError(f"metadata LGs without a label in science markers.lg_labels: {unknown_lg}")
    inconsistent, non_acgt, uniq = [], [], []
    for s, ids in by_seq.items():
        lgs = {meta[m][0] for m in ids}
        poss = {float(meta[m][1]) for m in ids}
        if len(lgs) > 1 or len(poss) > 1:
            inconsistent.append((s, ids))
            continue
        if set(s) - set("ACGT"):
            non_acgt.append((s, ids))
            if sci.get("drop_non_acgt", True):
                continue
        lg, pos = meta[ids[0]]
        uniq.append((lg_sort_key(lg), float(pos), s, lg, pos, sorted(ids)))
    if inconsistent and sci.get("on_inconsistent_duplicate", "fail") == "fail":
        raise StageError(f"{len(inconsistent)} sequences carry different LGs or map positions")
    uniq.sort(key=lambda t: (t[0], t[1], t[2]))
    width = max(7, len(str(len(uniq))))
    dropped = {m for _, ids in inconsistent for m in ids}
    if sci.get("drop_non_acgt", True):
        dropped |= {m for _, ids in non_acgt for m in ids}
    with open(os.path.join(wd, "unique_markers.fa"), "w") as fa, \
            open(os.path.join(wd, "unique_markers.tsv"), "w") as tsv, \
            gzip.open(os.path.join(wd, "marker_to_unique.tsv.gz"), "wt") as m2u:
        tsv.write("unique_id\tlg\tmap_pos\tn_marker_ids\tseq_sha1\tfirst_marker_id\n")
        m2u.write("marker_id\tunique_id\tlg\tmap_pos\tstatus\n")
        for i, (_, _, s, lg, pos, ids) in enumerate(uniq, 1):
            uid = f"um{i:0{width}d}"
            fa.write(f">{uid}\n{s}\n")
            tsv.write(f"{uid}\t{lg}\t{pos}\t{len(ids)}\t{hashlib.sha1(s.encode()).hexdigest()}\t{ids[0]}\n")
            for m in ids:
                m2u.write(f"{m}\t{uid}\t{lg}\t{pos}\tmapped_as_unique\n")
        for m in sorted(dropped):
            lg, pos = meta[m]
            why = "inconsistent_duplicate" if any(m in ids for _, ids in inconsistent) else "non_ACGT"
            m2u.write(f"{m}\t.\t{lg}\t{pos}\t{why}\n")
    positions = {(t[3], t[1]) for t in uniq}
    per_lg = collections.Counter(t[3] for t in uniq)
    pos_range = {}
    for t in uniq:
        lo, hi = pos_range.get(t[3], (float("inf"), float("-inf")))
        pos_range[t[3]] = (min(lo, t[1]), max(hi, t[1]))
    copies = collections.Counter(min(len(ids), 20) for ids in by_seq.values())
    val = {"marker_records": n_fa, "metadata_records": n_meta, "ids_match": True,
           "record_lengths": {str(k): v for k, v in sorted(lengths.items())},
           "distinct_sequences": len(by_seq), "redundant_records": n_fa - len(by_seq),
           "sequences_with_multiple_ids": sum(1 for v in by_seq.values() if len(v) > 1),
           "inconsistent_duplicate_sequences": len(inconsistent),
           "non_acgt_sequences": len(non_acgt), "non_acgt_dropped": bool(sci.get("drop_non_acgt", True)),
           "unique_markers_written": len(uniq), "marker_ids_represented": sum(len(t[5]) for t in uniq),
           "distinct_lg_positions": len(positions),
           "unique_markers_per_lg": {k: per_lg[k] for k in sorted(per_lg, key=lg_sort_key)},
           "map_position_range_per_lg": {k: list(pos_range[k]) for k in sorted(pos_range, key=lg_sort_key)},
           "ids_per_sequence_histogram_capped20": {str(k): v for k, v in sorted(copies.items())}}
    # independent re-read of what was written: every unique sequence once, every marker ID once
    seen_seq, n_written = set(), 0
    for name, seq in fasta_records(os.path.join(wd, "unique_markers.fa")):
        if seq in seen_seq:
            raise StageError(f"unique_markers.fa: sequence of {name} written twice")
        seen_seq.add(seq)
        n_written += 1
    with gzip.open(os.path.join(wd, "marker_to_unique.tsv.gz"), "rt") as fh:
        next(fh)
        ids_out = [l.split("\t", 1)[0] for l in fh]
    if n_written != len(uniq) or len(ids_out) != n_fa or len(set(ids_out)) != n_fa:
        raise StageError("marker tables failed re-read: counts differ from the input")
    write_json(os.path.join(wd, "marker_validation.json"), val)
    return {"commands": [], "validation": val,
            "outputs": ["unique_markers.fa", "unique_markers.tsv", "marker_to_unique.tsv.gz", "marker_validation.json"]}


def body_index(cfg, rt, st, wd):
    sd = cfg.samples[st.sample]
    fa = sd["fasta"]
    sha = file_digest(fa)
    exp = st.params.get("expected_sha256")
    if exp and exp != sha:
        raise StageError(f"{fa}: sha256 {sha} != expected {exp}")
    _f = Fai(fa)
    lengths = _f.lengths()
    _f.close()
    os.makedirs(os.path.join(wd, "idx"))
    prefix = os.path.join(wd, "idx", st.sample)
    img = Image(cfg, cfg.runtime["hisat2_image"])
    info = aw.run_cmd(cfg, img, ["hisat2-build", "-p", str(cfg.res["threads"]), fa, prefix], wd, "hisat2_build")
    # the index must describe exactly this FASTA: names and lengths from the index itself
    out = os.path.join(wd, "inspect.txt")
    # the inspect binary directly: the container's python wrapper (hisat2-inspect) needs the removed 'imp' module
    large = any(f.endswith(".ht2l") for f in os.listdir(os.path.join(wd, "idx")))
    aw.run_cmd(cfg, img, ["hisat2-inspect-l" if large else "hisat2-inspect-s", "-s", prefix], wd, "hisat2_inspect",
               stdout_file=out)
    seen = collections.OrderedDict()
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("Sequence-") and len(f) >= 3:
                seen[f[1].split()[0]] = int(f[2])
    if dict(seen) != dict(lengths):
        missing = set(lengths) - set(seen)
        diff = [k for k in lengths if k in seen and seen[k] != lengths[k]]
        raise StageError(f"index does not match the FASTA: {len(missing)} records missing, {len(diff)} lengths differ")
    parts = sorted(os.listdir(os.path.join(wd, "idx")))
    val = {"fasta": fa, "fasta_sha256": sha, "records": len(lengths), "bp": sum(lengths.values()),
           "index_prefix": os.path.join("idx", st.sample), "index_files": parts, "index_matches_fasta": True}
    write_json(os.path.join(wd, "index_identity.json"), val)
    return {"commands": [info], "validation": val,
            "outputs": [os.path.join("idx", p) for p in parts] + ["inspect.txt", "index_identity.json"]}


def classify_alignment(f, flt):
    """Classify the PRIMARY SAM record of one unique marker. Returns (class, fields)."""
    flag = int(f[1])
    if flag & 4:
        return "unmapped", None
    tags = {}
    for t in f[11:]:
        k, _, v = t.split(":", 2)
        tags[k] = v
    nh = int(tags.get("NH", 1))
    mapq = int(f[4])
    cigar = f[5]
    if "N" in cigar:
        return "spliced", None
    if "S" in cigar or "H" in cigar:
        return "clipped", None
    if nh > flt["max_nh"]:
        return "multi_nh", None
    if mapq < flt["min_mapq"]:
        return "low_mapq", None
    if flt.get("reject_zs", True) and "ZS" in tags:
        return "second_best", None
    nm = int(tags.get("NM", 0))
    if nm > flt["max_nm"]:
        return "divergent", None
    return "kept", {"contig": f[2], "pos": int(f[3]), "strand": "-" if flag & 16 else "+", "nm": nm,
                    "mapq": mapq, "nh": nh}


def body_map(cfg, rt, st, wd):
    sd = cfg.samples[st.sample]
    mdir = os.path.join(cfg.outdir, "markers")
    idx_dir = os.path.join(cfg.outdir, "samples", st.sample, "index")
    idx_rec = read_json(os.path.join(idx_dir, "VALIDATED.json"))
    prefix = os.path.join(idx_dir, "idx", st.sample)  # this sample's own index, never another stage's
    if idx_rec["validation"]["fasta_sha256"] != file_digest(sd["fasta"]):
        raise StageError("index was built from a FASTA whose sha256 differs from the configured FASTA")
    img = Image(cfg, cfg.runtime["hisat2_image"])
    sam = os.path.join(wd, "mapped.sam")
    argv = ["hisat2", "-p", str(cfg.res["threads"]), "-f", "-x", prefix, "-U",
            os.path.join(mdir, "unique_markers.fa")] + list(st.params["hisat2"]["args"]) + ["-S", sam]
    info = aw.run_cmd(cfg, img, argv, wd, "hisat2")
    meta = {}
    with open(os.path.join(mdir, "unique_markers.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            meta[f[0]] = (f[1], f[2])
    flt = st.params["filter"]
    cls = collections.Counter()
    nm_hist = collections.Counter()
    tag_counts = collections.Counter()
    n_records = n_secondary = n_supp = 0
    primary_seen = set()
    hits = []
    with open(sam) as fh, gzip.open(sam + ".gz", "wt") as gz:
        for line in fh:
            gz.write(line)
            if line.startswith("@"):
                continue
            n_records += 1
            f = line.rstrip("\n").split("\t")
            for t in f[11:]:
                k = t[:2]
                if k in ("XA", "SA", "ZS", "NH"):
                    tag_counts[k] += 1
            flag = int(f[1])
            if flag & 256:
                n_secondary += 1
                continue
            if flag & 2048:
                n_supp += 1
                continue
            if f[0] in primary_seen:
                raise StageError(f"two primary records for {f[0]}")
            primary_seen.add(f[0])
            c, h = classify_alignment(f, flt)
            cls[c] += 1
            for t in f[11:]:
                if t.startswith("NM:i:"):
                    nm_hist[min(int(t[5:]), 21)] += 1
            if h:
                lg, pos = meta[f[0]]
                hits.append((h["contig"], h["pos"], f[0], h["strand"], h["nm"], h["mapq"], lg, pos))
    os.unlink(sam)  # the identical gzip copy is kept (workflow intermediate replaced by its compressed form)
    if len(primary_seen) != len(meta):
        raise StageError(f"SAM has primary records for {len(primary_seen)} of {len(meta)} unique markers")
    hits.sort(key=lambda t: (t[0], t[1], t[2]))
    with open(os.path.join(wd, "hits.tsv"), "w") as out:
        out.write("contig\tpos\tunique_id\tstrand\tnm\tmapq\tlg\tmap_pos\n")
        for h in hits:
            out.write("\t".join(map(str, h)) + "\n")
    val = {"unique_markers": len(meta), "sam_records": n_records, "secondary_records": n_secondary,
           "supplementary_records": n_supp, "primary_classes": dict(cls), "kept": cls["kept"],
           "kept_frac": round(cls["kept"] / len(meta), 5), "tag_records": dict(tag_counts),
           "xa_sa_tags_present": bool(tag_counts["XA"] or tag_counts["SA"]),
           "nm_histogram_capped21": {str(k): v for k, v in sorted(nm_hist.items())},
           "contigs_with_kept_markers": len({h[0] for h in hits})}
    write_json(os.path.join(wd, "filter_stats.json"), val)
    return {"commands": [info], "validation": val,
            "outputs": ["mapped.sam.gz", "hits.tsv", "filter_stats.json", "hisat2.log"]}


# ---------------------------------------------------------------- scaffold helpers

def spearman(xs, ys):
    """Spearman rank correlation (average ranks for ties); None if undefined."""
    n = len(xs)
    if n < 2:
        return None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / (sxx * syy) ** 0.5


def quantile(sorted_vals, q):
    """Nearest-rank quantile of a sorted list; None if empty."""
    if not sorted_vals:
        return None
    return sorted_vals[min(len(sorted_vals) - 1, max(0, int(round(q * (len(sorted_vals) - 1)))))]


def overlap_flags(placed_rows, min_frac=0.5):
    """placed_rows: dicts with contig, object, order, length, q05, q95. Returns {contig: set(flags)}:
    OVERLAPS_LARGER  >= min_frac of the contig's robust map interval lies inside intervals of LONGER contigs
                     placed on the same object (a redundant or allelic copy candidate: the pseudomolecule would
                     hold that region twice)
    ORDER_AMBIGUOUS  the robust intervals of two adjacent contigs overlap by > min_frac of the narrower one"""
    out = collections.defaultdict(set)
    by = collections.defaultdict(list)
    for r in placed_rows:
        if r["q05"] is not None:
            by[r["object"]].append(r)
    for obj, rs in by.items():
        for r in rs:
            w = r["q95"] - r["q05"]
            if w <= 0:
                continue
            segs = sorted((max(r["q05"], o["q05"]), min(r["q95"], o["q95"])) for o in rs
                          if o["length"] > r["length"] and min(r["q95"], o["q95"]) > max(r["q05"], o["q05"]))
            cov, end = 0.0, float("-inf")
            for a, b in segs:  # union length of the overlaps
                a = max(a, end)
                if b > a:
                    cov += b - a
                    end = b
            if cov / w >= min_frac:
                out[r["contig"]].add("OVERLAPS_LARGER")
        rs = sorted(rs, key=lambda r: r["order"])
        for a, b in zip(rs, rs[1:]):
            ov = min(a["q95"], b["q95"]) - max(a["q05"], b["q05"])
            narrow = min(a["q95"] - a["q05"], b["q95"] - b["q05"])
            if ov > 0 and (narrow <= 0 or ov / narrow > min_frac):
                out[a["contig"]].add("ORDER_AMBIGUOUS")
                out[b["contig"]].add("ORDER_AMBIGUOUS")
    return out


def even_sample(items, k):
    """Deterministic, evenly spaced subset of a position-sorted list (all of it if len <= k)."""
    n = len(items)
    if n <= k:
        return list(items)
    return [items[min(n - 1, int((i + 0.5) * n / k))] for i in range(k)]


def conflict_blocks(markers, majority, min_n, min_span):
    """Runs of consecutive (contig-ordered) markers of one non-majority LG."""
    out, i = [], 0
    while i < len(markers):
        lg = markers[i][1]
        j = i
        while j + 1 < len(markers) and markers[j + 1][1] == lg:
            j += 1
        if lg != majority:
            n, span = j - i + 1, markers[j][0] - markers[i][0] + 1
            if n >= min_n and span >= min_span:
                out.append({"lg": lg, "start": markers[i][0], "end": markers[j][0], "n_markers": n,
                            "map_pos_min": min(m[2] for m in markers[i:j + 1]),
                            "map_pos_max": max(m[2] for m in markers[i:j + 1])})
        i = j + 1
    return out


def contig_support(lengths, hits, p):
    """Per-contig marker support, LG composition, conflict blocks and ALLMAPS eligibility."""
    by = collections.defaultdict(list)
    for h in hits:
        by[h["contig"]].append((h["pos"], h["lg"], h["map_pos"], h["uid"]))
    rows, selected = {}, {}
    for ctg, ln in lengths.items():
        ms = sorted(by.get(ctg, []))
        lgc = collections.Counter(m[1] for m in ms)
        majority, maj_n = (max(lgc.items(), key=lambda kv: (kv[1], -lg_sort_key(kv[0])[1])) if lgc else (None, 0))
        maj = [m for m in ms if m[1] == majority]
        distinct = len({m[2] for m in maj})
        other = len(ms) - maj_n
        blocks = conflict_blocks(ms, majority, p["conflict_min_markers"], p["conflict_min_span_bp"]) if ms else []
        flags = []
        if other >= p["mixed_lg_min_markers"] and other / max(1, len(ms)) >= p["mixed_lg_min_frac"]:
            flags.append("MIXED_LG")
        if blocks:
            flags.append("CONFLICT_BLOCK")
        if not ms:
            status = "NO_MARKERS"
        elif ln < p["min_contig_bp"]:
            status = "TOO_SHORT"
        elif maj_n < p["min_unique_markers"] or distinct < p["min_distinct_positions"]:
            status = "TOO_FEW_MARKERS"
        else:
            status = "ELIGIBLE"
        rho = spearman([m[0] for m in maj], [m[2] for m in maj]) if maj else None
        rows[ctg] = {"contig": ctg, "length": ln, "unique_markers": len(ms), "lgs": len(lgc),
                     "majority_lg": majority, "majority_markers": maj_n,
                     "majority_frac": round(maj_n / len(ms), 4) if ms else None,
                     "other_lg_markers": other,
                     "other_lgs": ";".join(f"{k}:{v}" for k, v in sorted(lgc.items(), key=lambda kv: lg_sort_key(kv[0]))
                                           if k != majority),
                     "distinct_map_positions": distinct,
                     "map_pos_min": min((m[2] for m in maj), default=None),
                     "map_pos_max": max((m[2] for m in maj), default=None),
                     "map_q05": quantile(sorted(m[2] for m in maj), 0.05),
                     "map_q95": quantile(sorted(m[2] for m in maj), 0.95),
                     "spearman_rho": None if rho is None else round(rho, 4),
                     "status": status, "flags": ",".join(flags), "conflict_blocks": blocks}
        if status == "ELIGIBLE":
            selected[ctg] = even_sample(maj, p["max_markers_per_contig"])
    return rows, selected


def write_map_csv(path, selected, labels):
    n = 0
    with open(path, "w") as fh:
        fh.write(MAP_HEADER + "\n")
        for ctg in sorted(selected):
            for pos, lg, mpos, uid in selected[ctg]:
                fh.write(f"{ctg},{pos},{labels[lg]},{mpos}\n")
                n += 1
    return n


def parse_agp(path):
    objs = collections.OrderedDict()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            o = objs.setdefault(f[0], [])
            rec = {"obj_beg": int(f[1]), "obj_end": int(f[2]), "part": int(f[3]), "type": f[4]}
            if f[4] in ("N", "U"):
                rec["gap"] = int(f[5])
            else:
                rec.update(comp=f[5], beg=int(f[6]), end=int(f[7]), ori=f[8])
            o.append(rec)
    return objs


def verify_agp(agp_path, out_fasta, src):
    """Independent reconstruction: rebuild every object from the source FASTA through the AGP and
    compare with the output FASTA record by record. Checks that every source contig is used exactly
    once, entirely, and that object coordinates are contiguous."""
    lengths = src.lengths()
    objs = parse_agp(agp_path)
    used = collections.Counter()
    problems, built = [], {}
    gap_n = gap_bp = 0
    for name, parts in objs.items():
        pos, h, total = 1, hashlib.sha256(), 0
        for k, r in enumerate(parts, 1):
            if r["part"] != k or r["obj_beg"] != pos:
                problems.append(f"{name}: part {r['part']} at {r['obj_beg']} (expected part {k} at {pos})")
            if "gap" in r:
                seg = "N" * r["gap"]
                gap_n += 1
                gap_bp += r["gap"]
            else:
                c = r["comp"]
                if c not in lengths:
                    problems.append(f"{name}: unknown component {c}")
                    continue
                if (r["beg"], r["end"]) != (1, lengths[c]):
                    problems.append(f"{name}: {c} used {r['beg']}-{r['end']} of {lengths[c]}")
                used[c] += 1
                seg = src.fetch(c)
                if r["ori"] == "-":
                    seg = revcomp(seg)
            if r["obj_end"] - r["obj_beg"] + 1 != len(seg):
                problems.append(f"{name}: part {k} spans {r['obj_end'] - r['obj_beg'] + 1} bp, content {len(seg)}")
            h.update(seg.encode())
            total += len(seg)
            pos += len(seg)
        built[name] = (h.hexdigest(), total)
    missing = [c for c in lengths if used[c] == 0]
    twice = [c for c, n in used.items() if n > 1]
    got = {}
    for name, seq in fasta_records(out_fasta):
        got[name] = (hashlib.sha256(seq.encode()).hexdigest(), len(seq))
    if set(got) != set(built):
        problems.append(f"output FASTA records != AGP objects ({len(set(got) ^ set(built))} differ)")
    mismatch = [n for n in built if n in got and got[n] != built[n]]
    src_bp = sum(lengths.values())
    out_bp = sum(v[1] for v in got.values())
    return {"objects": len(objs), "source_contigs": len(lengths), "contigs_used_once": len(lengths) - len(missing),
            "contigs_missing": missing[:20], "n_missing": len(missing), "contigs_used_twice": twice[:20],
            "n_used_twice": len(twice), "records_mismatch": mismatch[:20], "n_records_mismatch": len(mismatch),
            "source_bp": src_bp, "output_bp": out_bp, "gaps": gap_n, "gap_bp": gap_bp,
            "bp_balance_ok": out_bp == src_bp + gap_bp, "problems": problems[:50], "n_problems": len(problems),
            "ok": not (missing or twice or mismatch or problems) and out_bp == src_bp + gap_bp}, objs


def body_scaffold(cfg, rt, st, wd):
    sd = cfg.samples[st.sample]
    p = st.params["scaffold"]
    labels = st.params["lg_labels"]
    src_fa = sd["fasta"]
    src = Fai(src_fa)
    lengths = src.lengths()
    hits = []
    with open(os.path.join(cfg.outdir, "samples", st.sample, "map", "hits.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            hits.append({"contig": f[0], "pos": int(f[1]), "uid": f[2], "lg": f[6], "map_pos": float(f[7])})
    remap = None
    if st.kind == "scaffold_refpos":
        hits, remap = reference_positions(cfg, st, hits)
    unknown = {h["contig"] for h in hits} - set(lengths)
    if unknown:
        raise StageError(f"hits on {len(unknown)} contigs absent from the FASTA index (wrong index?)")
    rows, selected = contig_support(lengths, hits, p)
    # ALLMAPS inputs live in this stage directory; the source FASTA is reached through a symlink so that
    # nothing (e.g. a .sizes file) is ever written next to the source
    label = st.params["label"]
    link = os.path.join(wd, "source_contigs.fasta")
    os.symlink(os.path.realpath(src_fa), link)
    # ALLMAPS names the map after the CSV file and splits "<map>-<LG>:<pos>" on the first '-'. jcvi 1.6.6
    # replaces '-' in map names (filename_to_mapname); a fixed name keeps this independent of the jcvi version
    csv = os.path.join(wd, "pgmap.csv")
    n_used = write_map_csv(csv, selected, labels)
    if not selected:
        raise StageError("no contig is eligible for ALLMAPS (see contig support thresholds)")
    img = Image(cfg, cfg.runtime["jcvi_image"])
    bed = os.path.join(wd, f"{label}.allmaps.bed")
    cmds = [aw.run_cmd(cfg, img, ["python", "-m", "jcvi.assembly.allmaps", "merge", csv, "-o", bed], wd,
                       "allmaps_merge")]
    with open(bed) as fh:  # allmaps merge silently skips lines it cannot parse
        n_bed = sum(1 for l in fh if l.strip())
    if n_bed != n_used:
        raise StageError(f"allmaps merge wrote {n_bed} markers, the map CSV has {n_used}")
    pa = list(st.params["allmaps"]["path_args"]) + ["--cpus", str(st.estimate["threads"])]
    fmt = st.params["allmaps"].get("plot_format")
    pa += ["--format", fmt] if fmt else ["--noplot"]
    cmds.append(aw.run_cmd(cfg, img, ["python", "-m", "jcvi.assembly.allmaps", "path"] + pa + [bed, link], wd,
                           "allmaps_path"))
    pre = os.path.join(wd, f"{label}.allmaps")
    agp, fasta = pre + ".agp", pre + ".fasta"
    for f in (agp, fasta, pre + ".chr.agp", pre + ".unplaced.agp"):
        if not os.path.exists(f):
            raise StageError(f"ALLMAPS output missing: {f}")
    ver, objs = verify_agp(agp, fasta, src)
    src.close()
    if not ver["ok"]:
        write_json(os.path.join(wd, "agp_verification.json"), ver)
        raise StageError(f"AGP reconstruction failed: {ver['problems'][:3]} missing={ver['n_missing']} "
                         f"twice={ver['n_used_twice']} mismatch={ver['n_records_mismatch']}")
    write_json(os.path.join(wd, "agp_verification.json"), ver)
    pseudo = {"chr" + v for v in labels.values()}
    expected = {"chr" + labels[rows[c]["majority_lg"]] for c in selected}
    if not expected <= set(objs):
        raise StageError(f"ALLMAPS objects {sorted(set(objs) & pseudo)} lack expected pseudomolecules "
                         f"{sorted(expected - set(objs))}")
    # placement table: one row per source contig
    placed = {}
    chr_rows = []
    for name, parts in objs.items():
        if name not in pseudo:
            continue
        comps = [r for r in parts if "comp" in r]
        gaps = [r for r in parts if "gap" in r]
        chr_rows.append({"object": name, "length": parts[-1]["obj_end"], "contigs": len(comps), "gaps": len(gaps),
                         "gap_bp": sum(r["gap"] for r in gaps),
                         "contig_bp": sum(r["end"] - r["beg"] + 1 for r in comps),
                         "unoriented": sum(1 for r in comps if r["ori"] not in "+-")})
        for k, r in enumerate(comps, 1):
            placed[r["comp"]] = (name, k, r["obj_beg"], r["obj_end"], r["ori"])
    ovf = overlap_flags([{"contig": c, "object": pl[0], "order": pl[1], "length": lengths[c],
                          "q05": rows[c]["map_q05"], "q95": rows[c]["map_q95"]} for c, pl in placed.items()])
    for c, fl in ovf.items():
        rows[c]["flags"] = ",".join([f for f in rows[c]["flags"].split(",") if f] + sorted(fl))
    with open(os.path.join(wd, "placement.tsv"), "w") as out:
        cols = ["contig", "length", "status", "object", "order", "object_start", "object_end", "orientation",
                "orientation_support", "majority_lg", "lg_label_agrees", "unique_markers", "majority_markers",
                "majority_frac", "other_lg_markers", "other_lgs", "distinct_map_positions", "allmaps_markers",
                "map_pos_min", "map_pos_max", "map_q05", "map_q95", "spearman_rho", "flags"]
        out.write("\t".join(cols) + "\n")
        for ctg in lengths:
            r = rows[ctg]
            pl = placed.get(ctg)
            if pl:
                status = "PLACED_ORIENTED" if pl[4] in "+-" else "PLACED_UNORIENTED"
                rho = r["spearman_rho"]
                if (rho is None or abs(rho) < p["orientation_min_abs_rho"]
                        or r["distinct_map_positions"] < p["orientation_min_positions"]):
                    sup = "weak"
                else:
                    sup = "agrees" if (rho > 0) == (pl[4] == "+") else "DISAGREES"
                lg_ok = (pl[0] == "chr" + labels[r["majority_lg"]]) if r["majority_lg"] else False
            else:
                status = "UNPLACED_" + (r["status"] if r["status"] != "ELIGIBLE" else "BY_ALLMAPS")
                sup, lg_ok = "", ""
            vals = [ctg, r["length"], status, pl[0] if pl else "", pl[1] if pl else "", pl[2] if pl else "",
                    pl[3] if pl else "", pl[4] if pl else "", sup, r["majority_lg"] or "", lg_ok, r["unique_markers"],
                    r["majority_markers"], r["majority_frac"] if r["majority_frac"] is not None else "",
                    r["other_lg_markers"], r["other_lgs"], r["distinct_map_positions"], len(selected.get(ctg, [])),
                    "" if r["map_pos_min"] is None else r["map_pos_min"],
                    "" if r["map_pos_max"] is None else r["map_pos_max"],
                    "" if r["map_q05"] is None else r["map_q05"], "" if r["map_q95"] is None else r["map_q95"],
                    "" if r["spearman_rho"] is None else r["spearman_rho"], r["flags"]]
            out.write("\t".join(map(str, vals)) + "\n")
    with open(os.path.join(wd, "conflict_blocks.tsv"), "w") as out:
        out.write("contig\tcontig_length\tmajority_lg\tblock_lg\tstart\tend\tspan_bp\tn_markers\tmap_pos_min\tmap_pos_max\n")
        for ctg in lengths:
            for b in rows[ctg]["conflict_blocks"]:
                out.write(f"{ctg}\t{lengths[ctg]}\t{rows[ctg]['majority_lg']}\t{b['lg']}\t{b['start']}\t{b['end']}\t"
                          f"{b['end'] - b['start'] + 1}\t{b['n_markers']}\t{b['map_pos_min']}\t{b['map_pos_max']}\n")
    with open(os.path.join(wd, "pseudomolecules.tsv"), "w") as out:
        cols = ["object", "length", "contigs", "gaps", "gap_bp", "contig_bp", "unoriented"]
        out.write("\t".join(cols) + "\n")
        for r in sorted(chr_rows, key=lambda r: r["object"]):
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    st_counts = collections.Counter()
    st_bp = collections.Counter()
    with open(os.path.join(wd, "placement.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            st_counts[f[2]] += 1
            st_bp[f[2]] += int(f[1])
    total = sum(lengths.values())
    placed_bp = st_bp["PLACED_ORIENTED"] + st_bp["PLACED_UNORIENTED"]
    val = {"source_contigs": len(lengths), "source_bp": total, "output_bp": ver["output_bp"],
           "pseudomolecules": len(chr_rows), "placed_contigs": st_counts["PLACED_ORIENTED"] + st_counts["PLACED_UNORIENTED"],
           "oriented_contigs": st_counts["PLACED_ORIENTED"], "placed_bp": placed_bp,
           "placed_frac": round(placed_bp / total, 5), "oriented_bp": st_bp["PLACED_ORIENTED"],
           "status_counts": dict(st_counts), "status_bp": dict(st_bp),
           "allmaps_input_markers": n_used, "eligible_contigs": len(selected),
           "eligible_not_placed": st_counts["UNPLACED_BY_ALLMAPS"],
           "mixed_lg_contigs": sum(1 for r in rows.values() if "MIXED_LG" in r["flags"]),
           "overlaps_larger_contigs": sum(1 for c in placed if "OVERLAPS_LARGER" in rows[c]["flags"]),
           "overlaps_larger_bp": sum(lengths[c] for c in placed if "OVERLAPS_LARGER" in rows[c]["flags"]),
           "order_ambiguous_contigs": sum(1 for c in placed if "ORDER_AMBIGUOUS" in rows[c]["flags"]),
           "conflict_block_contigs": sum(1 for r in rows.values() if r["conflict_blocks"]),
           "gaps": ver["gaps"], "gap_bp": ver["gap_bp"], "agp_reconstruction_ok": ver["ok"]}
    if remap:
        val["position_source"] = remap
    write_json(os.path.join(wd, "scaffold_stats.json"), val)
    outs = sorted(f for f in os.listdir(wd) if os.path.isfile(os.path.join(wd, f)) and not os.path.islink(os.path.join(wd, f)))
    return {"commands": cmds, "validation": val, "outputs": outs}


def reference_positions(cfg, st, hits):
    """Replace each hit's map position by the marker's unique position (Mb) in the reference sample, kept only
    if the marker is unique there and lies on the sequence named for its LG."""
    ps = st.params["position_source"]
    ref = {}
    with open(os.path.join(cfg.outdir, "samples", ps["sample"], "map", "hits.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ref[f[2]] = (f[0], int(f[1]), f[6])
    out, why = [], collections.Counter()
    for h in hits:
        r = ref.get(h["uid"])
        if r is None:
            why["not_unique_in_reference"] += 1
            continue
        if r[0] != ps["lg_to_seq"].get(h["lg"]):
            why["reference_sequence_not_its_lg"] += 1
            continue
        out.append(dict(h, map_pos=round(r[1] / 1e6, 6)))
    rec = {"reference_sample": ps["sample"], "units": "Mb", "hits_in": len(hits), "hits_kept": len(out),
           "dropped": dict(why)}
    return out, rec


BODIES = {"markers": body_markers, "index": body_index, "map": body_map, "scaffold": body_scaffold,
          "scaffold_refpos": body_scaffold}


# ---------------------------------------------------------------- execution

def run_stage(cfg, st, vers, rerun_stale):
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
    problems, snapshot = aw.preflight(cfg, st)
    if problems:
        log(f"{st.sid}: preflight refused: " + "; ".join(problems))
        return False
    os.makedirs(st.partial)
    rec = {"schema": SCHEMA, "stage": st.sid, "start": now(), "host": socket.gethostname(), "user": getpass.getuser(),
           "pid": os.getpid(), "run_id": cfg.run_id, "config": cfg.path, "profile": cfg.profile,
           "science_file": cfg.science_file, "params": st.params, "params_hash": short_hash(st.params),
           "code_version": CODE_VERSIONS[st.kind], "pipeline_git": git_state(cfg.pipeline_dir),
           "inputs": file_signature(st.inputs), "resources": cfg.res, "memory_cap": cfg.memory_cap,
           "estimate": st.estimate, "preflight": snapshot, "runtime": vers}
    log(f"{st.sid}: running in {st.partial}")
    t0 = time.time()
    try:
        result = BODIES[st.kind](cfg, None, st, st.partial)
    except StageError as e:
        rec.update(status="failed", end=now(), elapsed_s=round(time.time() - t0, 1), error=str(e.args[0]),
                   command=e.args[1] if len(e.args) > 1 else None)
        write_json(os.path.join(st.partial, "FAILED.json"), rec)
        log(f"{st.sid}: FAILED: {e.args[0]}")
        return False
    except (aw.Interrupted, KeyboardInterrupt) as e:
        why = e.args[0] if e.args else "SIGINT"
        rec.update(status="interrupted", end=now(), elapsed_s=round(time.time() - t0, 1),
                   error=f"runner received {why}; stage process group stopped")
        write_json(os.path.join(st.partial, "FAILED.json"), rec)
        log(f"{st.sid}: INTERRUPTED ({why}); attempt kept in {st.partial}")
        raise
    rec.update(status="validated", end=now(), elapsed_s=round(time.time() - t0, 1),
               commands=result["commands"], validation=result["validation"],
               outputs=[[f, os.path.getsize(os.path.join(st.partial, f)),
                         file_digest(os.path.join(st.partial, f))] for f in result["outputs"]])
    write_json(os.path.join(st.partial, "VALIDATED.json"), rec)
    os.rename(st.partial, st.outdir)
    log(f"{st.sid}: validated in {rec['elapsed_s'] / 60:.1f} min: {json.dumps(result['validation'])[:400]}")
    return True


def cmd_check_env(cfg, args):
    ok = True
    for label, p in ([("marker fasta", cfg.markers["fasta"]), ("marker metadata", cfg.markers["metadata"]),
                      ("hisat2 image", cfg.runtime["hisat2_image"]), ("jcvi image", cfg.runtime["jcvi_image"])]
                     + [(s, sd["fasta"]) for s, sd in cfg.samples.items()]):
        e = os.path.exists(p) and ("image" in label or "marker" in label or os.path.exists(p + ".fai"))
        print(f"{'ok ' if e else 'MISSING'} {label}: {p}")
        ok &= e
    v = versions(cfg)
    print(json.dumps(v, indent=2))
    ok &= v["hisat2_exit"] == 0 and v["jcvi_exit"] == 0
    exp = cfg.runtime.get("expect", {})
    for k, val in exp.items():
        if val and v.get(k) != val:
            print(f"MISMATCH {k}: expected {val}, got {v.get(k)}")
            ok = False
    print(f"memory cap: {aw.check_memory_cap(cfg.res) or 'none (watchdog only)'}")
    return 0 if ok else 1


def cmd_plan(cfg, args):
    for st in select(build_stages(cfg), args.stages, args.samples):
        e = st.estimate
        print(f"{st.state(cfg):11s} {st.sid:30s} ~{e['hours']:.1f} h  {e['ram_gb']:.0f} GB RAM  "
              f"{e['disk_gb']:.1f} GB disk  -> {st.outdir}")
    return 0


def cmd_status(cfg, args):
    stages = build_stages(cfg)
    for st in stages:
        line = f"{st.state(cfg):11s} {st.sid}"
        v = os.path.join(st.outdir, "VALIDATED.json")
        if os.path.exists(v):
            r = read_json(v)
            line += f"  ({r['elapsed_s'] / 60:.1f} min, end {r['end']})"
        print(line)
    for line in aw.orphan_report(stages, cfg):
        print("ORPHANED " + line)
    return 0


def cmd_run(cfg, args):
    os.makedirs(cfg.outdir, exist_ok=True)
    os.makedirs(cfg.tmpdir, exist_ok=True)
    with open(os.path.join(cfg.outdir, ".runner.lock"), "w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("another g2f_allmaps run holds the lock")
        cfg.holding_lock = True
        for s in (signal.SIGTERM, signal.SIGHUP):
            signal.signal(s, aw._on_signal)
        try:
            return _run_locked(cfg, args)
        except aw.Interrupted as e:
            log(f"run stopped: {e.args[0]}")
            return 143


def _run_locked(cfg, args):
    every = {s.sid: s for s in build_stages(cfg)}
    orphans = aw.orphan_report(every.values(), cfg)
    if orphans:
        raise SystemExit("refusing to run: an earlier runner left processes running:\n  " + "\n  ".join(orphans))
    cfg.memory_cap = aw.check_memory_cap(cfg.res)
    vers = versions(cfg)
    for k, val in cfg.runtime.get("expect", {}).items():
        if val and vers.get(k) != val:
            raise SystemExit(f"runtime {k}: expected {val}, got {vers.get(k)}")
    todo = select(list(every.values()), args.stages, args.samples)
    cfg.run_id = f"{datetime.datetime.now():%Y%m%dT%H%M%S.%f}_{socket.gethostname().split('.')[0]}_{os.getpid()}"
    d = os.path.join(cfg.outdir, "runs", cfg.run_id)
    os.makedirs(d)
    shutil.copy2(cfg.path, os.path.join(d, "config.yaml"))
    shutil.copy2(cfg.science_file, os.path.join(d, "science.yaml"))
    run_rec = {"schema": "g2f-allmaps-run/1", "run_id": cfg.run_id, "start": now(), "host": socket.gethostname(),
               "user": getpass.getuser(), "argv": sys.argv, "cwd": os.getcwd(), "config": cfg.path,
               "config_sha256": file_digest(cfg.path), "science_sha256": file_digest(cfg.science_file),
               "pipeline_git": git_state(cfg.pipeline_dir), "code_versions": CODE_VERSIONS, "runtime": vers,
               "memory_cap": cfg.memory_cap, "resources": cfg.res, "markers": cfg.markers,
               "samples": cfg.samples, "selected_stages": [s.sid for s in todo], "rerun_stale": args.rerun_stale}
    run_json = os.path.join(d, "RUN.json")
    write_json(run_json, run_rec)
    log(f"run {cfg.run_id}: record {run_json}")
    chosen = {s.sid for s in todo}
    ok_all = True
    try:
        for st in todo:
            unmet = [x for x in st.deps if x not in chosen and every[x].state(cfg) != "validated"]
            failed = [x for x in st.deps if x in chosen and every[x].state(cfg) != "validated"]
            if unmet or failed:
                log(f"{st.sid}: skipped, dependencies not validated: {unmet + failed}")
                ok_all = False
                continue
            ok_all &= run_stage(cfg, st, vers, args.rerun_stale)
    except BaseException as e:
        run_rec.update(end=now(), status=f"aborted: {type(e).__name__} {e}")
        write_json(run_json, run_rec)
        raise
    rc = 0 if ok_all else 1
    run_rec.update(end=now(), status="all selected stages validated" if ok_all else "some stages not validated",
                   exit_code=rc, stage_states={s.sid: s.state(cfg) for s in todo})
    write_json(run_json, run_rec)
    log(f"run finished: {run_rec['status']}")
    return rc


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
