"""Deterministic discovery and lightweight identity inspection of read inputs.

Discovery is cheap: it stats files and reads only the first few records (and
the first MiB for a fingerprint). It never hashes whole files. Full-file
checks (checksums, read counts, complete movie / read-ID census) run inside
the workflow and gate assembly.

Layout:  <data_dir>/<SAMPLE>/{hifi,ont}/...  plus optional read-only
external references from config/local/local.yaml.
"""
import gzip
import hashlib
import os
import re
import struct
from dataclasses import asdict, dataclass, field

from . import archive

TECHS = ("hifi", "ont")
READ_EXTS = {
    ".bam": "bam",
    ".fastq": "fastq", ".fq": "fastq",
    ".fasta": "fasta", ".fa": "fasta", ".fna": "fasta",
}
RAW_SIGNAL_EXTS = (".pod5", ".fast5", ".blow5", ".slow5")
UNSUPPORTED_EXTS = (".sam", ".cram", ".h5", ".bax.h5", ".bas.h5")

# PacBio read names: <movie>/<zmw>/ccs | /ccs/fwd | /ccs/rev | /<start>_<end>
PACBIO_NAME_RE = re.compile(
    r"^(m\d+[A-Za-z]?_\d{6}_\d{6}(?:_s\d+)?)/(\d+)/(ccs(?:/fwd|/rev)?|\d+_\d+)(?:\s|$)")
UUID_RE = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
FINGERPRINT_BYTES = 1 << 20
MAX_SAMPLE_BYTES = 64 << 20


@dataclass
class Issue:
    severity: str          # "error" | "warning"
    code: str
    sample: str
    message: str
    action: str = ""
    paths: list = field(default_factory=list)

    def render(self):
        s = f"[{self.severity.upper()}] {self.code} ({self.sample}): {self.message}"
        for p in self.paths:
            s += f"\n    - {p}"
        if self.action:
            s += f"\n    Action: {self.action}"
        return s


@dataclass
class InputFile:
    sample: str
    tech: str
    path: str
    key: str               # relative path within the tech dir, or basename (external)
    source: str            # data_dir | external
    status: str = "use"    # use | excluded | ignored
    format: str = ""
    compressed: bool = False
    size: int = 0
    fingerprint: str = ""
    input_id: str = ""
    identity: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)


# --------------------------------------------------------------------------
# Low-level inspection
# --------------------------------------------------------------------------

def _is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def fingerprint(path, size):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read(FINGERPRINT_BYTES))
    h.update(str(size).encode())
    return "sha256-head1MiB:" + h.hexdigest()[:32]


def _strip_exts(name):
    if archive.is_tar_name(name):
        return name[: -len(".tar")], "tar"
    base = name
    if base.endswith(".gz"):
        base = base[:-3]
    for ext in READ_EXTS:
        if base.lower().endswith(ext):
            return base[: -len(ext)], READ_EXTS[ext]
    return base, None


def parse_bam_head(path, n_records):
    """Return (header_text, [(read_name, read_len), ...]) from a BAM file."""
    with gzip.open(path, "rb") as fh:
        if fh.read(4) != b"BAM\x01":
            raise ValueError("not a BAM file (bad magic)")
        (l_text,) = struct.unpack("<i", fh.read(4))
        text = fh.read(l_text).decode("ascii", "replace").rstrip("\x00")
        (n_ref,) = struct.unpack("<i", fh.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", fh.read(4))
            fh.read(l_name + 4)
        recs = []
        while len(recs) < n_records:
            raw = fh.read(4)
            if len(raw) < 4:
                break
            (block_size,) = struct.unpack("<i", raw)
            block = fh.read(block_size)
            l_read_name = block[8]
            flag = struct.unpack("<H", block[14:16])[0]
            (l_seq,) = struct.unpack("<i", block[16:20])
            name = block[32:32 + l_read_name].rstrip(b"\x00").decode("ascii", "replace")
            if not flag & 0x900:
                recs.append((name, l_seq))
    return text, recs


def parse_sam_header(text):
    out = {"RG": [], "PG": [], "SQ": 0}
    for line in text.splitlines():
        if line.startswith("@SQ"):
            out["SQ"] += 1
        elif line.startswith(("@RG", "@PG")):
            tags = {}
            for tok in line.split("\t")[1:]:
                if ":" in tok:
                    k, v = tok.split(":", 1)
                    tags[k] = v
            out[line[1:3]].append(tags)
    return out


def parse_fastx_head(path, compressed, n_records):
    """Return (format, [(header_line, read_len, mean_q_or_None), ...])."""
    opener = gzip.open if compressed else open
    with opener(path, "rt", errors="replace") as fh:
        return parse_fastx_stream(fh, n_records)


def parse_fastx_stream(fh, n_records):
    """parse_fastx_head on an open text stream (e.g. a member of a tar archive)."""
    recs, read_bytes = [], 0
    first = fh.readline()
    read_bytes += len(first)
    if first.startswith("@"):
        fmt, header = "fastq", first
        while header and len(recs) < n_records:
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()
            read_bytes += len(seq) + len(qual) + len(plus)
            if not plus.startswith("+") or len(qual) != len(seq):
                raise ValueError("malformed FASTQ record near the start of the file")
            q = (sum(ord(c) - 33 for c in qual) / len(qual)) if qual else None
            recs.append((header[1:].rstrip("\n"), len(seq), q))
            header = fh.readline()
            if read_bytes > MAX_SAMPLE_BYTES:
                break
    elif first.startswith(">"):
        fmt, header, length = "fasta", first, 0
        for line in fh:
            read_bytes += len(line)
            if line.startswith(">"):
                recs.append((header[1:].rstrip("\n"), length, None))
                if len(recs) >= n_records:
                    header = None
                    break
                header, length = line, 0
            else:
                length += len(line.strip())
            if read_bytes > MAX_SAMPLE_BYTES:
                raise ValueError("first records exceed 64 MiB; this does not look "
                                 "like a read file (assembly?)")
        if header is not None and length:
            recs.append((header[1:].rstrip("\n"), length, None))
    else:
        raise ValueError("file does not start with '@' (FASTQ) or '>' (FASTA)")
    return fmt, recs


def _ont_header_fields(header):
    toks = re.split(r"[ \t]+", header.strip())
    fields = {"read_id": toks[0]}
    for tok in toks[1:]:
        if "=" in tok:
            k, v = tok.split("=", 1)
            fields[k] = v
        elif re.match(r"^[A-Za-z][A-Za-z0-9]:[AifZHB]:", tok):
            k, _, v = tok.split(":", 2)
            if k == "RG":
                fields["RG"] = v
                fields.setdefault("runid", v.split("_")[0])
    return fields


# --------------------------------------------------------------------------
# Per-file inspection
# --------------------------------------------------------------------------

def inspect_file(inp, n_records, issues):
    path, name = inp.path, os.path.basename(inp.path)
    try:
        inp.size = os.path.getsize(path)
    except OSError as exc:
        issues.append(Issue("error", "UNREADABLE", inp.sample, f"cannot stat file: {exc}",
                            "Check the path and permissions.", [path]))
        return
    if inp.size == 0:
        issues.append(Issue("error", "EMPTY_FILE", inp.sample, "file is empty",
                            "Remove the file or replace it with the complete data.", [path]))
        return
    stem, ext_fmt = _strip_exts(name)
    inp.input_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)
    inp.compressed = _is_gzip(path)
    inp.fingerprint = fingerprint(path, inp.size)
    if ext_fmt == "bam":
        _inspect_bam(inp, n_records, issues)
    elif ext_fmt == "tar":
        _inspect_tar(inp, n_records, issues)
    else:
        try:
            fmt, recs = parse_fastx_head(path, inp.compressed, n_records)
        except (ValueError, OSError, EOFError) as exc:
            issues.append(Issue("error", "BAD_FORMAT", inp.sample, f"cannot parse reads: {exc}",
                                "Check the file is complete and in FASTQ/FASTA format.", [path]))
            return
        if fmt != ext_fmt:
            issues.append(Issue("error", "FORMAT_MISMATCH", inp.sample,
                                f"extension says {ext_fmt} but content is {fmt}",
                                "Rename the file with the correct extension.", [path]))
            return
        inp.format = fmt
        if inp.tech == "ont" and fmt == "fasta":
            issues.append(Issue("error", "UNSUPPORTED_ONT_FASTA", inp.sample,
                                "ONT reads must be FASTQ or BAM (base qualities are needed for QC)",
                                "Provide the basecalled FASTQ/BAM.", [inp.path]))
            return
        if not recs:
            issues.append(Issue("error", "NO_READS", inp.sample, "no reads found",
                                "Remove or replace the file.", [path]))
            return
        if inp.tech == "hifi":
            _identify_hifi_fastx(inp, recs, issues)
        else:
            _identify_ont(inp, [(_ont_header_fields(h), n, q) for h, n, q in recs], issues)


def _inspect_tar(inp, n_records, issues):
    """ONT FASTQ delivered as a tar archive (see g2f.archive): read in place."""
    if inp.tech != "ont":
        issues.append(Issue("error", "UNSUPPORTED_FORMAT", inp.sample,
                            "tar archives are supported only for ONT FASTQ deliveries",
                            "Provide HiFi reads as BAM or FASTQ(.gz).", [inp.path]))
        return
    try:
        listing = archive.list_members(inp.path)
    except archive.ArchiveError as exc:
        issues.append(Issue("error", "BAD_FORMAT", inp.sample, str(exc),
                            "Check the archive is complete (compare its size with the source).",
                            [inp.path]))
        return
    if listing["errors"]:
        issues.append(Issue("error", "BAD_ARCHIVE_CONTENT", inp.sample,
                            "archive contains members other than FASTQ reads: "
                            + "; ".join(listing["errors"][:5]),
                            "Provide an archive of basecalled FASTQ only, or list it under "
                            "exclude_inputs in samples.yaml.", [inp.path]))
        return
    if not listing["use"]:
        issues.append(Issue("error", "NO_READS", inp.sample,
                            "archive contains no usable (pass) FASTQ members",
                            "Remove or replace the file.", [inp.path]))
        return
    inp.format = "tar"
    inp.compressed = bool(listing["compressed"])
    inp.identity.update(archive.summary(listing))
    try:
        stream = archive.open_text(inp.path, listing)
        fmt, recs = parse_fastx_stream(_LineReader(stream), n_records)
    except (ValueError, OSError, EOFError) as exc:
        issues.append(Issue("error", "BAD_FORMAT", inp.sample,
                            f"cannot parse reads in archive: {exc}",
                            "Check the archive members are FASTQ.", [inp.path]))
        return
    if fmt != "fastq" or not recs:
        issues.append(Issue("error", "BAD_FORMAT", inp.sample,
                            "archive members must be FASTQ (base qualities are needed for QC)",
                            "Provide basecalled FASTQ.", [inp.path]))
        return
    _identify_ont(inp, [(_ont_header_fields(h), n, q) for h, n, q in recs], issues)


class _LineReader:
    """Minimal readline()/iteration adapter over a line iterator."""

    def __init__(self, lines):
        self._it = iter(lines)

    def readline(self):
        return next(self._it, "")

    def __iter__(self):
        return self._it


def _inspect_bam(inp, n_records, issues):
    inp.format = "bam"
    try:
        text, recs = parse_bam_head(inp.path, n_records)
    except (ValueError, OSError, EOFError, struct.error, IndexError) as exc:
        issues.append(Issue("error", "BAD_FORMAT", inp.sample, f"cannot parse BAM: {exc}",
                            "Check the BAM is complete (not truncated).", [inp.path]))
        return
    hdr = parse_sam_header(text)
    rgs = hdr["RG"]
    platforms = sorted({rg.get("PL", "").upper() for rg in rgs if rg.get("PL")})
    inp.identity["platforms"] = platforms
    inp.identity["n_sq"] = hdr["SQ"]
    inp.identity["programs"] = sorted({pg.get("PN", pg.get("ID", "")) for pg in hdr["PG"]})
    if not recs:
        issues.append(Issue("error", "NO_READS", inp.sample, "BAM contains no primary reads",
                            "Remove or replace the file.", [inp.path]))
        return
    if inp.tech == "hifi":
        readtypes = sorted({m.group(1) for rg in rgs
                            for m in [re.search(r"READTYPE=([A-Z_]+)", rg.get("DS", ""))] if m})
        inp.identity["readtypes"] = readtypes
        inp.identity["read_groups"] = sorted(
            f"{rg.get('PU', '?')}|{rg.get('BC', rg.get('ID', '?'))}" for rg in rgs)
        name = os.path.basename(inp.path).lower()
        if "PACBIO" not in platforms:
            issues.append(Issue("error", "NOT_PACBIO", inp.sample,
                                f"BAM in hifi/ does not have PL:PACBIO read groups (PL={platforms})",
                                "Move non-PacBio data out of hifi/.", [inp.path]))
            return
        if readtypes != ["CCS"]:
            issues.append(Issue("error", "NOT_HIFI", inp.sample,
                                f"PacBio BAM READTYPE={readtypes or 'missing'}; only CCS/HiFi "
                                "BAMs are accepted (subreads must be processed with ccs first)",
                                "Provide hifi_reads BAMs, or exclude this file in samples.yaml.",
                                [inp.path]))
            return
        if "fail_reads" in name or ".unassigned" in name or "unbarcoded" in name:
            issues.append(Issue("error", "NOT_HIFI", inp.sample,
                                "file name indicates failed-quality or unassigned-barcode reads",
                                "Move it out of hifi/ or list it under exclude_inputs.",
                                [inp.path]))
            return
        if hdr["SQ"]:
            issues.append(Issue("error", "ALIGNED_BAM", inp.sample,
                                "HiFi BAM has @SQ lines (aligned); unaligned hifi_reads BAMs expected",
                                "Provide the unaligned hifi_reads BAM.", [inp.path]))
            return
        _identify_hifi_names(inp, [(n, l, None) for n, l in recs], issues)
    else:
        if "PACBIO" in platforms:
            issues.append(Issue("error", "WRONG_TECH", inp.sample,
                                "PacBio BAM found in ont/", "Move it to hifi/.", [inp.path]))
            return
        if hdr["SQ"]:
            issues.append(Issue("warning", "ALIGNED_BAM", inp.sample,
                                "ONT BAM is aligned; only primary records will be used",
                                "", [inp.path]))
        runids = sorted({m.group(1) for rg in rgs
                         for m in [re.search(r"runid=(\S+)", rg.get("DS", ""))] if m})
        models = sorted({m.group(1) for rg in rgs
                         for m in [re.search(r"basecall_model=(\S+)", rg.get("DS", ""))] if m})
        recs2 = [({"read_id": n, **({"runid": runids[0]} if len(runids) == 1 else {})}, l, None)
                 for n, l in recs]
        _identify_ont(inp, recs2, issues)
        if models:
            inp.identity["basecall_models"] = models


def _identify_hifi_fastx(inp, recs, issues):
    _identify_hifi_names(inp, recs, issues)


def _identify_hifi_names(inp, recs, issues):
    movies, kinds, unparsed = set(), set(), 0
    for header, _, _ in recs:
        m = PACBIO_NAME_RE.match(header)
        if not m:
            unparsed += 1
            continue
        movies.add(m.group(1))
        kinds.add("subread" if "_" in m.group(3) else m.group(3))
    inp.identity.update({
        "sampled_movies": sorted(movies),
        "sampled_read_ids": [h.split()[0] for h, _, _ in recs[:5]],
        "sampled_n": len(recs),
        "sampled_len_min": min(n for _, n, _ in recs),
        "sampled_len_max": max(n for _, n, _ in recs),
    })
    if "subread" in kinds:
        issues.append(Issue("error", "NOT_HIFI", inp.sample,
                            "read names look like PacBio subreads (movie/zmw/start_end)",
                            "Provide CCS/HiFi reads.", [inp.path]))
    elif kinds & {"ccs/fwd", "ccs/rev"}:
        issues.append(Issue("error", "BYSTRAND_CCS", inp.sample,
                            "by-strand CCS reads (ccs/fwd, ccs/rev) would double-count ZMWs",
                            "Provide standard (not by-strand) HiFi reads.", [inp.path]))
    elif unparsed:
        if inp.notes and "ack:unverified_hifi" in inp.notes:
            return
        issues.append(Issue("error", "HIFI_IDENTITY_UNVERIFIED", inp.sample,
                            f"{unparsed}/{len(recs)} sampled read names are not PacBio "
                            "movie/zmw/ccs names; cannot verify these are HiFi reads or "
                            "detect duplicate copies",
                            "If you have verified the file, add under acknowledge: "
                            f"{{{inp.key}: {{unverified_hifi: <reason>}}}} in samples.yaml.",
                            [inp.path]))


def _identify_ont(inp, recs, issues):
    ids = [f["read_id"] for f, _, _ in recs]
    lens = [n for _, n, _ in recs]
    quals = [q for _, _, q in recs if q is not None]
    if any(PACBIO_NAME_RE.match(i) for i in ids):
        issues.append(Issue("error", "WRONG_TECH", inp.sample,
                            "PacBio read names found in ont/", "Move the file to hifi/.",
                            [inp.path]))
        return
    inp.identity.update({
        "sampled_read_ids": ids,
        "sampled_runids": sorted({f["runid"] for f, _, _ in recs if "runid" in f}),
        "sampled_flow_cells": sorted({f["flow_cell_id"] for f, _, _ in recs if "flow_cell_id" in f}),
        "sampled_sample_ids": sorted({f["sample_id"] for f, _, _ in recs if f.get("sample_id")}),
        "sampled_n": len(recs),
        "sampled_len_min": min(lens),
        "sampled_len_max": max(lens),
        "sampled_mean_q": round(sum(quals) / len(quals), 2) if quals else None,
    })
    models = sorted({f.get("basecall_model_version_id") or f.get("model_version_id")
                     for f, _, _ in recs
                     if f.get("basecall_model_version_id") or f.get("model_version_id")})
    if models:
        inp.identity["basecall_models"] = models
    if not all(UUID_RE.match(i) for i in ids):
        issues.append(Issue("warning", "ONT_IDS_NOT_UUID", inp.sample,
                            "ONT read IDs are not UUIDs; duplicate detection relies on IDs as given",
                            "", [inp.path]))


# --------------------------------------------------------------------------
# Sample- and project-level discovery
# --------------------------------------------------------------------------

def _walk(tech_dir):
    for root, dirs, files in os.walk(tech_dir, followlinks=True):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        for f in sorted(files):
            if not f.startswith("."):
                full = os.path.join(root, f)
                yield os.path.relpath(full, tech_dir), full


def _classify(key, path, sample, tech, aux_suffixes, issues):
    """Return status for a file name: 'read', 'ignored', or None (error recorded)."""
    low = os.path.basename(path).lower()
    if low.endswith(RAW_SIGNAL_EXTS):
        issues.append(Issue("error", "RAW_SIGNAL", sample,
                            "raw ONT signal (POD5/FAST5/SLOW5) found; these must be basecalled "
                            "first by a separately configured basecalling step",
                            "Basecall to FASTQ/BAM and place the basecalled reads in ont/; move "
                            "raw signal files out of the data directory.", [path]))
        return None
    if low.endswith(UNSUPPORTED_EXTS):
        issues.append(Issue("error", "UNSUPPORTED_FORMAT", sample,
                            "unsupported sequencing file format",
                            "Convert to BAM/FASTQ/FASTA or move it out of the data directory.",
                            [path]))
        return None
    if _strip_exts(os.path.basename(path))[1]:
        return "read"
    if any(low.endswith(s.lower()) for s in aux_suffixes):
        return "ignored"
    issues.append(Issue("error", "UNKNOWN_FILE", sample,
                        f"unrecognised file in {tech}/ (not a supported read file or a known "
                        "auxiliary file)",
                        "Move it out of the data directory, or list it under exclude_inputs "
                        "in samples.yaml with a reason.", [path]))
    return None


def discover(cfg):
    """Discover all inputs. Returns (inputs_by_sample, issues)."""
    sci = cfg["science"].get("discovery") or {}
    n_records = int(sci.get("sample_records", 200))
    aux = sci.get("auxiliary_suffixes") or []
    prefilter_len = sci.get("ont_prefilter_suspect_min_len")
    issues, inputs = [], {}

    candidates = []   # (sample, tech, key, path, source)
    data_dir = cfg["data_dir"]
    if os.path.isdir(data_dir):
        for sample in sorted(os.listdir(data_dir)):
            sdir = os.path.join(data_dir, sample)
            if sample.startswith(".") or not os.path.isdir(sdir):
                continue
            for extra in sorted(set(os.listdir(sdir)) - set(TECHS)):
                if not extra.startswith("."):
                    issues.append(Issue("error", "UNEXPECTED_ENTRY", sample,
                                        f"unexpected entry '{extra}' in sample directory "
                                        "(only hifi/ and ont/ are allowed)",
                                        "Move it into hifi/ or ont/, or out of data/.",
                                        [os.path.join(sdir, extra)]))
            for tech in TECHS:
                tdir = os.path.join(sdir, tech)
                if os.path.isdir(tdir):
                    for key, full in _walk(tdir):
                        candidates.append((sample, tech, key, full, "data_dir"))
    for sample, techs in sorted((cfg["local"].get("external_inputs") or {}).items()):
        for tech, paths in sorted((techs or {}).items()):
            if tech not in TECHS:
                issues.append(Issue("error", "BAD_EXTERNAL", sample,
                                    f"external_inputs technology must be hifi or ont, got {tech}"))
                continue
            for p in paths or []:
                p = os.path.expanduser(p)
                if not os.path.isabs(p):
                    issues.append(Issue("error", "BAD_EXTERNAL", sample,
                                        "external input paths must be absolute", "", [p]))
                    continue
                if os.path.isdir(p):
                    for key, full in _walk(p):
                        candidates.append((sample, tech, key, full, "external"))
                else:
                    candidates.append((sample, tech, os.path.basename(p), p, "external"))

    only = set(cfg.get("only_samples") or [])
    # Sample directories with no files yet (e.g. data still in transfer) are
    # skipped, but reported so a run never silently omits a sample.
    # A sample requested explicitly (only_samples) that has no reads is an
    # error, unless it has a registered historical assembly (imported_qc).
    with_files = {c[0] for c in candidates}
    if os.path.isdir(data_dir):
        for sample in sorted(os.listdir(data_dir)):
            if (sample.startswith(".") or sample in with_files or sample in only
                    or not os.path.isdir(os.path.join(data_dir, sample))
                    or (only and sample not in only)):
                continue
            issues.append(Issue("warning", "NO_INPUTS", sample,
                                "sample directory contains no files; sample skipped",
                                "Add reads to hifi/ (and ont/ if available).",
                                [os.path.join(data_dir, sample)]))
    imported = cfg.get("external_assemblies") or {}
    for sample in sorted(only - with_files):
        where = os.path.join(data_dir, sample)
        if sample in imported:
            issues.append(Issue("warning", "NO_INPUTS", sample,
                                "requested sample has no reads; only its registered "
                                "historical assemblies (imported_qc) can be processed",
                                "Add reads to hifi/ (and ont/ if available).", [where]))
        else:
            issues.append(Issue("error", "REQUESTED_SAMPLE_NO_INPUTS", sample,
                                "sample requested in only_samples has no read files "
                                + ("(directory is empty)" if os.path.isdir(where)
                                   else "(no such sample directory)"),
                                "Add reads to hifi/ (and ont/ if available), or remove it "
                                "from only_samples.", [where]))
    seen_paths = {}
    for sample, tech, key, path, source in candidates:
        if only and sample not in only:
            continue
        if not SAMPLE_ID_RE.match(sample):
            issues.append(Issue("error", "BAD_SAMPLE_ID", sample,
                                "sample directory names may contain only letters, digits, "
                                "'.', '_' and '-'", "Rename the directory."))
            continue
        if "," in path or any(c.isspace() for c in path):
            issues.append(Issue("error", "BAD_PATH", sample,
                                "read file paths may not contain commas or whitespace (hifiasm "
                                "takes comma-separated --ul lists)", "Rename the file or directory.",
                                [path]))
            continue
        real = os.path.realpath(path)
        if real in seen_paths:
            issues.append(Issue("error", "SAME_FILE_TWICE", sample,
                                "the same file is referenced more than once",
                                "Remove the duplicate reference.", [seen_paths[real], path]))
            continue
        seen_paths[real] = path
        entry = cfg["samples"].get(sample) or {}
        excluded = entry.get("exclude_inputs") or {}
        acks = entry.get("acknowledge") or {}
        inp = InputFile(sample=sample, tech=tech, path=path, key=key, source=source)
        bn = os.path.basename(path)
        if key in excluded or bn in excluded:
            inp.status = "excluded"
            inp.notes.append(f"excluded: {excluded.get(key) or excluded.get(bn)}")
            inputs.setdefault(sample, []).append(inp)
            continue
        kind = _classify(key, path, sample, tech, aux, issues)
        if kind is None:
            continue
        if kind == "ignored":
            inp.status = "ignored"
            inputs.setdefault(sample, []).append(inp)
            continue
        ack = acks.get(key) or acks.get(bn) or {}
        inp.notes.extend(f"ack:{k}" for k in sorted(ack))
        inspect_file(inp, n_records, issues)
        inputs.setdefault(sample, []).append(inp)
        if (tech == "ont" and prefilter_len and inp.identity.get("sampled_len_min", 0)
                >= prefilter_len and "ack:prefiltered_ont" not in inp.notes):
            issues.append(Issue("error", "ONT_POSSIBLY_PREFILTERED", sample,
                                f"all {inp.identity.get('sampled_n')} sampled reads are >= "
                                f"{prefilter_len} bp; this file looks already length-filtered, "
                                "so its yield cannot be reported as raw ONT yield",
                                "Prefer the raw basecalled reads. If only this file exists, add "
                                f"under acknowledge: {{{key}: {{prefiltered_ont: <reason>}}}} in "
                                "samples.yaml.", [path]))

    for sample in sorted(inputs):
        _check_sample(sample, inputs[sample], issues)
    _check_cross_sample(inputs, issues)
    return inputs, issues


def _check_sample(sample, files, issues):
    used = [f for f in files if f.status == "use" and f.format]
    ids = {}
    for f in used:
        uid = f"{f.tech}:{f.input_id}"
        if uid in ids:
            issues.append(Issue("error", "INPUT_ID_COLLISION", sample,
                                "two inputs map to the same input ID (same name, different "
                                "extension or directory)",
                                "Rename one file, or exclude the redundant copy.",
                                [ids[uid].path, f.path]))
        ids[uid] = f
    by_fp = {}
    for f in used:
        if f.fingerprint in by_fp:
            issues.append(Issue("error", "DUPLICATE_COPY", sample,
                                "files have identical size and leading content (duplicate copies)",
                                "Keep one copy; list the other under exclude_inputs.",
                                [by_fp[f.fingerprint].path, f.path]))
        by_fp.setdefault(f.fingerprint, f)
    hifi = [f for f in used if f.tech == "hifi"]
    for i, a in enumerate(hifi):
        for b in hifi[i + 1:]:
            shared = set(a.identity.get("sampled_movies", [])) & set(b.identity.get("sampled_movies", []))
            if a.format == b.format == "bam":
                rg_a, rg_b = set(a.identity.get("read_groups", [])), set(b.identity.get("read_groups", []))
                if rg_a & rg_b:
                    issues.append(Issue("error", "HIFI_DUPLICATE_RUN", sample,
                                        f"BAMs share read groups {sorted(rg_a & rg_b)}",
                                        "Keep one; list the other under exclude_inputs.",
                                        [a.path, b.path]))
                continue
            if shared:
                derived = [x for x in (a, b) if x.format != "bam"]
                hint = ("BAM is the authoritative source; a FASTQ/FASTA sharing its movies is "
                        "presumed derived from it." if len(derived) == 1 else "")
                exclude = derived[0] if len(derived) == 1 else b
                issues.append(Issue("error", "HIFI_POSSIBLE_DERIVED_COPY", sample,
                                    f"inputs share PacBio movie(s) {sorted(shared)}; using both "
                                    f"would double-count reads. {hint}",
                                    "Decide which file is authoritative and add, e.g., "
                                    f"exclude_inputs: {{{exclude.key}: <reason>}} to samples.yaml.",
                                    [a.path, b.path]))
    ont = [f for f in used if f.tech == "ont"]
    for i, a in enumerate(ont):
        for b in ont[i + 1:]:
            shared = set(a.identity.get("sampled_read_ids", [])) & set(b.identity.get("sampled_read_ids", []))
            if shared:
                issues.append(Issue("error", "ONT_DUPLICATE_READS", sample,
                                    f"ONT inputs share {len(shared)} sampled read IDs (duplicate "
                                    "or derived copies, e.g. raw and filtered versions)",
                                    "Keep the raw basecalled file; exclude the other in samples.yaml.",
                                    [a.path, b.path]))
    if any(f.tech == "hifi" for f in files) and not hifi:
        issues.append(Issue("error", "NO_USABLE_HIFI", sample,
                            "HiFi files are present but none are usable", "See errors above."))
    elif not any(f.tech == "hifi" for f in files):
        issues.append(Issue("error", "NO_HIFI", sample, "no HiFi inputs found",
                            f"Place HiFi reads in data/{sample}/hifi/ or remove the sample directory."))


def _check_cross_sample(inputs, issues):
    seen = {}
    for sample in sorted(inputs):
        for f in inputs[sample]:
            if f.status != "use" or not f.fingerprint:
                continue
            other = seen.get(f.fingerprint)
            if other and other.sample != sample:
                issues.append(Issue("error", "CROSS_SAMPLE_DUPLICATE", sample,
                                    f"identical file also assigned to sample {other.sample}",
                                    "Each read file must belong to exactly one sample.",
                                    [other.path, f.path]))
            seen.setdefault(f.fingerprint, f)


def select_method(sample, files, cfg):
    """Return (method_id, reason). Assumes discovery errors were already handled."""
    entry = cfg["samples"].get(sample) or {}
    policy = cfg["science"]["method_policy"]
    if entry.get("method_override"):
        return entry["method_override"], f"override: {entry['override_reason']}"
    ont_used = [f for f in files if f.tech == "ont" and f.status == "use"]
    ont_excluded = [f for f in files if f.tech == "ont" and f.status == "excluded"]
    if ont_used:
        return policy["ont_present"], (f"{len(ont_used)} usable ONT input(s) present "
                                       "(fully validated before assembly)")
    if ont_excluded:
        return policy["ont_absent"], (f"ONT inputs present but all {len(ont_excluded)} excluded "
                                      "by human decision in samples.yaml")
    return policy["ont_absent"], "no ONT inputs present"


def inventory_records(sample, files, cfg):
    method, reason = select_method(sample, files, cfg)
    return {
        "sample": sample,
        "aliases": (cfg["samples"].get(sample) or {}).get("aliases", []),
        "selected_method": method,
        "selection_reason": reason,
        "inputs": [asdict(f) for f in sorted(files, key=lambda x: (x.tech, x.key, x.path))],
    }


def stable_view(record):
    """Fields that define input identity; used to detect changed inputs."""
    return {
        "selected_method": record["selected_method"],
        "inputs": [[i["tech"], i["path"], i["status"], i["size"], i["fingerprint"], i["input_id"],
                    sorted(n for n in i.get("notes", []) if n.startswith(("ack:", "excluded:")))]
                   for i in record["inputs"]],
    }


def is_prefiltered(inp):
    notes = inp.notes if hasattr(inp, "notes") else inp.get("notes", [])
    return "ack:prefiltered_ont" in notes
