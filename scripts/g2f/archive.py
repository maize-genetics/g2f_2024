"""Read-only access to ONT FASTQ delivered as uncompressed tar archives.

Sequencing facilities often deliver MinKNOW output as one tar per run, e.g.

    <run>_fastq.tar
      ./fastq_pass/<flowcell>_pass_<protocol>_<acq>_<n>.fastq.gz   (4000 reads each)

The workflow never extracts these archives. Members are streamed in a fixed
(natural) order straight from the tar into the consumer (seqkit, sketching),
so no copy of the reads is written to disk.

Member policy (deterministic, recorded in the inventory):
  * FASTQ members (.fastq/.fq, optionally .gz) under a directory whose name
    contains "fail" (MinKNOW fastq_fail/) are EXCLUDED: they did not pass the
    basecaller's quality filter.
  * all other FASTQ members are USED;
  * directories are ignored;
  * anything else (summaries, FAST5/POD5, BAM, ...) is an error, because the
    archive's content is then not what the workflow assumes;
  * used members must all be gzip-compressed or all uncompressed, so that
    their concatenation is a valid (multi-member gzip) FASTQ stream.
"""
import gzip
import io
import re
import tarfile

FASTQ_MEMBER_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
COPY_CHUNK = 1 << 22


class ArchiveError(ValueError):
    pass


def is_tar_name(name):
    return name.lower().endswith(".tar")


def _natural_key(name):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def _is_fail(name):
    return any("fail" in part.lower() for part in name.split("/")[:-1])


def list_members(path):
    """Return {"use": [...], "excluded": [...], "errors": [...], "compressed": bool|None}.

    Each member entry is {"name", "size", "offset"}; only headers are read
    (the archive is seeked, not decompressed).
    """
    use, excluded, errors = [], [], []
    try:
        with tarfile.open(path, "r:") as tf:
            for m in tf:
                if m.isdir():
                    continue
                rec = {"name": m.name, "size": m.size, "offset": m.offset_data}
                if not m.isfile():
                    errors.append(f"{m.name}: not a regular file (link or device)")
                elif not FASTQ_MEMBER_RE.search(m.name):
                    errors.append(f"{m.name}: not a FASTQ file")
                elif _is_fail(m.name):
                    excluded.append(dict(rec, reason="fastq_fail directory"))
                else:
                    use.append(rec)
    except (tarfile.TarError, OSError, EOFError) as exc:
        raise ArchiveError(f"cannot read tar archive: {exc}")
    use.sort(key=lambda r: _natural_key(r["name"]))
    kinds = {r["name"].lower().endswith(".gz") for r in use}
    if len(kinds) > 1:
        errors.append("archive mixes gzip-compressed and uncompressed FASTQ members")
    return {"use": use, "excluded": excluded, "errors": errors,
            "compressed": kinds.pop() if len(kinds) == 1 else None}


def stream_raw(path, out, members=None):
    """Write the raw bytes of the used members, in order, to binary stream `out`.

    gzip members are passed through compressed (their concatenation is a valid
    multi-member gzip stream). Returns (n_members, n_bytes)."""
    listing = members or list_members(path)
    if listing["errors"]:
        raise ArchiveError("; ".join(listing["errors"]))
    n_bytes = 0
    with open(path, "rb") as fh:
        for m in listing["use"]:
            fh.seek(m["offset"])
            left = m["size"]
            while left:
                buf = fh.read(min(COPY_CHUNK, left))
                if not buf:
                    raise ArchiveError(f"{m['name']}: archive truncated "
                                       f"({left} of {m['size']} bytes missing)")
                out.write(buf)
                left -= len(buf)
            n_bytes += m["size"]
    return len(listing["use"]), n_bytes


def open_text(path, members=None):
    """Iterate decompressed text lines of all used members (for Python consumers)."""
    listing = members or list_members(path)
    if listing["errors"]:
        raise ArchiveError("; ".join(listing["errors"]))
    with open(path, "rb") as fh:
        for m in listing["use"]:
            fh.seek(m["offset"])
            raw = _Bounded(fh, m["size"])
            stream = gzip.GzipFile(fileobj=raw) if m["name"].lower().endswith(".gz") else raw
            yield from io.TextIOWrapper(stream, encoding="ascii", errors="replace")


class _Bounded(io.RawIOBase):
    """Read-only view of `size` bytes of `fh` from its current position."""

    def __init__(self, fh, size):
        self._fh, self._left = fh, size

    def readable(self):
        return True

    def readinto(self, b):
        if self._left <= 0:
            return 0
        n = self._fh.readinto(memoryview(b)[: min(len(b), self._left)])
        if not n:
            raise ArchiveError("archive truncated inside a member")
        self._left -= n
        return n


def head_lines(path, n_lines, members=None):
    out = []
    for line in open_text(path, members):
        out.append(line)
        if len(out) >= n_lines:
            break
    return out


def summary(listing):
    return {
        "n_members": len(listing["use"]),
        "member_bytes": sum(m["size"] for m in listing["use"]),
        "n_excluded_members": len(listing["excluded"]),
        "excluded_members": [m["name"] for m in listing["excluded"]][:20],
        "first_member": listing["use"][0]["name"] if listing["use"] else None,
        "last_member": listing["use"][-1]["name"] if listing["use"] else None,
        "members_compressed": listing["compressed"],
    }
