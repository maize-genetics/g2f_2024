"""Streaming FracMinHash sketches of sequencing reads for genotype fingerprinting.

One pass over (a bounded prefix of) one read file, without writing reads to
disk. Reads are taken from:
  * PacBio BAM (unaligned or aligned; primary records via `samtools fasta`),
  * FASTQ / FASTA, optionally gzip-compressed,
  * ONT tar archives of FASTQ members (g2f.archive; pass members only).

Sampling (all deterministic and recorded):
  * the first `max_bases` bases of the stream are used (0 = all). PacBio BAMs
    are in ZMW order and ONT output is in acquisition order; neither is
    related to genome position, so a prefix is an unbiased genomic sample;
  * optionally, a seeded read-level subsample (`read_fraction` < 1) decided by
    a hash of the read name, for inputs whose order might be genomic;
  * every used read is assigned to replicate half A or B by a seeded hash of
    its read name. Each half is sketched separately; the full sketch is their
    sum. Halves provide within-run technical replicates for free.

Sketches: sourmash FracMinHash (MurmurHash3 of canonical k-mers, seed 42),
keeping hashes below 2^64/scaled, with abundance tracking. One MinHash per
(k, half) is written to a sourmash-format signature file (.sig.gz) so that the
sketches can also be inspected with the sourmash CLI.

Only `sketch()` needs the sourmash library; the readers are standard library.
"""
import gzip
import hashlib
import json
import os
import subprocess
import time
import zlib

from . import archive

HALVES = ("A", "B")


# --------------------------------------------------------------------------
# Readers: yield (read_name, sequence)
# --------------------------------------------------------------------------

def detect_format(path):
    low = os.path.basename(path).lower()
    if archive.is_tar_name(low):
        return "tar"
    if low.endswith(".bam"):
        return "bam"
    with open(path, "rb") as fh:
        head = fh.read(2)
    opener = gzip.open if head == b"\x1f\x8b" else open
    with opener(path, "rt", errors="replace") as fh:
        c = fh.read(1)
    if c == "@":
        return "fastq"
    if c == ">":
        return "fasta"
    raise ValueError(f"cannot determine read format of {path}")


def _fastx_lines(lines):
    """Parse FASTQ or (multi-line) FASTA from a line iterator."""
    it = iter(lines)
    first = next(it, "")
    if not first:
        return
    if first.startswith("@"):
        header = first
        while header:
            seq = next(it, "").rstrip("\n")
            plus = next(it, "")
            next(it, "")
            if not plus.startswith("+"):
                raise ValueError("malformed FASTQ record")
            yield header[1:].split(None, 1)[0] if header[1:].strip() else "", seq
            header = next(it, "")
    elif first.startswith(">"):
        name, chunks = first[1:].split(None, 1)[0] if first[1:].strip() else "", []
        for line in it:
            if line.startswith(">"):
                yield name, "".join(chunks)
                name, chunks = line[1:].split(None, 1)[0] if line[1:].strip() else "", []
            else:
                chunks.append(line.rstrip("\n"))
        yield name, "".join(chunks)
    else:
        raise ValueError("not FASTQ/FASTA")


class ReadSource:
    """Iterator over (name, seq) for one input file; close() stops any subprocess."""

    def __init__(self, path, fmt=None, samtools="samtools", threads=2):
        self.path, self.fmt = path, fmt or detect_format(path)
        self._proc, self._fh = None, None
        if self.fmt == "bam":
            self._proc = subprocess.Popen(
                [*samtools.split(), "fasta", "-F", "0x900", "-@", str(max(1, threads - 1)), path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1 << 20)
            self._lines = self._proc.stdout
        elif self.fmt == "tar":
            self._lines = archive.open_text(path)
        elif self.fmt in ("fastq", "fasta"):
            with open(path, "rb") as fh:
                gz = fh.read(2) == b"\x1f\x8b"
            self._fh = (gzip.open if gz else open)(path, "rt", errors="replace")
            self._lines = self._fh
        else:
            raise ValueError(f"unsupported format {self.fmt}")

    def __iter__(self):
        return _fastx_lines(self._lines)

    def close(self, complete):
        """Stop reading. If the stream was read to the end, verify the producer succeeded."""
        if self._fh:
            self._fh.close()
        if self._proc:
            if not complete:
                self._proc.kill()
            self._proc.stdout.close()
            err = self._proc.stderr.read()
            rc = self._proc.wait()
            if complete and rc != 0:
                raise RuntimeError(f"samtools fasta failed (exit {rc}) on {self.path}: {err[-2000:]}")


# --------------------------------------------------------------------------
# Sampling helpers
# --------------------------------------------------------------------------

def _seeded(seed, salt):
    return zlib.crc32(f"{seed}:{salt}".encode())


def half_of(name, seed):
    """Replicate half for a read: stable for a given seed and read name."""
    return HALVES[zlib.crc32(name.encode(), _seeded(seed, "half")) & 1]


def keep_read(name, seed, fraction):
    if fraction >= 1:
        return True
    return (zlib.crc32(name.encode(), _seeded(seed, "subsample")) / 0xFFFFFFFF) < fraction


def _n50(lengths):
    total, acc = sum(lengths), 0
    for n in sorted(lengths, reverse=True):
        acc += n
        if acc * 2 >= total:
            return n
    return 0


def file_identity(path):
    st = os.stat(path)
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read(1 << 20))
    h.update(str(st.st_size).encode())
    return {"path": os.path.abspath(path), "size": st.st_size,
            "fingerprint": "sha256-head1MiB:" + h.hexdigest()[:32]}


# --------------------------------------------------------------------------
# Sketching
# --------------------------------------------------------------------------

def sketch(path, out_sig, out_json, label, tech, run_id, ksizes=(21,), scaled=2000,
           max_bases=0, read_fraction=1.0, min_read_len=0, seed=1, fmt=None,
           samtools="samtools", threads=2):
    """Sketch one read file. Writes a sourmash signature file and a JSON record."""
    import sourmash
    from sourmash import MinHash, SourmashSignature

    t0 = time.time()
    mhs = {(k, h): MinHash(n=0, ksize=k, scaled=scaled, track_abundance=True, seed=42)
           for k in ksizes for h in HALVES}
    src = ReadSource(path, fmt, samtools=samtools, threads=threads)
    stats = {h: {"reads": 0, "bases": 0} for h in HALVES}
    seen = skipped_short = skipped_subsample = 0
    used_bases, lengths, complete, finished = 0, [], True, False
    names_digest = hashlib.sha256()
    try:
        for name, seq in src:
            seen += 1
            if len(seq) < max(min_read_len, max(ksizes)):
                skipped_short += 1
                continue
            if not keep_read(name, seed, read_fraction):
                skipped_subsample += 1
                continue
            h = half_of(name, seed)
            for k in ksizes:
                mhs[(k, h)].add_sequence(seq, force=True)
            stats[h]["reads"] += 1
            stats[h]["bases"] += len(seq)
            used_bases += len(seq)
            lengths.append(len(seq))
            names_digest.update(name.encode() + b"\n")
            if max_bases and used_bases >= max_bases:
                complete = False
                break
        finished = True
    finally:
        # A read error leaves `finished` False: stop the producer without masking the error.
        src.close(complete and finished)

    sigs = [SourmashSignature(mhs[(k, h)], name=f"{label}|{tech}|{run_id}|k{k}|half{h}")
            for k in ksizes for h in HALVES]
    tmp_sig = out_sig + ".tmp"
    data = sourmash.save_signatures_to_json(sigs)
    with gzip.open(tmp_sig, "wb") as fh:
        fh.write(data if isinstance(data, bytes) else data.encode())
    os.replace(tmp_sig, out_sig)

    rec = {
        "label": label, "tech": tech, "run_id": run_id,
        "input": file_identity(path), "format": src.fmt,
        "signature": os.path.abspath(out_sig),
        "sourmash_version": sourmash.VERSION,
        "params": {"ksizes": list(ksizes), "scaled": scaled, "hash_seed": 42,
                   "max_bases": max_bases, "read_fraction": read_fraction,
                   "min_read_len": min_read_len, "sampling_seed": seed},
        "sampling": {
            "reads_seen": seen, "reads_used": sum(s["reads"] for s in stats.values()),
            "bases_used": used_bases, "skipped_short": skipped_short,
            "skipped_subsample": skipped_subsample,
            "stopped_at_max_bases": not complete,
            "read_len_mean": round(used_bases / len(lengths), 1) if lengths else 0,
            "read_len_n50": _n50(lengths),
            "halves": stats,
            "used_read_names_sha256": names_digest.hexdigest(),
        },
        "hashes": {f"k{k}_half{h}": len(mhs[(k, h)]) for k in ksizes for h in HALVES},
        "runtime_s": round(time.time() - t0, 1),
    }
    tmp = out_json + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, out_json)
    return rec
