"""GFA segment / FASTA consistency checks for hifiasm contig outputs."""
import hashlib
import re

VALID_SEQ = re.compile(rb"^[ACGTNacgtn]*$")


def gfa_segments(path):
    """Yield (name, length, sha1) for each S line; checks LN:i: tag if present."""
    with open(path, "rb") as fh:
        for line in fh:
            if not line.startswith(b"S\t"):
                continue
            parts = line.rstrip(b"\n").split(b"\t")
            name, seq = parts[1].decode(), parts[2]
            if seq == b"*":
                raise ValueError(f"segment {name} has no sequence ('*'); need a full GFA")
            for tag in parts[3:]:
                if tag.startswith(b"LN:i:") and int(tag[5:]) != len(seq):
                    raise ValueError(f"segment {name}: LN tag {int(tag[5:])} != sequence length {len(seq)}")
            yield name, len(seq), hashlib.sha1(seq.upper()).hexdigest()


def fasta_records(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            line = line.rstrip(b"\n")
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode() if line[1:].strip() else "", []
            else:
                if name is None:
                    raise ValueError("FASTA does not start with a '>' header")
                chunks.append(line)
    if name is not None:
        yield name, b"".join(chunks)


def compare(gfa_path, fasta_path):
    """Return a validation dict; 'ok' is True only if everything matches."""
    errors = []
    gfa = list(gfa_segments(gfa_path))
    fa, names = [], set()
    for name, seq in fasta_records(fasta_path):
        if not name:
            errors.append("empty FASTA header")
        if name in names:
            errors.append(f"duplicate FASTA name {name}")
        names.add(name)
        if not seq:
            errors.append(f"empty sequence {name}")
        if not VALID_SEQ.match(seq):
            errors.append(f"non-ACGTN characters in {name}")
        fa.append((name, len(seq), hashlib.sha1(seq.upper()).hexdigest()))
    if len(gfa) != len(fa):
        errors.append(f"segment count {len(gfa)} != FASTA record count {len(fa)}")
    for g, f in zip(gfa, fa):
        if g != f:
            errors.append(f"mismatch: GFA {g[0]} ({g[1]} bp) vs FASTA {f[0]} ({f[1]} bp)")
            break
    return {
        "ok": not errors,
        "errors": errors[:20],
        "n_segments": len(gfa),
        "n_fasta_records": len(fa),
        "gfa_total_bp": sum(x[1] for x in gfa),
        "fasta_total_bp": sum(x[1] for x in fa),
    }
