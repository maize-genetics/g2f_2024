"""Aggregate per-read records from `seqkit fx2tab -n -l [-q]` into read QC.

One streaming pass gives yield, length distribution, read N50, quality
summary and an identity census (PacBio movies / ONT run IDs, duplicate read
IDs). Names may contain tabs (dorado FASTQ tags), so columns are parsed from
the right.
"""
import gzip
import re
from array import array

from .discovery import PACBIO_NAME_RE

LEN_BIN_EDGES_KB = [0, 1, 2, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 1000]
QUAL_BIN_EDGES = [0, 5, 7, 10, 12, 15, 20, 25, 30, 40]
RUNID_RE = re.compile(r"(?:^|\s)runid=(\S+)")
RG_RE = re.compile(r"(?:^|\s)RG:Z:([^_\s]+)")


def n50(sorted_desc, total):
    half, acc = total / 2, 0
    for i, n in enumerate(sorted_desc, 1):
        acc += n
        if acc >= half:
            return n, i
    return 0, 0


def _bin(value, edges):
    lo = edges[0]
    for hi in edges[1:]:
        if value < hi:
            return f"[{lo},{hi})"
        lo = hi
    return f"[{edges[-1]},inf)"


def aggregate(lines, tech, has_qual, genome_size=None, ids_out=None):
    lengths = array("Q")
    len_hist, q_hist = {}, {}
    census, kinds = {}, {}
    q_sum, q_n, q20, q30 = 0.0, 0, 0, 0
    seen_zmw, dup_ids, n_unparsed = {}, 0, 0
    seen_ids = set()
    idfh = gzip.open(ids_out, "wt") if ids_out else None
    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.rsplit("\t", 2 if has_qual else 1)
        name = parts[0]
        length = int(parts[1])
        lengths.append(length)
        b = _bin(length / 1000, LEN_BIN_EDGES_KB)
        len_hist[b] = len_hist.get(b, 0) + 1
        if has_qual:
            q = float(parts[2])
            q_sum += q
            q_n += 1
            q20 += q >= 20
            q30 += q >= 30
            qb = _bin(q, QUAL_BIN_EDGES)
            q_hist[qb] = q_hist.get(qb, 0) + 1
        rid = name.split(None, 1)[0] if name else ""
        if tech == "hifi":
            m = PACBIO_NAME_RE.match(name)
            if not m:
                n_unparsed += 1
                kinds["unparsed"] = kinds.get("unparsed", 0) + 1
                continue
            movie, zmw, kind = m.group(1), int(m.group(2)), m.group(3)
            kind = "subread" if "_" in kind else kind
            kinds[kind] = kinds.get(kind, 0) + 1
            census[movie] = census.get(movie, 0) + 1
            zs = seen_zmw.setdefault(movie, set())
            if zmw in zs:
                dup_ids += 1
            zs.add(zmw)
        else:
            m = RUNID_RE.search(name) or RG_RE.search(name)
            run = m.group(1) if m else "unknown"
            census[run] = census.get(run, 0) + 1
            key = rid.encode()
            if key in seen_ids:
                dup_ids += 1
            seen_ids.add(key)
            if idfh:
                idfh.write(rid + "\n")
    if idfh:
        idfh.close()
    srt = sorted(lengths, reverse=True)
    total = sum(srt)
    rn50, rl50 = n50(srt, total)
    out = {
        "tech": tech,
        "n_reads": len(srt),
        "total_bases": total,
        "min_len": srt[-1] if srt else 0,
        "max_len": srt[0] if srt else 0,
        "mean_len": round(total / len(srt), 1) if srt else 0,
        "read_n50": rn50,
        "read_l50": rl50,
        "length_hist_kb": {k: len_hist[k] for k in sorted(len_hist, key=_bin_key)},
        "duplicate_read_ids": dup_ids,
        "census": dict(sorted(census.items())),
        "census_kind": "pacbio_movie" if tech == "hifi" else "ont_runid",
    }
    if tech == "hifi":
        out["name_kinds"] = dict(sorted(kinds.items()))
        out["unparsed_names"] = n_unparsed
    if has_qual:
        out["mean_read_qual"] = round(q_sum / q_n, 3) if q_n else None
        out["frac_q20"] = round(q20 / q_n, 5) if q_n else None
        out["frac_q30"] = round(q30 / q_n, 5) if q_n else None
        out["qual_hist"] = {k: q_hist[k] for k in sorted(q_hist, key=_bin_key)}
        out["qual_definition"] = "seqkit fx2tab -q average read quality"
    if genome_size:
        out["genome_size_bp"] = genome_size
        out["approx_coverage"] = round(total / genome_size, 2)
    return out


def _bin_key(label):
    return float(label[1:].split(",")[0])
