"""Minimal PAF reader and interval helpers."""
import gzip


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_paf(path, min_len=0, min_mapq=0, primary_only=True):
    """Yield dicts for alignments with block length >= min_len and MAPQ >= min_mapq."""
    with open_text(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            tags = {t[:2]: t[5:] for t in p[12:]}
            if primary_only and tags.get("tp", "P") not in ("P",):
                continue
            a = {"q": p[0], "qlen": int(p[1]), "qs": int(p[2]), "qe": int(p[3]), "strand": p[4],
                 "t": p[5], "tlen": int(p[6]), "ts": int(p[7]), "te": int(p[8]),
                 "match": int(p[9]), "alen": int(p[10]), "mapq": int(p[11])}
            if "de" in tags:
                a["identity"] = 1 - float(tags["de"])  # gap-compressed divergence (needs -c)
            else:
                a["identity"] = a["match"] / a["alen"] if a["alen"] else 0.0
            if a["qe"] - a["qs"] >= min_len and a["mapq"] >= min_mapq:
                yield a


def merge_intervals(iv):
    """Union of half-open intervals [(s, e)], sorted."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def covered_bp(iv):
    return sum(e - s for s, e in merge_intervals(iv))


def subtract(a, b):
    """Parts of merged intervals a not covered by merged intervals b."""
    out, j = [], 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                out.append([cur, b[k][0]])
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            out.append([cur, e])
    return out
