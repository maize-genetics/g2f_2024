"""Genotype fingerprinting: which HiFi genotype is each query read set most consistent with?

Standard library only. Operates on FracMinHash sketches written by
g2f.fpsketch (sourmash signature JSON; hash -> abundance, per replicate half).

Method (docs/fingerprinting.md has the full rationale)
------------------------------------------------------
Reference panel: for every HiFi genotype G, the abundances of all its HiFi
runs are summed.
  * present(G): every hash seen at least once (includes HiFi error k-mers);
  * solid(G):  hashes whose abundance lies in [lo, hi], where lo is the first
    minimum of the abundance spectrum (end of the error peak, at least
    `min_count_floor`) and hi = `max_copy` x the main (single-copy) peak.
    This drops sequencing errors as well as high-copy repeats, organelle and
    high-coverage contaminant k-mers;
  * private(G): solid(G) hashes absent from present(G') of every other
    genotype G' (strict absence: not seen even once);
  * core: hashes solid in every genotype.

  * distinct(G vs H): solid(G) hashes absent from present(H): G's alleles
    (and G-specific sequence) where the two genotypes differ.

Query (an ONT run, a pool of runs, or a HiFi control) = the set of hashes it
contains. Sequencing errors in the query can only remove k-mers or add novel
ones; they cannot create a genotype's k-mers except by an error that happens
to reproduce another allele, which is measured by the empirical null.
  sensitivity s(Q)      = |core ∩ Q| / |core|   (coverage x per-k-mer accuracy)
  has(Q; G vs H)        = |distinct(G vs H) ∩ Q| / |distinct(G vs H)| / s(Q)
                          ~1 if Q is G; ~0 (error null) if Q is H
  private score(Q,G)    = |private(G) ∩ Q| / |private(G)| / s(Q)   (similarity
                          matrix: one number per query x genotype)
  containment(Q,G)      = |solid(G) ∩ Q| / |solid(G)| / s(Q)  (genome-wide
                          shared k-mer fraction; descriptive)
Q is consistent with G when, against every other genotype H, Q carries G's
distinct alleles (has(G vs H) high) and lacks H's (has(H vs G) low). Pairwise
comparisons also resolve close relatives that share most k-mers, and a query
carrying both G's and H's distinct alleles is reported as a possible mixture.

Empirical calibration (no fixed similarity threshold):
  * positive controls: HiFi runs scored against a panel in which their own
    genotype is represented only by OTHER runs (leave-one-run-out:
    independent SMRT cells), and replicate half A of each genotype scored
    against half B (same reads' library; technical replicate only);
  * negative controls: the same control queries scored against every other
    genotype (truth known), plus each ONT query's scores for every genotype
    except its best match (platform-matched null; label-free);
  * threshold tau (for pairwise scores) = midpoint between the largest
    negative and the smallest positive control score. If they overlap, the data cannot
    discriminate and every query is reported as ambiguous.
A query is then classified using a binomial (Wilson) interval of each score.
"""
import gzip
import json
import math
import statistics

CLASSES = ("supported", "ambiguous", "potentially_mismatched")


# --------------------------------------------------------------------------
# Loading sketches
# --------------------------------------------------------------------------

def load_sig(path, k):
    """Return {"A": {hash: abundance}, "B": {...}} for k-mer size k."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        data = json.load(fh)
    out = {}
    for entry in data:
        name = entry.get("name", "")
        half = name.rsplit("|half", 1)[1] if "|half" in name else "A"
        for s in entry["signatures"]:
            if s["ksize"] != k:
                continue
            ab = s.get("abundances") or [1] * len(s["mins"])
            d = out.setdefault(half, {})
            for h, a in zip(s["mins"], ab):
                d[h] = d.get(h, 0) + a
    if not out:
        raise ValueError(f"{path}: no sketch with k={k}")
    return out


def add_counts(dst, src):
    for h, a in src.items():
        dst[h] = dst.get(h, 0) + a
    return dst


# --------------------------------------------------------------------------
# Abundance spectrum -> solid band
# --------------------------------------------------------------------------

def spectrum(counts):
    hist = {}
    for a in counts.values():
        hist[a] = hist.get(a, 0) + 1
    return hist


def solid_band(hist, min_count_floor=2, max_copy=3.0, expected_depth=None):
    """Return dict(lo, peak, hi, note) from an abundance histogram.

    expected_depth (bases sketched / genome size) anchors the search for the
    single-copy peak to [0.4, 1.6] x the expected depth, so that at low depth
    a repeat peak is never mistaken for it."""
    if not hist:
        return {"lo": min_count_floor, "peak": None, "hi": None, "note": "empty"}
    top = max(hist)
    h = [hist.get(c, 0) for c in range(top + 2)]
    sm = [0.0] * len(h)
    for c in range(1, len(h)):
        win = [h[x] for x in (c - 1, c, c + 1) if 1 <= x < len(h)]
        sm[c] = sum(win) / len(win)
    valley = None
    for c in range(1, len(h) - 1):
        if sm[c] <= sm[c + 1]:
            valley = c
            break
    lo = max(min_count_floor, valley or min_count_floor)
    if lo > top:
        return {"lo": lo, "peak": None, "hi": None,
                "note": "no abundance peak: coverage too low to separate errors"}
    a, b = lo, top
    note = "" if valley else "no error/coverage valley found; floor used"
    if expected_depth:
        a = max(lo, int(0.4 * expected_depth))
        b = min(top, max(a, int(math.ceil(1.6 * expected_depth))))
        if valley and valley > 1.6 * expected_depth:
            # The "valley" lies beyond the expected peak: errors and the
            # single-copy peak are not separated at this depth.
            lo = a = max(min_count_floor, int(0.4 * expected_depth))
            note = "error valley not resolved at this depth; lower bound from expected depth"
        if expected_depth < 5:
            note = (note + "; " if note else "") + f"low depth ({expected_depth:.1f}x)"
    peak = max(range(a, b + 1), key=lambda c: (sm[c], -c))
    if peak <= lo:
        note = (note + "; " if note else "") + "main peak not above the error valley"
    hi = max(lo, int(math.floor(max_copy * peak)))
    return {"lo": lo, "peak": peak, "hi": hi, "note": note}


def solid_set(counts, band):
    lo, hi = band["lo"], band["hi"] or float("inf")
    return {h for h, a in counts.items() if lo <= a <= hi}


# --------------------------------------------------------------------------
# Panel
# --------------------------------------------------------------------------

class Reference:
    """One genotype's reference (full panel member or a control variant)."""

    def __init__(self, name, counts, params, runs, kind="full", bases=None):
        self.name, self.kind, self.runs = name, kind, list(runs)
        depth = bases / params["genome_size"] if bases and params.get("genome_size") else None
        self.expected_depth = depth
        self.band = solid_band(spectrum(counts), params["min_count_floor"], params["max_copy"],
                               depth)
        self.present = set(counts)
        self.solid = solid_set(counts, self.band)
        self.n_present = len(self.present)


def core_set(refs):
    sets = sorted((r.solid for r in refs.values()), key=len)
    if not sets:
        return set()
    core = set(sets[0])
    for s in sets[1:]:
        core &= s
    return core


class Panel:
    """References of all genotypes, partitioned by presence pattern.

    mask[h] (h in the union of solid sets) has bit i set when genotype i's
    reference contains h at all (present). part[G][m] = solid(G) hashes with
    mask m. Then, for genotypes G != H:
      private(G)       = part[G][bit(G)]
      distinct(G vs H) = union of part[G][m] over masks m without bit(H)
                         (G's solid k-mers never seen in H's reads)."""

    def __init__(self, refs, _mask=None):
        self.refs = refs
        self.genotypes = sorted(refs)
        self.bit = {g: 1 << i for i, g in enumerate(self.genotypes)}
        if _mask is None:
            universe = set()
            for g in self.genotypes:
                universe |= refs[g].solid
            pres = [(self.bit[g], refs[g].present) for g in self.genotypes]
            _mask = {h: sum(b for b, p in pres if h in p) for h in universe}
        self.mask = _mask
        self.part = {}
        for g in self.genotypes:
            d = {}
            for h in refs[g].solid:
                d.setdefault(self.mask[h], set()).add(h)
            self.part[g] = d
        self.private = {g: self.part[g].get(self.bit[g], set()) for g in self.genotypes}
        self.core = core_set(refs)
        self.universe = set(self.mask)
        self.distinct_n = {g: {h: sum(len(v) for m, v in self.part[g].items()
                                      if not m & self.bit[h])
                               for h in self.genotypes if h != g} for g in self.genotypes}

    def with_variant(self, genotype, variant):
        """Same panel with one genotype's reference replaced (for controls)."""
        refs = dict(self.refs)
        refs[genotype] = variant
        bit = self.bit[genotype]
        others = [(self.bit[g], refs[g].present) for g in self.genotypes if g != genotype]
        universe = set()
        for g in self.genotypes:
            universe |= refs[g].solid
        mask = {}
        for h in universe:
            m = self.mask.get(h)
            if m is None:
                m = sum(b for b, p in others if h in p)
            mask[h] = (m & ~bit) | (bit if h in variant.present else 0)
        return Panel(refs, _mask=mask)


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def wilson(hits, n, z):
    if n == 0:
        return 0.0, 1.0
    p = hits / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def _norm(hits, n, s, z):
    lo, hi = wilson(hits, n, z)
    if not s:
        return {"n": n, "hits": hits, "score": 0.0, "lo": 0.0, "hi": float("inf")}
    return {"n": n, "hits": hits, "score": (hits / n / s) if n else 0.0,
            "lo": lo / s, "hi": hi / s}


def score(panel, query, z):
    """Score a query hash set against every genotype (and genotype pair) of a panel."""
    q = query if isinstance(query, (set, frozenset)) else set(query)
    n_core = len(panel.core)
    core_hits = len(panel.core & q)
    s = core_hits / n_core if n_core else 0.0
    rows, pairs = {}, {}
    for g in panel.genotypes:
        hits_m = {m: len(v & q) for m, v in panel.part[g].items()}
        n_solid = sum(len(v) for v in panel.part[g].values())
        bg = panel.bit[g]
        p = _norm(hits_m.get(bg, 0), len(panel.private[g]), s, z)
        rows[g] = {"n_private": p["n"], "private_hits": p["hits"],
                   "recovery": p["hits"] / p["n"] if p["n"] else 0.0,
                   "score": p["score"], "score_lo": p["lo"], "score_hi": p["hi"],
                   "containment": (sum(hits_m.values()) / n_solid / s) if (s and n_solid) else 0.0}
        for h in panel.genotypes:
            if h == g:
                continue
            bh = panel.bit[h]
            hits = sum(c for m, c in hits_m.items() if not m & bh)
            pairs[(g, h)] = _norm(hits, panel.distinct_n[g][h], s, z)
    return {"sensitivity": s, "core_hits": core_hits, "n_core": n_core,
            "genotypes": rows, "pairs": pairs}


def support(result, g, min_n):
    """Pairwise evidence for genotype g in a query.

    has(g,h):   normalized recovery of g-vs-h distinct k-mers (query carries g's
                alleles where g and h differ); ~1 if the query is g.
    lacks(g,h): recovery of h-vs-g distinct k-mers; ~null if the query is g."""
    out = {"has_min": float("inf"), "has_min_lo": float("inf"), "lacks_max": 0.0,
           "lacks_max_hi": 0.0, "unresolved": [], "weakest": None, "strongest_other": None}
    for (a, b), v in result["pairs"].items():
        if a == g:
            if v["n"] < min_n or result["pairs"][(b, a)]["n"] < min_n:
                out["unresolved"].append(b)
                continue
            if v["score"] < out["has_min"]:
                out["has_min"], out["weakest"] = v["score"], b
            out["has_min_lo"] = min(out["has_min_lo"], v["lo"])
        elif b == g:
            if v["n"] < min_n or result["pairs"][(b, a)]["n"] < min_n:
                continue
            if v["score"] > out["lacks_max"]:
                out["lacks_max"], out["strongest_other"] = v["score"], a
            out["lacks_max_hi"] = max(out["lacks_max_hi"], v["hi"])
    if out["has_min"] == float("inf"):
        out["has_min"] = out["has_min_lo"] = 0.0
    out["net"] = out["has_min"] - out["lacks_max"]
    out["unresolved"] = sorted(set(out["unresolved"]))
    return out


def ranked(result, min_n):
    """Genotypes ordered by net pairwise support (best first)."""
    sup = {g: support(result, g, min_n) for g in result["genotypes"]}
    return sorted(sup.items(), key=lambda kv: -kv[1]["net"])


# --------------------------------------------------------------------------
# Calibration and classification
# --------------------------------------------------------------------------

def _pair_scores(result, g, min_n):
    """(positives, negatives) of pairwise scores if the query truly is g."""
    pos, neg = [], []
    for (a, b), v in result["pairs"].items():
        if v["n"] < min_n or result["pairs"][(b, a)]["n"] < min_n:
            continue
        if a == g:
            pos.append(v["score"])
        elif b == g:
            neg.append(v["score"])
    return pos, neg


def calibrate(controls, ont_results, min_sensitivity, min_n):
    """Decision threshold from empirical positive and negative pairwise scores.

    controls: [(truth, result)] (HiFi, truth known) define the threshold.
    ont_results: ONT run-level results, scored against their best-supported
    genotype (label-free). Their null is reported and compared with the control
    null, but not used for the threshold, because a genuinely mixed or
    contaminated ONT sample would inflate it."""
    pos, neg_ctrl, neg_ont, pos_ont = [], [], [], []
    for t, r in controls:
        if r["sensitivity"] < min_sensitivity or r["n_core"] < min_n or not r["sensitivity"]:
            continue
        p, n = _pair_scores(r, t, min_n)
        pos += p
        neg_ctrl += n
    for r in ont_results:
        if r["sensitivity"] < min_sensitivity or not r["sensitivity"] or r["n_core"] < min_n:
            continue
        best = ranked(r, min_n)[0][0]
        p, n = _pair_scores(r, best, min_n)
        pos_ont += p
        neg_ont += n
    neg = neg_ctrl
    med = statistics.median(neg) if neg else None
    out = {
        "n_positive": len(pos), "n_negative_controls": len(neg_ctrl),
        "n_negative_ont": len(neg_ont),
        "positive_min": min(pos) if pos else None,
        "positive_median": statistics.median(pos) if pos else None,
        "ont_best_positive_min": min(pos_ont) if pos_ont else None,
        "ont_best_positive_median": statistics.median(pos_ont) if pos_ont else None,
        "negative_max": max(neg) if neg else None,
        "negative_max_controls": max(neg_ctrl) if neg_ctrl else None,
        "negative_max_ont": max(neg_ont) if neg_ont else None,
        "negative_median": med,
        "negative_mad": statistics.median([abs(x - med) for x in neg]) if neg else None,
        "ont_negative_median": statistics.median(neg_ont) if neg_ont else None,
        "ont_negative_p95": (sorted(neg_ont)[int(0.95 * (len(neg_ont) - 1))]
                             if neg_ont else None),
    }
    if not pos or not neg:
        out.update(discriminating=False, threshold=None,
                   reason="no positive controls" if not pos else "no negative scores")
    elif out["positive_min"] <= out["negative_max"]:
        out.update(discriminating=False, threshold=None,
                   reason="positive and negative score distributions overlap")
    else:
        out.update(discriminating=True,
                   threshold=(out["positive_min"] + out["negative_max"]) / 2,
                   reason="positive controls separated from all negative scores")
    if pos and neg:
        out["separation"] = out["positive_min"] - out["negative_max"]
    return out


def classify(result, expected, calib, min_n):
    """Assignment of one query: supported | potentially_mismatched | ambiguous
    (| assigned, when no expected label is given)."""
    order = ranked(result, min_n)
    (best, bs), (second, ss) = order[0], (order[1] if len(order) > 1 else (None, None))
    rows = result["genotypes"]
    tau = calib.get("threshold")
    neg_med, neg_mad = calib.get("negative_median"), calib.get("negative_mad")
    es = support(result, expected, min_n) if expected in rows else None
    out = {
        "best": best, "best_net": bs["net"], "best_has_min": bs["has_min"],
        "best_lacks_max": bs["lacks_max"],
        "second": second, "second_net": ss["net"] if ss else None,
        "margin": bs["net"] - ss["net"] if ss else None,
        "best_private_score": rows[best]["score"],
        "expected_rank": [g for g, _ in order].index(expected) + 1 if es else None,
        "expected_has_min": es["has_min"] if es else None,
        "expected_lacks_max": es["lacks_max"] if es else None,
        "expected_private_score": rows[expected]["score"] if es else None,
        "expected_z": ((es["has_min"] - neg_med) / (1.4826 * neg_mad)
                       if es and neg_mad else None),
        "sensitivity": result["sensitivity"], "threshold": tau,
        "unresolved_pairs": ",".join(f"{best}~{u}" for u in bs["unresolved"]),
    }
    # Alleles of another genotype present clearly above the control null
    # (interval entirely above the largest negative control score).
    null_max = calib.get("negative_max")
    admix = [] if null_max is None or not result["sensitivity"] else sorted(
        h for (h, b), v in result["pairs"].items()
        if b == best and v["n"] >= min_n and result["pairs"][(b, h)]["n"] >= min_n
        and v["lo"] > null_max)
    out["possible_admixture"] = ",".join(admix)

    def consistent(sp):          # all pairwise intervals clear of the threshold
        return sp["has_min_lo"] > tau and sp["lacks_max_hi"] < tau

    def point_consistent(sp):
        return sp["has_min"] > tau and sp["lacks_max"] < tau

    if expected and expected not in rows:
        cls, why = "ambiguous", f"expected genotype {expected} is not in the HiFi panel"
    elif not calib.get("discriminating"):
        cls, why = "ambiguous", "data do not discriminate genotypes: " + calib.get("reason", "")
    elif result["n_core"] < min_n:
        cls, why = "ambiguous", (f"panel core too small ({result['n_core']} k-mers) to normalize "
                                 "scores; check the HiFi reference depths (panel.tsv)")
    elif not result["sensitivity"]:
        cls, why = "ambiguous", "insufficient data: no panel core k-mers found in the query"
    elif consistent(bs) and admix:
        cls, why = "ambiguous", (f"best match {best}, but the reads also carry alleles of "
                                 f"{','.join(admix)} above the control null (mixture or "
                                 "contamination?)")
    elif consistent(bs):
        if not expected:
            cls, why = "assigned", f"consistent only with {best} (no expected label)"
        elif best == expected:
            cls, why = "supported", f"consistent only with the expected genotype {best}"
        else:
            cls, why = "potentially_mismatched", (f"consistent only with {best}, not with the "
                                                  f"expected {expected}")
        if bs["unresolved"]:
            why += f"; cannot be distinguished from: {','.join(bs['unresolved'])}"
    else:
        mixed = admix if bs["has_min"] > tau else []
        excluded = bool(es) and es["has_min_lo"] is not None and \
            any(result["pairs"][(expected, h)]["hi"] < tau and result["pairs"][(h, expected)]["lo"] > tau
                for h in rows if h != expected
                and result["pairs"][(expected, h)]["n"] >= min_n
                and result["pairs"][(h, expected)]["n"] >= min_n)
        if mixed:
            why = (f"reads carry alleles of both {best} and {','.join(mixed)} "
                   "(mixture, contamination, or a genotype intermediate between them)")
        elif point_consistent(bs):
            why = f"insufficient data: best match {best} is not significant at this coverage"
        else:
            why = "no panel genotype is consistent with the reads (genotype absent from panel?)"
        if excluded and not mixed:
            cls = "potentially_mismatched"
            why = f"expected {expected} is contradicted by the reads; " + why
        else:
            cls = "ambiguous"
    out.update(classification=cls, reason=why)
    return out
