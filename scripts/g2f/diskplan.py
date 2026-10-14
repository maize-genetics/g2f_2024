"""Disk projection and start gate for the full-scale launcher.

Estimates (bytes) per genotype come from its input sizes, using factors
measured on a ~2.3 Gb inbred maize genome (HiFi + ultra-long ONT):
  conv    converted HiFi FASTQ        1.1 x HiFi BAM bytes (0 for FASTQ/FASTA used in place)
  bins    hifiasm .bin intermediates  0.44 B per HiFi base + 13 GB (ONT index)
  mer     meryl DB + Merqury          0.30 B per HiFi base + 12 GB
  out     GFAs, FASTAs, BUSCO, logs   15 GB
  ont     filtered >= 50 kb ONT FASTQ 0.15 x ONT archive/FASTQ bytes (0 if used in place)
  pbi     PacBio index of each BAM    0.01 x HiFi BAM bytes
  peak         = conv + bins + mer + out + ont + pbi    (during hifiasm)
  persist_keep = conv + mer + out + ont + pbi           (after, prepared HiFi FASTQ kept)
  persist_temp = mer + out + ont + pbi                  (after, FASTQ deleted after use)

Gate for genotype g (never more than `jobs` genotypes running at once; all
running genotypes' remaining growth is reserved):
  committed = sum over running r of max(0, (peak_r if r not yet assembled
              else persist_r) - measured usage of results/<r>)
  available = free - committed - reserve
  start       if available >= peak_g
  keep FASTQ  (mode auto) if, after reserving peak_g, the space left still
              covers every not-yet-started genotype's persist_temp plus the
              (jobs - 1) largest of their transient needs (peak - persist_temp).
"""
import json
import os

GB = 1e9
BASES_PER_BYTE = {"bam": 84.3e9 / 50806416640, "fasta": 84.3e9 / 29188889732, "fastq": 1.5}


def estimate(inventory_record):
    hifi_bam = hifi_bases = ont_bytes = 0
    for i in inventory_record["inputs"]:
        if i["status"] != "use":
            continue
        if i["tech"] == "hifi":
            hifi_bases += i["size"] * BASES_PER_BYTE.get(i["format"], 1.5)
            if i["format"] == "bam":
                hifi_bam += i["size"]
        elif "ack:prefiltered_ont" not in i.get("notes", []):
            ont_bytes += i["size"]
    e = {"hifi_bam_bytes": hifi_bam, "hifi_bases_est": hifi_bases,
         "conv": 1.1 * hifi_bam, "bins": 0.44 * hifi_bases + 13 * GB,
         "mer": 0.30 * hifi_bases + 12 * GB, "out": 15 * GB, "ont": 0.15 * ont_bytes,
         "pbi": 0.01 * hifi_bam}
    e["persist_temp"] = e["mer"] + e["out"] + e["ont"] + e["pbi"]
    e["persist_keep"] = e["persist_temp"] + e["conv"]
    e["peak"] = e["persist_keep"] + e["bins"]
    return {k: int(v) for k, v in e.items()}


def decide(sample, est, free, running, pending, reserve, jobs, mode="auto", own=None):
    """running: {s: {"used": bytes, "assembled": bool, "retain": bool}}; pending: [s].
    own: {"used", "assembled"} of this genotype's existing results (a resumed run
    only needs what is not there yet)."""
    committed = 0
    for r, st in running.items():
        e = est[r]
        target = e["peak"] if not st.get("assembled") else (
            e["persist_keep"] if st.get("retain", True) else e["persist_temp"])
        committed += max(0, target - st.get("used", 0))
    avail = free - committed - reserve
    g = dict(est[sample])
    if own:
        target = g["persist_keep"] if own.get("assembled") else g["peak"]
        g["peak"] = max(0, target - own.get("used", 0))
    fut = [est[f] for f in pending if f != sample]
    trans = sorted((f["peak"] - f["persist_temp"] for f in fut), reverse=True)
    future = sum(f["persist_temp"] for f in fut) + sum(trans[:max(0, jobs - 1)])
    start = avail >= g["peak"]
    if mode == "keep":
        retain = True
    elif mode == "temp":
        retain = False
    else:
        retain = start and avail - g["peak"] >= future
    return {"sample": sample, "start": start, "retain_prepared_hifi": retain,
            "free_gb": round(free / GB, 1), "committed_gb": round(committed / GB, 1),
            "reserve_gb": round(reserve / GB, 1), "available_gb": round(avail / GB, 1),
            "peak_gb": round(g["peak"] / GB, 1), "future_need_gb": round(future / GB, 1),
            "mode": mode}


def du_bytes(path):
    """Allocated bytes under path (what `df` counts), without following links."""
    total = 0
    for root, dirs, files in os.walk(path):
        for f in files + dirs:
            try:
                st = os.lstat(os.path.join(root, f))
            except OSError:
                continue
            total += getattr(st, "st_blocks", 0) * 512
    return total


def assembled(results_dir, sample):
    adir = os.path.join(results_dir, sample, "assembly")
    if not os.path.isdir(adir):
        return False
    return any(os.path.exists(os.path.join(adir, m, f"{sample}.asm.bp.p_ctg.gfa"))
               for m in os.listdir(adir))


def gate(state_dir, sample, results_dir, reserve, jobs, mode, free):
    """Read launcher state files, decide, and record the decision if starting."""
    with open(os.path.join(state_dir, "estimates.json")) as fh:
        est = json.load(fh)
    running, started = {}, set()
    for f in os.listdir(state_dir):
        s, _, kind = f.rpartition(".")
        if kind in ("running", "done", "failed"):
            started.add(s)
        if kind == "running" and s != sample:
            with open(os.path.join(state_dir, f)) as fh:
                st = json.load(fh)
            running[s] = {"used": du_bytes(os.path.join(results_dir, s)),
                          "assembled": assembled(results_dir, s),
                          "retain": st.get("retain_prepared_hifi", True)}
    pending = [s for s in est if s not in started]
    own = {"used": du_bytes(os.path.join(results_dir, sample)),
           "assembled": assembled(results_dir, sample)}
    d = decide(sample, est, free, running, pending, reserve, jobs, mode, own)
    d["running"] = sorted(running)
    if len(running) >= jobs:
        # Hard concurrency limit, independent of how the launcher schedules:
        # never more than `jobs` genotypes running (1 = strictly sequential).
        d["start"] = False
        d["reason"] = f"{len(running)} genotype(s) already running, limit {jobs}"
    if d["start"]:
        with open(os.path.join(state_dir, f"{sample}.running"), "w") as fh:
            json.dump(d, fh, indent=2)
    return d


def has_checkpoint(results_dir, sample):
    """hifiasm error-correction checkpoint (<sample>.asm.ec.bin) left by an
    interrupted attempt in any method directory."""
    adir = os.path.join(results_dir, sample, "assembly")
    if not os.path.isdir(adir):
        return False
    return any(os.path.exists(os.path.join(adir, m, f"{sample}.asm.ec.bin"))
               for m in os.listdir(adir))


def start_order(est, results_dir):
    """Sequential start order that finishes assemblies soonest and exposes the
    least work to an outage: (0) genotypes already assembled (QC only), then
    (1) those with a hifiasm checkpoint to resume from, then (2) the rest by
    estimated HiFi bases, smallest (shortest assembly) first."""
    def key(s):
        rank = 0 if assembled(results_dir, s) else 1 if has_checkpoint(results_dir, s) else 2
        return (rank, est[s]["hifi_bases_est"], s)
    return sorted(est, key=key)
