#!/usr/bin/env python3
"""Telomere identification and telomere-based assembly-integrity QC for contig and scaffold assemblies.

Subcommands (each writes to a new directory, built as <dir>.partial and renamed after its outputs validate; never
overwrites):
  check-env  report tidk / minimap2 / samtools versions (tidk >= 0.3.0 is required for `tidk ends`)
  scan       one FASTA: tidk ends (primary call + sensitivity grid), tidk search (window counts) and a position-level
             array caller (every telomeric-repeat array, both strands, with coordinates); records, gaps, per-record MD5
  integrate  contigs + scaffolds of one assembly: end classes, contig arrays lifted through the AGP and checked against the
             independent scaffold scan, chromosome-end inventory, telomeres placed at internal joins, interstitial arrays
  tips       end identity against a reference (ideally telomere-to-telomere): each telomere-bearing contig end and each
             non-telomeric chromosome end, oriented end-first, aligned with minimap2; outer and inner window calls
  closure    existing long-read alignments to the contigs: do reads stop at a telomere-bearing end, and is any sequence
             clipped beyond it itself telomeric?
  compare    two scaffold scans of the same assembly (e.g. before/after a scaffold revision): end classes and internal arrays
  plot       chromosome maps (one per integrate directory) and a cross-assembly end heatmap (needs matplotlib)

Conventions (as tidk): a telomere reads C-rich (e.g. CCCTAAA) along the forward strand at a sequence START and G-rich
(TTTAGGG) at a sequence END. The motif is a scientific parameter (config/telomere_science.yaml; default the plant
TTTAGGG). Telomeric repeats inside a chromosome are evidence to examine, never a misassembly call on their own.

  g2f_telomere.py check-env
  g2f_telomere.py scan      --fasta X.fa[.gz] --label NAME --out DIR [--science YAML] [--tidk PATH] [--explore]
  g2f_telomere.py integrate --contigs-scan DIR --scaffolds-scan DIR --agp X.agp --label NAME --out DIR
                            [--join-table TSV] [--split-ends TSV] [--science YAML]
  g2f_telomere.py tips      --integrated DIR --contigs FA --scaffolds FA --reference REF.fa --out DIR [--threads N]
  g2f_telomere.py closure   --integrated DIR --alignments X.cram|bam --contigs FA --out DIR [--reference-fasta FA]
  g2f_telomere.py compare   --before SCAN_DIR --after SCAN_DIR --label NAME --out DIR
  g2f_telomere.py plot      --integrated DIR [DIR ...] --out DIR
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

SCHEMA = "g2f-telomere/1"
with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()
RC = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")

# Defaults = config/telomere_science.yaml (tests/test_telomere.py keeps them equal).
DEFAULTS = {
    "version": 1,
    "motif": "TTTAGGG",
    "tidk": {
        "primary": {"window": 10000, "min_length": 200, "max_offset": 1000, "exact": False},
        "length_window": 100000,
        "sensitivity": [
            {"name": "l500", "window": 100000, "min_length": 500, "max_offset": 1000, "exact": False},
            {"name": "l1000", "window": 100000, "min_length": 1000, "max_offset": 1000, "exact": False},
            {"name": "l3000", "window": 100000, "min_length": 3000, "max_offset": 1000, "exact": False},
            {"name": "o10k", "window": 100000, "min_length": 200, "max_offset": 10000, "exact": False},
            {"name": "exact", "window": 10000, "min_length": 200, "max_offset": 1000, "exact": True},
        ],
        "search_window": 10000,
    },
    "arrays": {"max_gap": 50, "min_hits": 5},
    "classes": {"kb_bp": 1000, "short_min_bp": 200, "terminal_window": 100000, "at_end": 1000},
    "chromosome_regex": r"_chr(\d+)(?:_part(\d+))?$",
    "tips": {"segment_bp": 1000000, "outer_bp": 300000, "inner_bp": [300000, 1000000], "near_end_bp": 4000000,
             "min_mapq": 5, "preset": "asm20"},
    "closure": {"min_mapq": 20, "outer_bp": 1000, "min_clip_bp": 500},
}


# ---------------------------------------------------------------- generic helpers
def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg):
    print(f"[{now()}] {msg}", file=sys.stderr, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def opener(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def write_tsv(path, rows, keys=None):
    if keys is None:
        keys = []
        for r in rows:
            keys += [k for k in r if k not in keys]
    with open(path, "w") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def revcomp(s):
    return s.translate(RC)[::-1]


def load_science(path):
    sci = json.loads(json.dumps(DEFAULTS))
    if path:
        import yaml
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}

        def merge(a, b):
            for k, v in b.items():
                if isinstance(v, dict) and isinstance(a.get(k), dict):
                    merge(a[k], v)
                else:
                    a[k] = v
        merge(sci, user)
    return sci


class Staged:
    """Write into <out>.partial; rename to <out> only via commit(). Refuses an existing <out> or <out>.partial."""

    def __init__(self, out):
        self.out = os.path.abspath(out)
        self.partial = self.out + ".partial"
        for p in (self.out, self.partial):
            if os.path.exists(p):
                sys.exit(f"refusing to overwrite {p} (inspect it; nothing was changed)")
        os.makedirs(self.partial)

    def path(self, *parts):
        p = os.path.join(self.partial, *parts)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        return p

    def dir(self, *parts):
        p = os.path.join(self.partial, *parts)
        os.makedirs(p, exist_ok=True)
        return p

    def commit(self, record):
        files = {}
        for d, _, fs in os.walk(self.partial):
            for f in sorted(fs):
                p = os.path.join(d, f)
                files[os.path.relpath(p, self.partial)] = sha256_file(p)
        record.update(schema=SCHEMA, script=os.path.abspath(__file__), script_sha256=SCRIPT_SHA256,
                      python=sys.version.split()[0], validated=now(), outputs_sha256=files)
        with open(os.path.join(self.partial, "VALIDATED.json"), "w") as fh:
            json.dump(record, fh, indent=1)
        os.rename(self.partial, self.out)
        log(f"validated: {self.out}")


def require_validated(d):
    v = os.path.join(d, "VALIDATED.json")
    if not os.path.exists(v):
        sys.exit(f"{d}: no VALIDATED.json; run the producing step first")
    with open(v) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- tools
def tidk_version(tidk):
    try:
        out = subprocess.run([tidk, "--version"], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return None
    m = re.search(r"tidk (\d+)\.(\d+)\.(\d+)", out)
    return tuple(int(x) for x in m.groups()) if m else None


def find_tool(explicit, env, name):
    return explicit or os.environ.get(env) or shutil.which(name) or name


def cmd_check_env(a):
    tidk = find_tool(a.tidk, "G2F_TIDK", "tidk")
    tv = tidk_version(tidk)
    print(f"tidk:     {tidk}: {'.'.join(map(str, tv)) if tv else 'NOT FOUND'}"
          + ("" if (tv and tv >= (0, 3, 0)) else "  (scan needs tidk >= 0.3.0 for `tidk ends`)"))
    for name, env, flag in (("minimap2", "G2F_MINIMAP2", "--version"), ("samtools", "G2F_SAMTOOLS", "--version")):
        exe = find_tool(None, env, name)
        try:
            v = subprocess.run([exe, flag], capture_output=True, text=True).stdout.splitlines()[0]
        except (FileNotFoundError, IndexError):
            v = "NOT FOUND"
        print(f"{name:9s} {exe}: {v}  ({'tips' if name == 'minimap2' else 'tips, closure'} only)")
    return 0 if tv and tv >= (0, 3, 0) else 1


# ---------------------------------------------------------------- array caller
def chains(seq, motif, max_gap, min_hits):
    """Exact, non-overlapping hits of `motif` chained while the gap between hits is <= max_gap. -> [(start0, end, hits)]"""
    out, cur, m = [], None, len(motif)
    for x in re.finditer(re.escape(motif), seq):
        p = x.start()
        if cur is not None and p - cur[1] <= max_gap:
            cur[1] = p + m
            cur[2] += 1
        else:
            if cur is not None and cur[2] >= min_hits:
                out.append(tuple(cur))
            cur = [p, p + m, 1]
    if cur is not None and cur[2] >= min_hits:
        out.append(tuple(cur))
    return out


def fasta_records(path, digest):
    name, buf = None, []
    with opener(path) as fh:
        for line in fh:
            digest.update(line)
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(buf).decode().upper()
                name, buf = line[1:].split()[0].decode(), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, b"".join(buf).decode().upper()


def scan_arrays(fasta, outdir, label, motif, max_gap, min_hits):
    motif = motif.upper()
    rc = revcomp(motif)
    if motif == rc:
        sys.exit("motif is its own reverse complement; strands cannot be distinguished")
    sha = hashlib.sha256()
    n = bp = 0
    hits = {"G-rich": 0, "C-rich": 0}
    with open(os.path.join(outdir, f"{label}.arrays.tsv"), "w") as fa, \
            open(os.path.join(outdir, f"{label}.gaps.tsv"), "w") as fg, \
            open(os.path.join(outdir, f"{label}.records.tsv"), "w") as fr:
        fa.write("seq\tseq_length\tstart0\tend\tspan\thits\tmotif_bp\tdensity\tstrand\tdist_start\tdist_end\n")
        fg.write("seq\tstart0\tend\tlength\n")
        fr.write("seq\tlength\tmd5_upper\tn_bp\thits_G\thits_C\n")
        for name, seq in fasta_records(fasta, sha):
            L = len(seq)
            n += 1
            bp += L
            hg, hc = seq.count(motif), seq.count(rc)
            hits["G-rich"] += hg
            hits["C-rich"] += hc
            arr = [(s, e, h, "G-rich") for s, e, h in chains(seq, motif, max_gap, min_hits)]
            arr += [(s, e, h, "C-rich") for s, e, h in chains(seq, rc, max_gap, min_hits)]
            for s, e, h, st in sorted(arr):
                fa.write(f"{name}\t{L}\t{s}\t{e}\t{e - s}\t{h}\t{len(motif) * h}\t{len(motif) * h / (e - s):.3f}\t{st}\t{s}\t{L - e}\n")
            for x in re.finditer("N{10,}", seq):
                fg.write(f"{name}\t{x.start()}\t{x.end()}\t{x.end() - x.start()}\n")
            fr.write(f"{name}\t{L}\t{hashlib.md5(seq.encode()).hexdigest()}\t{seq.count('N')}\t{hg}\t{hc}\n")
    return {"records": n, "bp": bp, "sha256_decompressed": sha.hexdigest(), "exact_hits": hits}


# ---------------------------------------------------------------- scan
def tidk_ends(tidk, motif, fasta, outdir, name, cfg):
    cmd = [tidk, "ends", "-s", motif, "-w", str(cfg["window"]), "--min-length", str(cfg["min_length"]),
           "--max-offset", str(cfg["max_offset"]), "-o", name, "-d", outdir, fasta]
    if cfg.get("exact"):
        cmd.insert(-1, "--exact")
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"tidk ends failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return " ".join(cmd)


def cmd_scan(a):
    sci = load_science(a.science)
    tidk = find_tool(a.tidk, "G2F_TIDK", "tidk")
    tv = tidk_version(tidk)
    if not tv or tv < (0, 3, 0):
        sys.exit(f"tidk >= 0.3.0 required (found {tv} at {tidk}); `tidk ends` does not exist before 0.3.0")
    if not os.path.exists(a.fasta):
        sys.exit(f"missing FASTA {a.fasta}")
    fsha = sha256_file(a.fasta)
    if a.expected_sha256 and fsha != a.expected_sha256:
        sys.exit(f"{a.fasta}: sha256 {fsha} != expected {a.expected_sha256}")
    st = Staged(a.out)
    motif = sci["motif"].upper()
    cmds = []
    t = sci["tidk"]
    log(f"{a.label}: tidk ends (primary + {len(t['sensitivity'])} sensitivity settings)")
    cmds.append(tidk_ends(tidk, motif, a.fasta, st.dir("ends"), "primary", t["primary"]))
    lw = dict(t["primary"], window=t["length_window"])
    cmds.append(tidk_ends(tidk, motif, a.fasta, st.dir("ends"), "length", lw))
    for s in t["sensitivity"]:
        cmds.append(tidk_ends(tidk, motif, a.fasta, st.dir("ends"), s["name"], s))
    log(f"{a.label}: tidk search")
    scmd = [tidk, "search", "-s", motif, "-w", str(t["search_window"]), "-o", "search", "-d", st.dir("search"), a.fasta]
    p = subprocess.run(scmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"tidk search failed: {p.stderr[-2000:]}")
    cmds.append(" ".join(scmd))
    win = st.path("search", "search_telomeric_repeat_windows.tsv")
    with open(win, "rb") as fi, gzip.open(win + ".gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(win)                                   # our own uncompressed intermediate inside the .partial directory
    if a.explore:
        log(f"{a.label}: tidk explore")
        with open(st.path("explore.tsv"), "w") as fo:
            ecmd = [tidk, "explore", "-m", "5", "-x", "12", a.fasta]
            subprocess.run(ecmd, stdout=fo, stderr=subprocess.DEVNULL, check=True)
        cmds.append(" ".join(ecmd) + " > explore.tsv")
    log(f"{a.label}: position-level array scan")
    arr = scan_arrays(a.fasta, st.dir("arrays"), a.label, motif, sci["arrays"]["max_gap"], sci["arrays"]["min_hits"])
    # validation: tidk saw every record; tidk search exact counts agree with the scanner up to window-straddling hits
    prim = read_tsv(st.path("ends", "primary.ends.tsv"))
    if len(prim) != 2 * arr["records"]:
        sys.exit(f"tidk ends reported {len(prim)} ends for {arr['records']} records")
    fwd = rev = 0
    with gzip.open(win + ".gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fwd += int(r["forward_repeat_number"])
            rev += int(r["reverse_repeat_number"])
    tol = max(10, 0.01 * (arr["exact_hits"]["G-rich"] + arr["exact_hits"]["C-rich"]))
    if not (0 <= arr["exact_hits"]["G-rich"] - fwd <= tol and 0 <= arr["exact_hits"]["C-rich"] - rev <= tol):
        sys.exit(f"array scanner and tidk search disagree beyond window-boundary losses: scanner {arr['exact_hits']}, "
                 f"tidk forward {fwd} reverse {rev}")
    st.commit({"step": "scan", "label": a.label, "fasta": os.path.abspath(a.fasta), "fasta_sha256": fsha,
               "tidk": tidk, "tidk_version": ".".join(map(str, tv)), "science": sci, "commands": cmds,
               "arrays": arr, "tidk_search_totals": {"forward": fwd, "reverse": rev}})
    return 0


# ---------------------------------------------------------------- integrate
def scan_tables(d):
    rec = require_validated(d)
    ends = {}
    for f in os.listdir(os.path.join(d, "ends")):
        if f.endswith(".ends.tsv"):
            ends[f[:-9]] = {(r["id"], "L" if r["side"] == "start" else "R"): r for r in read_tsv(os.path.join(d, "ends", f))}
    lab = rec["label"]
    arrays = collections.defaultdict(list)
    for r in read_tsv(os.path.join(d, "arrays", f"{lab}.arrays.tsv")):
        for k in ("seq_length", "start0", "end", "span", "hits", "motif_bp", "dist_start", "dist_end"):
            r[k] = int(r[k])
        r["density"] = float(r["density"])
        arrays[r["seq"]].append(r)
    records = {r["seq"]: r for r in read_tsv(os.path.join(d, "arrays", f"{lab}.records.tsv"))}
    return rec, ends, arrays, records


def tier(bp, sci):
    c = sci["classes"]
    return "kb" if bp >= c["kb_bp"] else ("short" if bp >= c["short_min_bp"] else "trace")


def classify_end(seq, L, side, ends, arrays, sci):
    c = sci["classes"]
    exp = "C-rich" if side == "L" else "G-rich"
    p, ln = ends["primary"][(seq, side)], ends["length"][(seq, side)]
    dist = (lambda r: r["dist_start"]) if side == "L" else (lambda r: r["dist_end"])
    near = [r for r in arrays.get(seq, []) if dist(r) <= min(c["terminal_window"], L // 2) and r["motif_bp"] >= c["short_min_bp"]]
    near_exp = [r for r in near if r["strand"] == exp]
    near_opp = [r for r in near if r["strand"] != exp]
    length = 0
    if p["status"] == "present":
        length = int(ln["telomere_length"]) if ln["status"] == "present" else int(p["telomere_length"])
        cls = "TEL_KB" if length >= c["kb_bp"] else "TEL_SHORT"
    elif near_exp:
        cls = "SUBTERMINAL"
    elif p["status"] == "wrong_strand" or near_opp:
        cls = "WRONG_STRAND"
    else:
        cls = "ABSENT"
    best = lambda xs: max(xs, key=lambda r: r["motif_bp"]) if xs else None
    fmt = lambda r: f"{r['start0'] + 1}-{r['end']} ({r['motif_bp']} bp, {dist(r)} bp from end)" if r else ""
    out = {"class": cls, "tidk_status": p["status"], "telomere_length": length if p["status"] == "present" else "",
           "offset": p["offset"] if p["status"] == "present" else "",
           "window_limited": ln["window_limited"] if ln["status"] == "present" else "",
           "expected_strand": exp, "expected_strand_array": fmt(best(near_exp)), "opposite_strand_array": fmt(best(near_opp))}
    for k in ends:
        if k not in ("primary", "length"):
            out[f"tidk_{k}"] = ends[k][(seq, side)]["status"]
    return out


def read_agp(path):
    objs = collections.OrderedDict()
    with open(path) as fh:
        for l in fh:
            if l.startswith("#") or not l.strip():
                continue
            x = l.rstrip("\n").split("\t")
            if x[4] in ("N", "U"):
                objs.setdefault(x[0], []).append({"s": int(x[1]), "e": int(x[2]), "gap": True})
            else:
                objs.setdefault(x[0], []).append({"obj": x[0], "s": int(x[1]), "e": int(x[2]), "ctg": x[5],
                                                  "cs": int(x[6]), "ce": int(x[7]), "o": "-" if x[8] == "-" else "+"})
    return objs


def lift(c, s0, e):
    """component 0-based half-open -> object 0-based half-open (component used from cs..ce)"""
    a, b = s0 - (c["cs"] - 1), e - (c["cs"] - 1)
    if c["o"] == "+":
        return c["s"] - 1 + a, c["s"] - 1 + b
    return c["e"] - b, c["e"] - a


def cmd_integrate(a):
    sci = load_science(a.science)
    crec, cends, carr, crecords = scan_tables(a.contigs_scan)
    srec, sends, sarr, srecords = scan_tables(a.scaffolds_scan)
    objs = read_agp(a.agp)
    split = {}
    if a.split_ends:
        for r in read_tsv(a.split_ends):
            split[(r["contig"], r["side"])] = r.get("note", "curation split")
    joins = {}
    if a.join_table:
        for r in read_tsv(a.join_table):
            joins[(r["contig_A"], r["contig_B"])] = joins[(r["contig_B"], r["contig_A"])] = r
    comp_of = {}
    for o, cl in objs.items():
        for i, c in enumerate(cl):
            if not c.get("gap"):
                comp_of[c["ctg"]] = (o, i, cl)
    missing = set(crecords) - set(comp_of)
    if missing:
        sys.exit(f"{len(missing)} contigs are not in the AGP (e.g. {sorted(missing)[:3]}); every contig must be placed or unplaced")
    if set(objs) != set(srecords):
        sys.exit("AGP objects differ from the scaffold FASTA records")
    st = Staged(a.out)
    cr = re.compile(sci["chromosome_regex"])         # chromosome objects; every other object is an unplaced record

    def join_info(ctg, other):
        j = joins.get((ctg, other)) if other else None
        return (f"{j.get('class', '')}:{j.get('evidence_call', '')}" if j else ("" if not other else "no join record"))

    # 1. contig ends
    cend = {}
    for seq, r in crecords.items():
        L = int(r["length"])
        o, i, cl = comp_of[seq]
        placed = cr.search(o) is not None            # on a chromosome object; otherwise its own unplaced record
        for side in "LR":
            d = classify_end(seq, L, side, cends, carr, sci)
            d.update(contig=seq, contig_length=L, side=side, split_end=split.get((seq, side), ""),
                     object=o, placed=placed)
            cend[(seq, side)] = d
    # 2. scaffold ends and chromosome inventory
    parts = collections.defaultdict(list)
    for o in objs:
        m = cr.search(o)
        if m:
            parts[int(m.group(1))].append((int(m.group(2) or 0), o))
    send, inv = [], []
    for chrn in sorted(parts):
        olist = [o for _, o in sorted(parts[chrn])]
        for pi, o in enumerate(olist):
            L = int(srecords[o]["length"])
            cl = [c for c in objs[o] if not c.get("gap")]
            for side in "LR":
                d = classify_end(o, L, side, sends, sarr, sci)
                c = cl[0] if side == "L" else cl[-1]
                out_side = ("L" if c["o"] == "+" else "R") if side == "L" else ("R" if c["o"] == "+" else "L")
                in_side = "R" if out_side == "L" else "L"
                nb = (cl[1] if side == "L" else cl[-2]) if len(cl) > 1 else None
                ce, ci = cend[(c["ctg"], out_side)], cend[(c["ctg"], in_side)]
                d.update(chromosome=f"chr{chrn:02d}", end="start" if side == "L" else "end",
                         is_chromosome_end=(side == "L" and pi == 0) or (side == "R" and pi == len(olist) - 1),
                         scaffold=o, scaffold_length=L, terminal_contig=c["ctg"], terminal_contig_orientation=c["o"],
                         terminal_contig_length=c["ce"] - c["cs"] + 1, outward_contig_end=out_side,
                         contig_end_class=ce["class"], consistent_with_contig=ce["class"] == d["class"],
                         inward_contig_end_class=ci["class"], split_end=ce["split_end"],
                         first_join=f"{c['ctg']}|{nb['ctg']}" if nb else "", first_join_class=join_info(c["ctg"], nb and nb["ctg"]))
                send.append(d)
    # 3. contig arrays >= short_min_bp lifted to scaffold coordinates and checked against the scaffold scan
    sidx = {(r["seq"], r["start0"], r["end"], r["strand"]) for v in sarr.values() for r in v}
    rows, lifted, found = [], set(), 0
    c_ = sci["classes"]
    for seq, alist in carr.items():
        for r in alist:
            if r["motif_bp"] < c_["short_min_bp"]:
                continue
            o, i, cl = comp_of[seq]
            c = cl[i]
            dmin = min(r["dist_start"], r["dist_end"])
            near = "L" if r["dist_start"] <= r["dist_end"] else "R"
            pos = "terminal" if dmin <= c_["at_end"] else ("subterminal" if dmin <= c_["terminal_window"] else "interstitial")
            s0, e = lift(c, r["start0"], r["end"])
            strand = r["strand"] if c["o"] == "+" else ("G-rich" if r["strand"] == "C-rich" else "C-rich")
            key = (o, s0, e, strand)
            lifted.add(key)
            found += key in sidx
            SL = int(srecords[o]["length"])
            scaf_side = near if c["o"] == "+" else ("R" if near == "L" else "L")
            nb = None
            if scaf_side == "L":
                nb = next((x for x in reversed(cl[:i]) if not x.get("gap")), None)
            else:
                nb = next((x for x in cl[i + 1:] if not x.get("gap")), None)
            gaps = [g for g in cl if g.get("gap")]
            gd = min([min(abs(g["s"] - 1 - e), abs(s0 - g["e"])) for g in gaps], default="")
            rows.append(dict(contig=seq, contig_length=r["seq_length"], contig_start=r["start0"] + 1, contig_end=r["end"],
                             motif_bp=r["motif_bp"], span=r["span"], density=r["density"], tier=tier(r["motif_bp"], sci),
                             contig_strand=r["strand"], position_in_contig=pos, dist_to_contig_end=dmin, nearest_contig_end=near,
                             expected_strand_for_that_end=r["strand"] == ("C-rich" if near == "L" else "G-rich"),
                             split_end=split.get((seq, near), "") if pos != "interstitial" else "",
                             scaffold=o, scaffold_start=s0 + 1, scaffold_end=e, scaffold_strand=strand,
                             dist_to_scaffold_end=min(s0, SL - e), dist_to_nearest_gap=gd, contig_orientation=c["o"],
                             faces="scaffold_end" if nb is None else "join", neighbour_across_join=nb["ctg"] if nb else "",
                             join_class=join_info(seq, nb and nb["ctg"]) if nb else "",
                             object_is_chromosome=cr.search(o) is not None))
    orphans = [k for k, v in sarr.items() for r in v if r["motif_bp"] >= c_["short_min_bp"]
               and (k, r["start0"], r["end"], r["strand"]) not in lifted]
    mism = [d for d in send if not d["consistent_with_contig"]]
    checks = {"contig_arrays_lifted": len(rows), "found_identically_in_scaffold_scan": found,
              "scaffold_arrays_without_contig_source": len(orphans), "scaffold_vs_contig_end_class_mismatches": len(mism)}
    if found != len(rows) or orphans or mism:
        sys.exit(f"self-check failed (AGP does not describe these FASTAs, or a scan is stale): {checks}")
    candidates = [r for r in rows if r["position_in_contig"] == "interstitial"
                  or (r["faces"] == "join" and r["tier"] == "kb")]
    chrom_ends = [d for d in send if d["is_chromosome_end"]]
    cnt = collections.Counter(d["class"] for d in chrom_ends)
    tel = [d for d in cend.values() if d["class"] in ("TEL_KB", "TEL_SHORT")]
    both = collections.Counter(d["contig"] for d in tel)
    summary = {"label": a.label, "chromosome_ends": len(chrom_ends), "telomere_kb": cnt["TEL_KB"],
               "telomere_short": cnt["TEL_SHORT"], "non_telomeric": len(chrom_ends) - cnt["TEL_KB"] - cnt["TEL_SHORT"],
               "contigs": len(crecords), "telomeric_contig_ends": len(tel),
               "contigs_telomeric_both_ends": sorted(s for s, n in both.items() if n == 2),
               "interstitial_arrays": sum(1 for r in rows if r["position_in_contig"] == "interstitial"),
               "telomeres_at_internal_joins_kb": sum(1 for r in rows if r["position_in_contig"] != "interstitial"
                                                     and r["faces"] == "join" and r["tier"] == "kb"),
               "unplaced_telomeric_contigs": sorted({d["contig"] for d in tel if not d["placed"]}), "self_checks": checks}
    write_tsv(st.path("contig_ends.tsv"), list(cend.values()))
    write_tsv(st.path("scaffold_ends.tsv"), send)
    write_tsv(st.path("chromosome_ends.tsv"), chrom_ends)
    write_tsv(st.path("contig_arrays.tsv"), rows)
    write_tsv(st.path("internal_candidates.tsv"), candidates, list(rows[0]) if rows else ["contig"])
    with open(st.path("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    st.commit({"step": "integrate", "label": a.label, "contigs_scan": os.path.abspath(a.contigs_scan),
               "scaffolds_scan": os.path.abspath(a.scaffolds_scan), "agp": os.path.abspath(a.agp), "agp_sha256": sha256_file(a.agp),
               "contigs_fasta_sha256": crec["fasta_sha256"], "scaffolds_fasta_sha256": srec["fasta_sha256"],
               "join_table": a.join_table, "split_ends": a.split_ends, "science": sci, "summary": summary})
    print(json.dumps(summary, indent=1))
    return 0


# ---------------------------------------------------------------- tips
def require_index(fa):
    """samtools faidx would silently write a missing index next to the input; refuse instead."""
    need = [fa + ".fai"] + ([fa + ".gzi"] if fa.endswith(".gz") else [])
    miss = [x for x in need if not os.path.exists(x)]
    if miss:
        sys.exit(f"missing index {miss}: create it yourself (samtools faidx); this tool writes nothing next to its inputs")


def faidx(fa, region, samtools="samtools"):
    out = subprocess.run([samtools, "faidx", fa, region], capture_output=True, text=True, check=True).stdout
    return "".join(out.splitlines()[1:])


def tip_call(alignments, ref_len, outer0, outer1, near, min_mapq):
    """alignments: PAF rows of one segment. Overlap-weighted window [outer0, outer1). -> dict"""
    A = [x for x in alignments if int(x[11]) >= min_mapq and int(x[2]) < outer1 and int(x[3]) > outer0]
    if not A:
        return {"call": "NO_ALIGNMENT"}
    w = lambda x: min(int(x[3]), outer1) - max(int(x[2]), outer0)
    by = collections.Counter()
    for x in A:
        by[x[5]] += w(x)
    chrom, bp = by.most_common(1)[0]
    bs = collections.Counter()
    for x in A:
        if x[5] == chrom:
            bs[x[4]] += w(x)
    strand = bs.most_common(1)[0][0]
    anc = min((x for x in A if x[5] == chrom and x[4] == strand), key=lambda x: int(x[2]))
    qs = int(anc[2])
    ref0 = int(anc[7]) - qs if strand == "+" else int(anc[8]) + qs
    CL = ref_len[chrom]
    call = "REF_START" if strand == "+" and ref0 <= near else ("REF_END" if strand == "-" and ref0 >= CL - near else "INTERIOR_or_INVERTED")
    return {"call": call, "ref_chr": chrom, "ref_strand": strand, "frac_on_chr": round(bp / sum(by.values()), 3),
            "ref_position_of_end": ref0, "ref_chr_length": CL, "dist_to_ref_end": min(abs(ref0), abs(CL - ref0))}


def cmd_tips(a):
    sci = load_science(a.science)
    t = sci["tips"]
    rec = require_validated(a.integrated)
    mm = find_tool(a.minimap2, "G2F_MINIMAP2", "minimap2")
    for fa in (a.contigs, a.scaffolds):
        require_index(fa)
    st = Staged(a.out)
    seg = st.path("segments.fa")
    jobs = []
    for r in read_tsv(os.path.join(a.integrated, "contig_ends.tsv")):
        if r["class"] != "ABSENT":
            jobs.append(("ctg", a.contigs, r["contig"], int(r["contig_length"]), r["side"], r["class"]))
    for r in read_tsv(os.path.join(a.integrated, "chromosome_ends.tsv")):
        if r["class"] not in ("TEL_KB", "TEL_SHORT"):
            jobs.append(("scf", a.scaffolds, r["scaffold"], int(r["scaffold_length"]), "L" if r["end"] == "start" else "R", r["class"]))
    with open(seg, "w") as fo:
        for kind, fa, seq, L, side, cls in jobs:
            n = min(t["segment_bp"], L)
            s = faidx(fa, f"{seq}:1-{n}") if side == "L" else revcomp(faidx(fa, f"{seq}:{L - n + 1}-{L}"))
            fo.write(f">{kind}|{seq}|{side}\n")
            for i in range(0, len(s), 100):
                fo.write(s[i:i + 100] + "\n")
    log(f"tips: {len(jobs)} segments -> minimap2 {t['preset']}")
    cmd = [mm, "-t", str(a.threads), "-x", t["preset"], "-c", "--secondary=no", a.reference, seg]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed: {p.stderr[-2000:]}")
    paf = [l.split("\t") for l in p.stdout.splitlines()]
    with open(st.path("segments.paf"), "w") as fh:
        fh.write("\n".join("\t".join(x[:12]) for x in paf) + ("\n" if paf else ""))
    os.remove(seg)                                   # regenerable from the inputs; keep only the PAF and calls
    ref_len = {x[5]: int(x[6]) for x in paf}          # target lengths from the PAF (no reference index needed)
    by = collections.defaultdict(list)
    for x in paf:
        by[x[0]].append(x)
    rows = []
    for kind, fa, seq, L, side, cls in jobs:
        q = f"{kind}|{seq}|{side}"
        o = tip_call(by[q], ref_len, 0, t["outer_bp"], t["near_end_bp"], t["min_mapq"])
        i = tip_call(by[q], ref_len, t["inner_bp"][0], t["inner_bp"][1], t["near_end_bp"], t["min_mapq"])
        row = {"segment": q, "kind": "contig_end" if kind == "ctg" else "chromosome_end_without_telomere", "sequence": seq,
               "side": side, "telomere_class": cls}
        row.update({f"outer_{k}": v for k, v in o.items()})
        row.update({f"inner_{k}": v for k, v in i.items()})
        rows.append(row)
    write_tsv(st.path("tip_calls.tsv"), rows)
    st.commit({"step": "tips", "integrated": os.path.abspath(a.integrated), "reference": os.path.abspath(a.reference),
               "reference_sha256": sha256_file(a.reference), "minimap2": mm, "command": " ".join(cmd), "science": sci,
               "segments": len(jobs), "label": rec.get("label")})
    return 0


# ---------------------------------------------------------------- closure
def cmd_closure(a):
    sci = load_science(a.science)
    c = sci["closure"]
    require_validated(a.integrated)
    sam = find_tool(None, "G2F_SAMTOOLS", "samtools")
    ref = ["-T", a.reference_fasta] if a.reference_fasta else []
    hdr = subprocess.run([sam, "view", "-H", *ref, a.alignments], capture_output=True, text=True, check=True).stdout
    sq = {}
    for l in hdr.splitlines():
        if l.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in l.split("\t")[1:] if ":" in x)
            sq[f["SN"]] = f
    # the alignment reference must be the contig set: same names, lengths and, where present, MD5
    md5 = {}
    with opener(a.contigs) as fh:
        name, h = None, None
        for line in fh:
            if line.startswith(b">"):
                if name:
                    md5[name] = h.hexdigest()
                name, h = line[1:].split()[0].decode(), hashlib.md5()
            else:
                h.update(line.strip().upper())
        if name:
            md5[name] = h.hexdigest()
    st = Staged(a.out)
    rows = []
    ends = [r for r in read_tsv(os.path.join(a.integrated, "contig_ends.tsv")) if r["class"] in ("TEL_KB", "TEL_SHORT")]
    for r in ends:
        ctg, side, L = r["contig"], r["side"], int(r["contig_length"])
        s = sq.get(ctg)
        if not s or int(s["LN"]) != L or ("M5" in s and s["M5"] != md5.get(ctg)):
            rows.append(dict(contig=ctg, side=side, status="REFERENCE_MISMATCH_OR_ABSENT"))
            continue
        lo, hi = (1, c["outer_bp"]) if side == "L" else (L - c["outer_bp"] + 1, L)
        out = subprocess.run([sam, "view", *ref, "-F", "0x904", "-q", str(c["min_mapq"]), a.alignments, f"{ctg}:{lo}-{hi}"],
                             capture_output=True, text=True, check=True).stdout
        reach = clipped = 0
        fracs = []
        for l in out.splitlines():
            x = l.split("\t")
            pos, cig, seq = int(x[3]), x[5], x[9]
            span = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cig))
            lc = re.match(r"^(\d+)S", cig)
            rc_ = re.search(r"(\d+)S$", cig)
            if side == "L" and pos <= c["outer_bp"]:
                reach += 1
                n = int(lc.group(1)) if lc else 0
                part = seq[:n]
            elif side == "R" and pos + span - 1 >= L - c["outer_bp"] + 1:
                reach += 1
                n = int(rc_.group(1)) if rc_ else 0
                part = seq[len(seq) - n:] if n else ""
            else:
                continue
            if n >= c["min_clip_bp"]:
                clipped += 1
                m = sci["motif"].upper()
                fracs.append(round(len(m) * (part.upper().count(m) + part.upper().count(revcomp(m))) / n, 2))
        status = ("NO_READS" if reach == 0 else "CLOSED" if clipped == 0 else
                  "CLIPPED_TELOMERIC" if all(f >= 0.5 for f in fracs) else "CLIPPED_NON_TELOMERIC")
        rows.append(dict(contig=ctg, side=side, telomere_class=r["class"], placed=r["placed"], reads_reaching_end=reach,
                         outward_clipped=clipped, clipped_motif_fraction=",".join(map(str, fracs)), status=status))
    write_tsv(st.path("closure.tsv"), rows)
    st.commit({"step": "closure", "integrated": os.path.abspath(a.integrated), "alignments": os.path.abspath(a.alignments),
               "contigs": os.path.abspath(a.contigs), "science": sci,
               "status_counts": dict(collections.Counter(r["status"] for r in rows)),
               "note": "CLOSED / CLIPPED_TELOMERIC: reads stop at the telomere (a terminus). CLIPPED_NON_TELOMERIC: sequence continues "
                       "past the array in some reads (inspect: interstitial array or a different haplotype). Read counts at "
                       "telomeric ends are usually small; treat as supporting evidence."})
    return 0


# ---------------------------------------------------------------- compare
def chrom_end_classes(d, label, sci):
    rec, ends, arrays, records = scan_tables(d)
    cr = re.compile(sci["chromosome_regex"])
    parts = collections.defaultdict(list)
    for o in records:
        m = cr.search(o)
        if m:
            parts[f"chr{int(m.group(1)):02d}"].append((int(m.group(2) or 0), o))
    out = {}
    for ch, v in parts.items():
        ol = [o for _, o in sorted(v)]
        for end, o, side in (("start", ol[0], "L"), ("end", ol[-1], "R")):
            out[(ch, end)] = classify_end(o, int(records[o]["length"]), side, ends, arrays, sci)
    internal = []
    c = sci["classes"]
    keep = {o for v in parts.values() for _, o in v}
    for o in keep:
        for r in arrays.get(o, []):
            if r["motif_bp"] >= c["kb_bp"] and r["dist_start"] > c["terminal_window"] and r["dist_end"] > c["terminal_window"]:
                internal.append(r)
    return out, internal


def cmd_compare(a):
    sci = load_science(a.science)
    b, ib = chrom_end_classes(a.before, a.label, sci)
    f, ia = chrom_end_classes(a.after, a.label, sci)
    st = Staged(a.out)
    tel = lambda c: c in ("TEL_KB", "TEL_SHORT")
    rows = []
    for k in sorted(set(b) | set(f)):
        x, y = b.get(k, {}).get("class", "MISSING"), f.get(k, {}).get("class", "MISSING")
        rows.append(dict(chromosome=k[0], end=k[1], before=x, after=y,
                         change="improved" if tel(y) and not tel(x) else "REGRESSED" if tel(x) and not tel(y) else "unchanged" if x == y else "class changed"))
    key = lambda r: (r["seq"], r["start0"], r["end"])
    kb = {key(r) for r in ib}
    arr = [dict(version="before", **{k: r[k] for k in ("seq", "start0", "end", "motif_bp", "strand")}) for r in ib] + \
          [dict(version="after", new=key(r) not in kb, **{k: r[k] for k in ("seq", "start0", "end", "motif_bp", "strand")}) for r in ia]
    write_tsv(st.path("chromosome_ends.tsv"), rows)
    write_tsv(st.path("internal_arrays.tsv"), arr, ["version", "seq", "start0", "end", "motif_bp", "strand", "new"])
    s = {"telomeric_ends_before": sum(tel(r["before"]) for r in rows), "telomeric_ends_after": sum(tel(r["after"]) for r in rows),
         "improved": sum(r["change"] == "improved" for r in rows), "regressed": sum(r["change"] == "REGRESSED" for r in rows),
         "internal_kb_arrays_before": len(ib), "internal_kb_arrays_after": len(ia),
         "new_internal_arrays_after": sum(1 for r in ia if key(r) not in kb)}
    with open(st.path("summary.json"), "w") as fh:
        json.dump(s, fh, indent=1)
    st.commit({"step": "compare", "label": a.label, "before": os.path.abspath(a.before), "after": os.path.abspath(a.after),
               "science": sci, "summary": s})
    print(json.dumps(s, indent=1))
    return 0 if s["regressed"] == 0 and s["new_internal_arrays_after"] == 0 else 3


# ---------------------------------------------------------------- plot
def cmd_plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle
    KB, SHORT, NONE_F, NONE_E, CRIT, INK, INK2, GRID, SURF = ("#1c5cab", "#86b6ef", "#f0efec", "#8a8984", "#d03b3b",
                                                              "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb")
    plt.rcParams.update({"font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42, "figure.facecolor": SURF,
                         "axes.facecolor": SURF, "savefig.facecolor": SURF})
    st = Staged(a.out)
    color = {"TEL_KB": KB, "TEL_SHORT": SHORT}

    def wedge(ax, x, y, h, w, facing, fc, ec=None):
        pts = [(x, y), (x + w, y + h / 2), (x + w, y - h / 2)] if facing == "L" else [(x, y), (x - w, y + h / 2), (x - w, y - h / 2)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=fc, edgecolor=ec or fc, lw=0.8, zorder=4))

    labels, mats = [], {}
    for d in a.integrated:
        rec = require_validated(d)
        lab = rec["label"]
        labels.append(lab)
        ends = read_tsv(os.path.join(d, "chromosome_ends.tsv"))
        cand = read_tsv(os.path.join(d, "internal_candidates.tsv"))
        mats[lab] = {(r["chromosome"], r["end"]): r["class"] for r in ends}
        chroms = sorted({r["chromosome"] for r in ends})
        L = {}
        for r in read_tsv(os.path.join(d, "scaffold_ends.tsv")):
            L.setdefault(r["chromosome"], {})[r["scaffold"]] = int(r["scaffold_length"])
        fig, ax = plt.subplots(figsize=(10, 0.45 * len(chroms) + 1.2))
        for i, ch in enumerate(chroms):
            y = len(chroms) - i
            off = 0
            objs = sorted(L[ch])
            for o in objs:
                ax.add_patch(Rectangle((off / 1e6, y - 0.16), L[ch][o] / 1e6, 0.32, facecolor="#e6e5e1", edgecolor="none"))
                off += L[ch][o] + (2_000_000 if o != objs[-1] else 0)
            for end, x, facing in (("start", 0, "L"), ("end", off / 1e6, "R")):
                c = mats[lab].get((ch, end))
                wedge(ax, x, y, 0.42, 2.2, facing, color.get(c, NONE_F), None if c in color else NONE_E)
            for r in cand:
                if r["scaffold"].startswith(tuple(objs)) and r["scaffold"] in objs and r["tier"] == "kb":
                    base = sum(L[ch][o] + 2_000_000 for o in objs[:objs.index(r["scaffold"])])
                    facing = "L" if r["scaffold_strand"] == "C-rich" else "R"
                    x = (base + (int(r["scaffold_start"]) if facing == "L" else int(r["scaffold_end"]))) / 1e6
                    wedge(ax, x, y + 0.36, 0.34, 1.8, facing, CRIT)
        ax.set_yticks(range(len(chroms), 0, -1))
        ax.set_yticklabels(chroms)
        ax.set_xlabel("scaffold position (Mb)")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.set_title(f"{lab}: telomeres at chromosome ends (wedge points the way the telomere faces; red = not at a chromosome end)",
                     loc="left", fontsize=8.5, color=INK)
        for ext in ("svg", "pdf"):
            fig.savefig(st.path(f"map_{lab}.{ext}"), bbox_inches="tight")
        plt.close(fig)
    cols = sorted({k for m in mats.values() for k in m}, key=lambda k: (k[0], k[1] != "start"))
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(cols) + 3), 0.32 * len(labels) + 1.5))
    for i, lab in enumerate(labels):
        y = len(labels) - i - 1
        for j, k in enumerate(cols):
            c = mats[lab].get(k, "ABSENT")
            fc = color.get(c, NONE_F)
            ax.add_patch(Rectangle((j + 0.05, y + 0.05), 0.9, 0.9, facecolor=fc, edgecolor=NONE_E if fc == NONE_F else SURF))
        ax.text(-0.3, y + 0.5, lab, ha="right", va="center", color=INK)
        ax.text(len(cols) + 0.2, y + 0.5, f"{sum(mats[lab].get(k) in color for k in cols)}/{len(cols)}", va="center", color=INK2)
    for j, k in enumerate(cols):
        ax.text(j + 0.5, len(labels) + 0.1, f"{k[0]}{'S' if k[1] == 'start' else 'E'}", rotation=90, ha="center", va="bottom", fontsize=6)
    ax.set_xlim(-6, len(cols) + 2)
    ax.set_ylim(-0.5, len(labels) + 2)
    ax.axis("off")
    h = [Rectangle((0, 0), 1, 1, facecolor=KB), Rectangle((0, 0), 1, 1, facecolor=SHORT), Rectangle((0, 0), 1, 1, facecolor=NONE_F, edgecolor=NONE_E)]
    ax.legend(h, ["telomere >= 1 kb", "200-999 bp", "none"], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.05))
    for ext in ("svg", "pdf"):
        fig.savefig(st.path(f"end_heatmap.{ext}"), bbox_inches="tight")
    plt.close(fig)
    st.commit({"step": "plot", "integrated": [os.path.abspath(d) for d in a.integrated]})
    return 0


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check-env")
    p.add_argument("--tidk")
    p.set_defaults(func=cmd_check_env)
    p = sub.add_parser("scan")
    p.add_argument("--fasta", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--science")
    p.add_argument("--tidk")
    p.add_argument("--expected-sha256")
    p.add_argument("--explore", action="store_true", help="also run a limited de novo tidk explore (-m 5 -x 12)")
    p.set_defaults(func=cmd_scan)
    p = sub.add_parser("integrate")
    for k in ("--contigs-scan", "--scaffolds-scan", "--agp", "--label", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--join-table", help="TSV with contig_A, contig_B[, class, evidence_call]")
    p.add_argument("--split-ends", help="TSV with contig, side (L/R)[, note]: contig ends made by curation splits")
    p.add_argument("--science")
    p.set_defaults(func=cmd_integrate)
    p = sub.add_parser("tips")
    for k in ("--integrated", "--contigs", "--scaffolds", "--reference", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--minimap2")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--science")
    p.set_defaults(func=cmd_tips)
    p = sub.add_parser("closure")
    for k in ("--integrated", "--alignments", "--contigs", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--reference-fasta", help="for CRAM decoding, if the header path is not usable")
    p.add_argument("--science")
    p.set_defaults(func=cmd_closure)
    p = sub.add_parser("compare")
    for k in ("--before", "--after", "--label", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--science")
    p.set_defaults(func=cmd_compare)
    p = sub.add_parser("plot")
    p.add_argument("--integrated", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_plot)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
