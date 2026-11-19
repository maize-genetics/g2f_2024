#!/usr/bin/env python3
"""Foreign-contamination and organelle screens of an assembly, and a candidate exclusion plan for review.

Subcommands (each writes to a new directory, built as <dir>.partial and renamed after its outputs validate; never
overwrites; inputs are only read):
  check-env      report minimap2, FCS and container-engine availability
  fcs            run NCBI FCS-adaptor and FCS-GX on one FASTA, then summarise the calls (as fcs-summarize)
  fcs-summarize  summarise existing FCS-adaptor / FCS-GX output directories (e.g. a screen run elsewhere)
  organelle      window-based organelle screen: every record (or every unplaced record, with --agp) is cut into
                 windows, each window aligned to organelle references with minimap2; organelle-dominant records get a
                 nuclear test against a nuclear reference; classes CONFIDENT_ORGANELLE / ORGANELLE_DOMINANT_MIXED
  plan           combine fcs and organelle results of one or more assemblies into a CANDIDATE exclusion plan for
                 scripts/g2f_release_exclude.py, plus a review table of everything that is not a whole-record call

  g2f_decontam.py check-env     [--fcs-dir DIR]
  g2f_decontam.py fcs           --fasta X.fa[.gz] --label UNIT --tax-id N --gx-db DIR/all --out DIR
                                [--agp X.agp] [--fcs-dir DIR] [--container-engine singularity] [--threads 48]
                                [--skip-gx] [--science YAML]
  g2f_decontam.py fcs-summarize --fasta X.fa[.gz] --label UNIT --adaptor-dir D --gx-dir D --out DIR [--agp X.agp]
  g2f_decontam.py organelle     --fasta X.fa[.gz] --label UNIT --organelles ORG.fa --out DIR [--agp X.agp]
                                [--all-records] [--plastid-ids ID[,ID...]] [--nuclear NUC.fa] [--threads 16]
  g2f_decontam.py plan          [--fcs DIR ...] [--organelle DIR ...] --out DIR

Scientific parameters: config/decontam_science.yaml (built-in defaults below; tests keep the two identical).

What it decides, and what it does not:
  - Calls come from the tools' own output (FCS-GX actions, FCS-adaptor actions, alignment coverage). Nothing is
    inferred from GC content, depth, repeats or placement.
  - Only whole-record calls of unplaced singleton records become plan rows: FCS-GX EXCLUDE (FOREIGN_CONTAMINANT) and
    CONFIDENT_ORGANELLE (ORGANELLE_PLASTID / ORGANELLE_MITOCHONDRION / ORGANELLE). Partial calls (GX TRIM / FIX, adaptor
    trims), ambiguous GX calls, mixed organelle records, conflicts and calls on placed contigs go to REVIEW.tsv.
  - The plan is a candidate. Excluding a record is a human decision; review the plan and the review table, then apply
    it with g2f_release_exclude.py and its independent verifier.
"""
import argparse
import collections
import csv
import datetime
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

SCHEMA = "g2f-decontam/1"
with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()

# Defaults = config/decontam_science.yaml (tests/test_decontam.py keeps them equal).
DEFAULTS = {
    "version": 1,
    "fcs": {
        "adaptor_mode": "euk",
        "contaminant_actions": ["EXCLUDE", "TRIM", "FIX"],
        "review_actions": ["REVIEW", "REVIEW_RARE", "INFO"],
        "gx_organelle_dominant_frac": 0.5,
    },
    "organelle": {
        "window_bp": 5000,
        "minimap2_preset": "asm5",
        "minimap2_extra": ["-N", "50", "-p", "0.1"],
        "min_identity": 0.90,
        "confident_frac": 0.90,
        "dominant_frac": 0.50,
        "max_nuclear_bp": 5000,
        "min_nonorganelle_segment_bp": 200,
        "nuclear_preset": "asm5",
    },
    "chromosome_regex": r"_chr(\d+)(?:_part(\d+))?$",
}
PLAN_COLS = ["unit", "record", "reason_class", "evidence"]


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


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def write_tsv(path, rows, keys, comment=None):
    with open(path, "w") as fh:
        if comment:
            for c in comment:
                fh.write(f"# {c}\n")
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


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


def find_tool(explicit, env, name):
    return explicit or os.environ.get(env) or shutil.which(name) or name


def merge_iv(iv):
    """Union of half-open [b, e) intervals."""
    out = []
    for b, e in sorted(iv):
        if out and b <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([b, e])
    return out


def iv_bp(iv):
    return sum(e - b for b, e in merge_iv(iv))


def fasta_records(path, keep=None):
    """Yield (id, sequence) from a FASTA (plain or gzip); with keep, only those IDs are assembled."""
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif keep is None or name in keep:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        yield name, "".join(buf)


def fasta_lengths(path):
    """Record lengths; from a samtools .fai next to the FASTA when one exists and is newer, else by reading it."""
    fai = path + ".fai"
    if os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(path):
        with open(fai) as fh:
            lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fh if l.strip()}
        if lens:
            return lens
    op = gzip.open if path.endswith(".gz") else open
    lens, name = {}, None
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in lens:
                    sys.exit(f"{path}: duplicate record {name}")
                lens[name] = 0
            else:
                lens[name] += len(line.strip())
    if not lens:
        sys.exit(f"{path}: no FASTA records")
    return lens


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

    def commit(self, record):
        files = {}
        for d, _, fs in os.walk(self.partial):
            for f in sorted(fs):
                p = os.path.join(d, f)
                if not os.path.islink(p):
                    files[os.path.relpath(p, self.partial)] = sha256_file(p)
        record.update(schema=SCHEMA, script=os.path.abspath(__file__), script_sha256=SCRIPT_SHA256,
                      python=sys.version.split()[0], validated=now(), outputs_sha256=files)
        with open(os.path.join(self.partial, "VALIDATED.json"), "w") as fh:
            json.dump(record, fh, indent=1)
        os.rename(self.partial, self.out)
        log(f"validated: {self.out}")


def require_validated(d, step):
    v = os.path.join(d, "VALIDATED.json")
    if not os.path.exists(v):
        sys.exit(f"{d}: no VALIDATED.json; run the producing step first")
    with open(v) as fh:
        j = json.load(fh)
    if j.get("schema") != SCHEMA or j.get("step") != step:
        sys.exit(f"{d}: not a g2f_decontam {step} directory")
    return j


# ---------------------------------------------------------------- placement
def placement(agp, lens, chrom_regex):
    """record -> (status, object). status: unplaced (the sole component of a non-chromosome object of the same name,
    used whole in + orientation: what g2f_release_exclude.py can remove), unplaced_nonstandard (sole component of a
    non-chromosome object, but renamed, partial or not +), placed, or absent (not in the AGP)."""
    objs = collections.defaultdict(list)
    with open(agp) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            objs[f[0]].append(f)
    rx = re.compile(chrom_regex)
    out = {}
    for obj, fs in objs.items():
        for f in fs:
            if f[4] in ("N", "U"):
                continue
            rec = f[5]
            if rec not in lens:
                sys.exit(f"{agp}: component {rec} is not in the FASTA (wrong AGP?)")
            if len(fs) == 1 and not rx.search(obj):
                std = obj == rec and int(f[6]) == 1 and int(f[7]) == lens[rec] and f[8] == "+"
                st = "unplaced" if std else "unplaced_nonstandard"
            else:
                st = "placed"
            if rec in out:
                sys.exit(f"{agp}: component {rec} is used more than once")
            out[rec] = (st, obj)
    for rec in lens:
        out.setdefault(rec, ("absent", ""))
    return out


# ---------------------------------------------------------------- FCS parsing
def header_json(path):
    with open(path) as fh:
        first = fh.readline()
    if first.startswith("##"):
        try:
            return json.loads(first[2:])[1]
        except (ValueError, IndexError, KeyError):
            return {}
    return {}


def one(pattern, what):
    hits = sorted(glob.glob(pattern))
    if len(hits) != 1:
        sys.exit(f"expected exactly one {what} ({pattern}), found {len(hits)}")
    return hits[0]


def parse_adaptor(report):
    calls = []
    with open(report) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            acc, ln, act, rng, name = (line.rstrip("\n").split("\t") + [""] * 5)[:5]
            if act not in ("ACTION_EXCLUDE", "ACTION_TRIM"):
                continue
            ivs = [(1, int(ln))] if act == "ACTION_EXCLUDE" else \
                [tuple(map(int, x.split(".."))) for x in rng.split(",") if ".." in x]
            for s, e in ivs:
                calls.append(dict(seq_id=acc, seq_len=int(ln), start=s, end=e, tool="FCS-adaptor", action=act,
                                  call_class="adaptor_vector", whole_record=act == "ACTION_EXCLUDE", div="",
                                  agg_cont_cov="", top_tax_name=name))
    return calls


def parse_gx(report, sci):
    con, rev = set(sci["fcs"]["contaminant_actions"]), set(sci["fcs"]["review_actions"])
    calls = []
    with open(report) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            sid, s, e, ln, act, div, cov, tax = (line.rstrip("\n").split("\t") + [""] * 8)[:8]
            s, e, ln = int(s), int(e), int(ln)
            if act == "EXCLUDE":
                s, e = 1, ln
            cls = "foreign_contaminant" if act in con else "gx_review" if act in rev else "gx_other"
            calls.append(dict(seq_id=sid, seq_len=ln, start=s, end=e, tool="FCS-GX", action=act, call_class=cls,
                              whole_record=act == "EXCLUDE", div=div, agg_cont_cov=cov, top_tax_name=tax))
    return calls


def parse_gx_taxonomy(rpt, sci):
    """Organelle-dominant sequences from the GX taxonomy report ((xp,lc,co,n,mt,pt,pm)-len column)."""
    rows = []
    frac = sci["fcs"]["gx_organelle_dominant_frac"]
    with open(rpt) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                ln = int(f[1])
                parts = [int(x) for x in f[2].split(",")]
            except (ValueError, IndexError):
                continue
            mt, pt = parts[4], parts[5]
            if ln and (mt + pt) / ln >= frac:
                rows.append(dict(seq_id=f[0], seq_len=ln, mt_bp=mt, pt_bp=pt, organelle_frac=f"{(mt + pt) / ln:.4f}",
                                 top_tax_name=f[5] if len(f) > 5 else ""))
    return rows


FCS_CALL_COLS = ["seq_id", "seq_len", "placement", "object", "start", "end", "tool", "action", "call_class",
                 "whole_record", "div", "agg_cont_cov", "top_tax_name"]


def summarize_fcs(st, a, sci, adaptor_dir, gx_dir, extra):
    lens = fasta_lengths(a.fasta)
    place = placement(a.agp, lens, sci["chromosome_regex"]) if a.agp else {r: ("no_agp", "") for r in lens}
    rep = os.path.join(adaptor_dir, "fcs_adaptor_report.txt")
    if not os.path.isfile(rep):
        sys.exit(f"missing {rep}")
    calls = parse_adaptor(rep)
    gx_meta, org = {}, []
    if gx_dir:
        grep = one(os.path.join(gx_dir, "*.fcs_gx_report.txt"), "FCS-GX report")
        trep = one(os.path.join(gx_dir, "*.taxonomy.rpt"), "GX taxonomy report")
        gx_meta = header_json(grep)
        calls += parse_gx(grep, sci)
        org = parse_gx_taxonomy(trep, sci)
    for c in calls + org:
        if c["seq_id"] not in lens:
            sys.exit(f"FCS output names {c['seq_id']}, which is not in {a.fasta} (wrong FASTA?)")
        if int(c["seq_len"]) != lens[c["seq_id"]]:
            sys.exit(f"FCS length of {c['seq_id']} ({c['seq_len']}) differs from the FASTA ({lens[c['seq_id']]})")
        c["placement"], c["object"] = place[c["seq_id"]]
    total = sum(lens.values())
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for c in calls:
        by[c["call_class"]][c["seq_id"]].append((c["start"] - 1, c["end"]))
    union = collections.defaultdict(list)
    for k in ("foreign_contaminant", "adaptor_vector"):
        for s, iv in by[k].items():
            union[s] += iv

    def bp(d):
        return sum(iv_bp(v) for v in d.values())
    divs = collections.Counter()
    for c in calls:
        if c["call_class"] == "foreign_contaminant":
            divs[c["div"]] += c["end"] - c["start"] + 1
    ri = gx_meta.get("run-info", {})
    summary = dict(unit=a.label, records=len(lens), total_bp=total, gx_run=bool(gx_dir),
                   foreign_bp=bp(by["foreign_contaminant"]), foreign_records=len(by["foreign_contaminant"]),
                   foreign_by_div=dict(sorted(divs.items())),
                   gx_review_bp=bp(by["gx_review"]), gx_review_records=len(by["gx_review"]),
                   adaptor_bp=bp(by["adaptor_vector"]), adaptor_records=len(by["adaptor_vector"]),
                   union_flagged_bp=bp(union),
                   union_flagged_pct=round(100 * bp(union) / total, 6),
                   gx_organelle_dominant_records=len(org),
                   gx_organelle_dominant_bp=sum(r["seq_len"] for r in org),
                   asserted_div=ri.get("asserted-div", ""), inferred_primary_divs=ri.get("inferred-primary-divs", []),
                   gx_version=gx_meta.get("git-rev", ""), gx_db_build=gx_meta.get("db", {}).get("build-date", ""))
    write_tsv(st.path("fcs_calls.tsv"), calls, FCS_CALL_COLS)
    write_tsv(st.path("gx_organelle_dominant.tsv"), org,
              ["seq_id", "seq_len", "placement", "object", "mt_bp", "pt_bp", "organelle_frac", "top_tax_name"])
    with open(st.path("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    rec = dict(step="fcs", label=a.label, fasta=os.path.abspath(a.fasta), fasta_sha256=sha256_file(a.fasta),
               agp=os.path.abspath(a.agp) if a.agp else None, agp_sha256=sha256_file(a.agp) if a.agp else None,
               science=sci, summary=summary)
    rec.update(extra)
    st.commit(rec)
    log(f"{a.label}: {summary['foreign_records']} foreign, {summary['adaptor_records']} adaptor, "
        f"{summary['gx_review_records']} GX review records")


# ---------------------------------------------------------------- subcommands
def tool_version(cmd):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
        return (p.stdout or p.stderr).strip().splitlines()[0] if (p.stdout or p.stderr).strip() else None
    except (FileNotFoundError, PermissionError):
        return None


def cmd_check_env(a):
    ok = True
    mm = find_tool(None, "G2F_MINIMAP2", "minimap2")
    v = tool_version([mm, "--version"])
    print(f"minimap2: {mm} {v or 'NOT FOUND'}")
    ok &= v is not None
    fcs = a.fcs_dir or os.environ.get("G2F_FCS_DIR")
    if fcs:
        for f in ("fcs.py", "run_fcsadaptor.sh", "fcs-gx.sif", "fcs-adaptor.sif"):
            p = os.path.join(fcs, f)
            print(f"FCS {f}: {'ok' if os.path.exists(p) else 'MISSING'} ({p})")
            ok &= os.path.exists(p)
    else:
        print("FCS: not checked (give --fcs-dir or G2F_FCS_DIR; needed only for `fcs`)")
    for e in ("singularity", "apptainer", "docker"):
        print(f"{e}: {shutil.which(e) or 'not found'}")
    sys.exit(0 if ok else 1)


def cmd_fcs(a):
    sci = load_science(a.science)
    fcs = a.fcs_dir or os.environ.get("G2F_FCS_DIR")
    if not fcs:
        sys.exit("give --fcs-dir (the directory with fcs.py, run_fcsadaptor.sh and the .sif images) or G2F_FCS_DIR")
    if not a.skip_gx and not glob.glob(a.gx_db + "*.gxi"):
        sys.exit(f"--gx-db {a.gx_db}: no {a.gx_db}*.gxi (give the database prefix, e.g. /path/gxdb/all)")
    if sci["fcs"]["adaptor_mode"] not in ("euk", "prok"):
        sys.exit("fcs.adaptor_mode must be euk or prok")
    fasta_lengths(a.fasta)                                   # fail early on an unreadable or empty FASTA
    st = Staged(a.out)
    env = dict(os.environ, NCBI_FCS_REPORT_ANALYTICS="0", DO_NOT_TRACK="1", GX_NUM_CORES=str(a.threads),
               FCS_DEFAULT_IMAGE=a.gx_image or os.path.join(fcs, "fcs-gx.sif"))
    cmds = []
    ad = os.path.join(st.partial, "adaptor")
    os.makedirs(ad)
    c = ["bash", os.path.join(fcs, "run_fcsadaptor.sh"), "--fasta-input", os.path.abspath(a.fasta),
         "--output-dir", ad, f"--{sci['fcs']['adaptor_mode']}", "--container-engine",
         a.container_engine, "--image", a.adaptor_image or os.path.join(fcs, "fcs-adaptor.sif")]
    log(f"{a.label}: FCS-adaptor")
    with open(st.path("adaptor.log"), "w") as lf:
        if subprocess.run(c, stdout=lf, stderr=subprocess.STDOUT, env=env).returncode:
            sys.exit(f"FCS-adaptor failed; see {st.path('adaptor.log')}")
    cmds.append(" ".join(c))
    gxd = None
    if not a.skip_gx:
        gxd = os.path.join(st.partial, "gx")
        os.makedirs(gxd)
        c = [sys.executable if a.fcs_python is None else a.fcs_python, os.path.join(fcs, "fcs.py"),
             "--no-report-analytics", "screen", "genome", "--fasta", os.path.abspath(a.fasta), "--out-dir", gxd,
             "--gx-db", a.gx_db, "--tax-id", str(a.tax_id)]
        log(f"{a.label}: FCS-GX (tax-id {a.tax_id}, {a.threads} cores)")
        with open(st.path("gx.log"), "w") as lf:
            if subprocess.run(c, stdout=lf, stderr=subprocess.STDOUT, env=env).returncode:
                sys.exit(f"FCS-GX failed; see {st.path('gx.log')}")
        cmds.append(" ".join(c))
    summarize_fcs(st, a, sci, ad, gxd,
                  dict(tax_id=a.tax_id, gx_db=a.gx_db, fcs_dir=os.path.abspath(fcs), commands=cmds,
                       gx_skipped=a.skip_gx, container_engine=a.container_engine))


def cmd_fcs_summarize(a):
    sci = load_science(a.science)
    st = Staged(a.out)
    summarize_fcs(st, a, sci, a.adaptor_dir, a.gx_dir,
                  dict(imported_from={"adaptor": os.path.abspath(a.adaptor_dir),
                                      "gx": os.path.abspath(a.gx_dir) if a.gx_dir else None},
                       gx_skipped=a.gx_dir is None))


def read_paf_blocks(paf_text, min_id):
    """query -> list of (qstart, qend, target) for blocks with identity >= min_id."""
    d = collections.defaultdict(list)
    for line in paf_text.splitlines():
        x = line.split("\t")
        if len(x) < 11 or int(x[10]) == 0 or int(x[9]) / int(x[10]) < min_id:
            continue
        d[x[0]].append((int(x[2]), int(x[3]), x[5]))
    return d


def organelle_class(L, org_bp, nuc_bp, nonorg_bp, o):
    """CONFIDENT_ORGANELLE / ORGANELLE_DOMINANT_MIXED / ''. Without a nuclear test, the whole non-organelle part
    counts as possible nuclear sequence (a conservative upper bound)."""
    f = org_bp / L if L else 0
    nuclear = nonorg_bp if nuc_bp is None else nuc_bp
    if f >= o["confident_frac"] and nuclear < o["max_nuclear_bp"]:
        return "CONFIDENT_ORGANELLE"
    if f >= o["dominant_frac"]:
        return "ORGANELLE_DOMINANT_MIXED"
    return ""


def cmd_organelle(a):
    sci = load_science(a.science)
    o = sci["organelle"]
    mm = find_tool(a.minimap2, "G2F_MINIMAP2", "minimap2")
    if tool_version([mm, "--version"]) is None:
        sys.exit(f"minimap2 not found ({mm}); set --minimap2 or G2F_MINIMAP2")
    org_ids = [r for r, _ in fasta_records(a.organelles)]
    plastid = set(x for x in (a.plastid_ids or "").split(",") if x)
    if plastid - set(org_ids):
        sys.exit(f"--plastid-ids not in {a.organelles}: {sorted(plastid - set(org_ids))}")
    lens = fasta_lengths(a.fasta)
    place = placement(a.agp, lens, sci["chromosome_regex"]) if a.agp else {r: ("no_agp", "") for r in lens}
    screen = [r for r in lens if a.all_records or not a.agp or place[r][0].startswith("unplaced")]
    if not screen:
        sys.exit("no records to screen")
    st = Staged(a.out)
    W = int(o["window_bp"])
    win = st.path("windows.fa")
    wmap, n = {}, 0
    keep = set(screen)
    with open(win, "w") as fo:
        for rec, s in fasta_records(a.fasta, keep):
            for i in range(0, len(s), W):
                wmap[f"w{n}"] = (rec, i)
                fo.write(f">w{n}\n{s[i:i + W]}\n")
                n += 1
    cmd1 = [mm, "-x", o["minimap2_preset"], "-c"] + [str(x) for x in o["minimap2_extra"]] + \
        ["-t", str(a.threads), a.organelles, win]
    log(f"{a.label}: {len(screen)} records, {n} windows vs {len(org_ids)} organelle references")
    paf1 = subprocess.run(cmd1, check=True, capture_output=True, text=True).stdout
    with open(st.path("windows_vs_organelles.paf"), "w") as fh:
        fh.write(paf1)
    os.remove(win)                                   # our own intermediate inside the .partial directory
    cov = collections.defaultdict(lambda: {"all": [], "pt": [], "mt": []})
    for q, blocks in read_paf_blocks(paf1, o["min_identity"]).items():
        rec, off = wmap[q]
        for b, e, t in blocks:
            iv = (off + b, off + e)
            cov[rec]["all"].append(iv)
            cov[rec]["pt" if t in plastid else "mt"].append(iv)
    # nuclear test of organelle-dominant records
    flagged = {r for r in screen if iv_bp(cov[r]["all"]) / lens[r] >= o["dominant_frac"]}
    nuc = collections.defaultdict(list)
    cmd2 = None
    if a.nuclear and flagged:
        seg = st.path("nonorganelle.fa")
        smap, k = {}, 0
        with open(seg, "w") as fo:
            for rec, s in fasta_records(a.fasta, flagged):
                pos = 0
                for b, e in merge_iv(cov[rec]["all"]) + [[len(s), len(s)]]:
                    if b - pos >= o["min_nonorganelle_segment_bp"]:
                        smap[f"s{k}"] = (rec, pos)
                        fo.write(f">s{k}\n{s[pos:b]}\n")
                        k += 1
                    pos = max(pos, e)
        if k:
            cmd2 = [mm, "-x", o["nuclear_preset"], "-c", "-t", str(a.threads), a.nuclear, seg]
            paf2 = subprocess.run(cmd2, check=True, capture_output=True, text=True).stdout
            with open(st.path("nonorganelle_vs_nuclear.paf"), "w") as fh:
                fh.write(paf2)
            for q, blocks in read_paf_blocks(paf2, o["min_identity"]).items():
                rec, off = smap[q]
                nuc[rec] += [(off + b, off + e) for b, e, _ in blocks]
        os.remove(seg)
    rows = []
    for rec in screen:
        L = lens[rec]
        ob, pt, mt = iv_bp(cov[rec]["all"]), iv_bp(cov[rec]["pt"]), iv_bp(cov[rec]["mt"])
        nb = iv_bp(nuc[rec]) if (a.nuclear and rec in flagged) else None
        cls = organelle_class(L, ob, nb, L - ob, o)
        kind = "" if not cls else ("ORGANELLE" if not plastid else
                                   "ORGANELLE_PLASTID" if pt >= mt else "ORGANELLE_MITOCHONDRION")
        rows.append(dict(record=rec, placement=place[rec][0], object=place[rec][1], length=L, organelle_bp=ob,
                         organelle_frac=f"{ob / L:.4f}", plastid_bp=pt if plastid else "",
                         mitochondrion_bp=mt if plastid else "", nonorganelle_bp=L - ob,
                         nonorganelle_nuclear_bp="" if nb is None else nb, window_class=cls, organelle_kind=kind))
    write_tsv(st.path("organelle_screen.tsv"), rows, list(rows[0]))
    summary = dict(unit=a.label, records_screened=len(rows), bp_screened=sum(r["length"] for r in rows), windows=n,
                   confident=sum(r["window_class"] == "CONFIDENT_ORGANELLE" for r in rows),
                   confident_bp=sum(r["length"] for r in rows if r["window_class"] == "CONFIDENT_ORGANELLE"),
                   dominant_mixed=sum(r["window_class"] == "ORGANELLE_DOMINANT_MIXED" for r in rows),
                   nuclear_test=bool(a.nuclear))
    with open(st.path("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    st.commit(dict(step="organelle", label=a.label, fasta=os.path.abspath(a.fasta), fasta_sha256=sha256_file(a.fasta),
                   agp=os.path.abspath(a.agp) if a.agp else None, agp_sha256=sha256_file(a.agp) if a.agp else None,
                   organelles=os.path.abspath(a.organelles), organelles_sha256=sha256_file(a.organelles),
                   organelle_ids=org_ids, plastid_ids=sorted(plastid),
                   nuclear=os.path.abspath(a.nuclear) if a.nuclear else None,
                   nuclear_sha256=sha256_file(a.nuclear) if a.nuclear else None,
                   minimap2=mm, minimap2_version=tool_version([mm, "--version"]),
                   commands=[" ".join(cmd1)] + ([" ".join(cmd2)] if cmd2 else []), science=sci, summary=summary))
    log(f"{a.label}: {summary['confident']} confident organelle, {summary['dominant_mixed']} mixed")


def cmd_plan(a):
    if not a.fcs and not a.organelle:
        sys.exit("give at least one --fcs or --organelle directory")
    cand = collections.defaultdict(list)     # (unit, record) -> [(reason_class, evidence, placement, object)]
    review = []
    inputs = []
    for d in a.fcs or []:
        v = require_validated(d, "fcs")
        u = v["label"]
        inputs.append(dict(step="fcs", label=u, dir=os.path.abspath(d), fasta_sha256=v["fasta_sha256"]))
        for c in read_tsv(os.path.join(d, "fcs_calls.tsv")):
            ev = (f"{c['tool']} {c['action']} {c['start']}-{c['end']} of {c['seq_len']} bp; div {c['div'] or '-'}; "
                  f"cov {c['agg_cont_cov'] or '-'}; {c['top_tax_name'] or '-'}")
            if c["call_class"] == "foreign_contaminant" and c["whole_record"] == "True":
                cand[(u, c["seq_id"])].append(("FOREIGN_CONTAMINANT", ev, c["placement"], c["object"]))
            else:
                review.append(dict(unit=u, record=c["seq_id"], placement=c["placement"], source=c["tool"],
                                   call=c["call_class"], reason="partial or ambiguous call; not a whole-record exclusion",
                                   evidence=ev))
    for d in a.organelle or []:
        v = require_validated(d, "organelle")
        u = v["label"]
        inputs.append(dict(step="organelle", label=u, dir=os.path.abspath(d), fasta_sha256=v["fasta_sha256"]))
        nref = len(v["organelle_ids"])
        for r in read_tsv(os.path.join(d, "organelle_screen.tsv")):
            if not r["window_class"]:
                continue
            nuc = r["nonorganelle_nuclear_bp"]
            ev = (f"organelle cov {r['organelle_frac']} ({nref} organelle refs, >= "
                  f"{v['science']['organelle']['min_identity']:.0%} id, {v['science']['organelle']['window_bp']} bp windows)"
                  + (f", plastid {r['plastid_bp']} bp, mitochondrion {r['mitochondrion_bp']} bp" if r["plastid_bp"] else "")
                  + f", non-organelle {r['nonorganelle_bp']} bp"
                  + (f" (nuclear-aligned {nuc} bp)" if nuc != "" else " (no nuclear test)"))
            if r["window_class"] == "CONFIDENT_ORGANELLE":
                cand[(u, r["record"])].append((r["organelle_kind"], ev, r["placement"], r["object"]))
            else:
                review.append(dict(unit=u, record=r["record"], placement=r["placement"], source="organelle screen",
                                   call=r["window_class"], reason="organelle-dominant but not confident", evidence=ev))
    plan = []
    for (u, rec), calls in sorted(cand.items()):
        classes = sorted({c[0] for c in calls})
        st, obj = calls[0][2], calls[0][3]
        ev = " | ".join(c[1] for c in calls)
        if len(classes) > 1:
            review.append(dict(unit=u, record=rec, placement=st, source="plan", call="+".join(classes),
                               reason="conflicting whole-record calls", evidence=ev))
        elif st in ("unplaced", "no_agp"):
            plan.append(dict(unit=u, record=rec, reason_class=classes[0], evidence=ev))
        else:
            why = {"placed": "placed inside a scaffold (exclusion would be a structural edit)",
                   "unplaced_nonstandard": "unplaced, but not a whole + record under its own name in the AGP",
                   "absent": "not in the AGP"}.get(st, st)
            review.append(dict(unit=u, record=rec, placement=st, source="plan", call=classes[0], reason=why, evidence=ev))
    st = Staged(a.out)
    write_tsv(st.path("EXCLUSION_PLAN.candidate.tsv"), plan, PLAN_COLS, comment=[
        "CANDIDATE exclusion plan from g2f_decontam.py plan. Every row is a proposal: excluding a record is a human",
        "decision. Review it with REVIEW.tsv, then apply with g2f_release_exclude.py check-plan / apply and verify."])
    write_tsv(st.path("REVIEW.tsv"), review, ["unit", "record", "placement", "source", "call", "reason", "evidence"])
    counts = collections.Counter((r["unit"], r["reason_class"]) for r in plan)
    summary = dict(plan_rows=len(plan), review_rows=len(review),
                   by_unit_class=[dict(unit=u, reason_class=c, records=n) for (u, c), n in sorted(counts.items())])
    with open(st.path("summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    st.commit(dict(step="plan", inputs=inputs, summary=summary))
    log(f"plan: {len(plan)} candidate rows, {len(review)} for review")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check-env")
    p.add_argument("--fcs-dir")
    p.set_defaults(fn=cmd_check_env)
    for name, fn in (("fcs", cmd_fcs), ("fcs-summarize", cmd_fcs_summarize)):
        p = sub.add_parser(name)
        p.add_argument("--fasta", required=True)
        p.add_argument("--label", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--agp", help="scaffold AGP: records are labelled placed / unplaced (needed for plan rows)")
        p.add_argument("--science")
        p.set_defaults(fn=fn)
        if name == "fcs":
            p.add_argument("--tax-id", type=int, required=True, help="NCBI taxonomy ID of the species (maize: 4577)")
            p.add_argument("--gx-db", help="FCS-GX database prefix, e.g. /path/gxdb/all (ideally a copy in /dev/shm)")
            p.add_argument("--fcs-dir", help="directory with fcs.py, run_fcsadaptor.sh and images (or G2F_FCS_DIR)")
            p.add_argument("--fcs-python", help="python for fcs.py (default: this interpreter)")
            p.add_argument("--adaptor-image")
            p.add_argument("--gx-image")
            p.add_argument("--container-engine", default="singularity")
            p.add_argument("--threads", type=int, default=48)
            p.add_argument("--skip-gx", action="store_true", help="FCS-adaptor only")
        else:
            p.add_argument("--adaptor-dir", required=True)
            p.add_argument("--gx-dir", help="FCS-GX output directory (omit if GX was not run)")
    p = sub.add_parser("organelle")
    p.add_argument("--fasta", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--organelles", required=True, help="FASTA of the species' plastid and mitochondrial genomes")
    p.add_argument("--plastid-ids", help="comma-separated record IDs in --organelles that are plastid (others: mt)")
    p.add_argument("--nuclear", help="nuclear reference FASTA (e.g. chromosomes only) for the nuclear test")
    p.add_argument("--agp", help="scaffold AGP: screen only unplaced singleton records (default: all records)")
    p.add_argument("--all-records", action="store_true", help="with --agp, screen placed records too")
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--minimap2")
    p.add_argument("--science")
    p.set_defaults(fn=cmd_organelle)
    p = sub.add_parser("plan")
    p.add_argument("--fcs", nargs="+")
    p.add_argument("--organelle", nargs="+")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_plan)
    a = ap.parse_args(argv)
    if a.cmd == "fcs" and not a.skip_gx and not a.gx_db:
        ap.error("fcs: --gx-db is required unless --skip-gx")
    a.fn(a)


if __name__ == "__main__":
    main()
