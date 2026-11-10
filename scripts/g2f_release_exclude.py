#!/usr/bin/env python3
"""Remove whole records (e.g. foreign-contaminant or organellar contigs) from a contig FASTA and its scaffold
representation, from an explicit, reviewed exclusion plan.

Only records that are UNPLACED SINGLETONS can be excluded: the record is the sole component of an AGP object of the
same name, used whole and in + orientation. Anything placed inside a scaffold is refused (that would be a structural
edit; see docs/workflow_guide.md, stage 8). Kept records are copied byte for byte, in their original order.

Subcommands:
  check-plan  --plan plan.tsv --agp X.agp --contigs X.contigs.fa.gz [--unit U]
              validate the plan; write nothing
  apply       --plan plan.tsv --unit U --contigs X.contigs.fa.gz --scaffolds X.scaffolds.fa.gz --agp X.agp
              --chain X.chain --map X.contig_to_scaffold.tsv --gaps-bed X.gaps.bed --out DIR
              --contigs-name NAME.fa.gz --scaffolds-name NAME.fa.gz --prefix NAME [--threads 8]
              [--agp-header TEMPLATE] [--relabel-linkage OLD=NEW]
              writes DIR/{contigs/,scaffolds_rg/,provenance/}, DIR/<prefix>.EXCLUDED_RECORDS.tsv and DIR/BUILD.json
              (built as DIR.partial, renamed on success). Verify with g2f_release_exclude_verify.py before any use.

AGP metadata (optional, both recorded in BUILD.json and checked by the verifier):
  --agp-header TEMPLATE   replace the comment lines at the top of the AGP with TEMPLATE. Its first line must be
                          '##agp-version 2.1' (tab or space). Placeholders: {contigs_name}, {contigs_sha256} (sha256 of
                          the uncompressed output contig FASTA), {scaffolds_name}, {excluded}. Without it, comments are
                          copied unchanged (and may then describe the previous version).
  --relabel-linkage OLD=NEW  on gap lines (N/U), replace the linkage-evidence token OLD by NEW (AGP 2.1 vocabulary).
                          Only that token changes; every other column and every component line is copied unchanged.

Plan TSV columns: unit, record, reason_class, evidence. The decision to exclude a record is a human decision; this
tool only applies and records it. The EXCLUDED_RECORDS table gives, for each removed record, its length and the
sha256 of its sequence (uppercase, no newlines), so downstream features on those sequence IDs can be found and dropped.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import subprocess
import sys

with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()
PLAN_COLS = ("unit", "record", "reason_class", "evidence")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def load_plan(path, unit):
    rows = read_tsv(path)
    if rows and any(c not in rows[0] for c in PLAN_COLS):
        sys.exit(f"plan must have columns {PLAN_COLS}")
    if unit is not None:
        rows = [r for r in rows if r["unit"] == unit]
    seen = set()
    for r in rows:
        k = (r["unit"], r["record"])
        if k in seen:
            sys.exit(f"duplicate plan row {k}")
        seen.add(k)
        if not r["reason_class"].strip():
            sys.exit(f"plan row {k} has no reason_class")
    return rows


def read_agp(path):
    with open(path) as fh:
        return [l for l in fh]


def agp_objects(lines):
    objs = {}
    for l in lines:
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        objs.setdefault(f[0], []).append(f)
    return objs


def fai_lengths(fa):
    fai = fa + ".fai"
    if not os.path.exists(fai):
        sys.exit(f"missing {fai}")
    with open(fai) as fh:
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in fh}


def problems(plan, agp_lines, contig_lens):
    probs = []
    objs = agp_objects(agp_lines)
    used_in = {}
    for o, rows in objs.items():
        for f in rows:
            if f[4] not in ("N", "U"):
                used_in.setdefault(f[5], []).append(o)
    for r in plan:
        c = r["record"]
        if c not in contig_lens:
            probs.append(f"{c}: not in the contig FASTA")
            continue
        where = used_in.get(c, [])
        if where != [c]:
            probs.append(f"{c}: not an unplaced singleton (component of {where or 'no object'}); refusing")
            continue
        rows = objs[c]
        f = rows[0]
        if len(rows) != 1 or f[4] in ("N", "U") or (f[6], f[7]) != ("1", str(contig_lens[c])) or f[8] != "+" \
                or (f[1], f[2]) != ("1", str(contig_lens[c])):
            probs.append(f"{c}: object {c} is not exactly one whole + component; refusing")
    return probs


def filter_fasta(src, dst, drop, threads, bgzip, samtools):
    """Stream src -> dst, dropping records in `drop`; kept lines are copied byte for byte.
    Returns {dropped_id: (length, sha256_upper)} and kept/dropped counts."""
    info, keep, cur, h, n = {}, True, None, None, 0
    kept = dropped = 0
    whole = hashlib.sha256()
    dec = subprocess.Popen([bgzip, "-dc", "-@", str(max(1, threads // 4)), src], stdout=subprocess.PIPE)
    with open(dst, "wb") as fo:
        enc = subprocess.Popen([bgzip, "-@", str(threads), "-c"], stdin=subprocess.PIPE, stdout=fo)

        def close_rec():
            if cur is not None and not keep:
                info[cur] = (n, h.hexdigest())
        for line in dec.stdout:
            if line.startswith(b">"):
                close_rec()
                cur = line[1:].split()[0].decode()
                keep = cur not in drop
                kept += keep
                dropped += not keep
                h, n = hashlib.sha256(), 0
            elif not keep:
                s = line.rstrip(b"\r\n").upper()
                h.update(s)
                n += len(s)
                continue
            if keep:
                enc.stdin.write(line)
                whole.update(line)
        close_rec()
        enc.stdin.close()
        if enc.wait() != 0 or dec.wait() != 0:
            sys.exit(f"bgzip failed on {src}")
    subprocess.run([samtools, "faidx", dst], check=True)
    missing = set(drop) - set(info)
    if missing:
        sys.exit(f"{src}: planned records not found: {sorted(missing)}")
    return info, kept, dropped, whole.hexdigest()


LINKAGE_EVIDENCE = {"na", "paired-ends", "align_genus", "align_xgenus", "align_trnscpt", "within_clone", "clone_contig",
                    "map", "pcr", "proximity_ligation", "strobe", "unspecified"}


def parse_relabel(spec):
    if spec is None:
        return None
    old, sep, new = spec.partition("=")
    if not sep or old not in LINKAGE_EVIDENCE or new not in LINKAGE_EVIDENCE or old == new or "na" in (old, new):
        sys.exit(f"--relabel-linkage {spec!r}: need OLD=NEW, two different AGP 2.1 linkage-evidence terms other than na")
    return old, new


def rewrite_agp(lines, header, relabel):
    """header: list of rendered comment lines or None; relabel: (old, new) or None. Returns (lines, n_relabelled)."""
    k = 0
    while k < len(lines) and lines[k].startswith("#"):
        k += 1
    if any(l.startswith("#") for l in lines[k:]):
        sys.exit("AGP comment lines must all precede the data lines")
    body = lines[k:]
    out, n = [], 0
    for l in body:
        f = l.rstrip("\n").split("\t")
        if relabel and len(f) == 9 and f[4] in ("N", "U"):
            ev = f[8].split(";")
            if relabel[0] in ev:
                f[8] = ";".join(relabel[1] if e == relabel[0] else e for e in ev)
                n += 1
                l = "\t".join(f) + "\n"
        out.append(l)
    head = header if header is not None else lines[:k]
    return head + out, n


def filter_lines(src, dst, drop_fn):
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for l in fi:
            if drop_fn(l):
                n += 1
            else:
                fo.write(l)
    return n


def filter_chain(src, dst, drop):
    n = 0
    with open(src) as fh:
        blocks = fh.read().split("\n\n")
    with open(dst, "w") as fo:
        out = []
        for b in blocks:
            if not b.strip():
                continue
            f = b.split("\n", 1)[0].split()
            if f[2] in drop:
                n += 1
                continue
            out.append(b.strip("\n"))
        fo.write("".join(x + "\n\n" for x in out))
    return n


def cmd_check_plan(a):
    plan = load_plan(a.plan, a.unit)
    probs = problems(plan, read_agp(a.agp), fai_lengths(a.contigs))
    print(f"{len(plan)} plan rows; " + ("OK: applies cleanly" if not probs else f"{len(probs)} problems:\n  " + "\n  ".join(probs)))
    return 0 if not probs else 1


def cmd_apply(a):
    out = os.path.abspath(a.out)
    part = out + ".partial"
    for p in (out, part):
        if os.path.exists(p):
            sys.exit(f"refusing to overwrite {p}")
    plan = load_plan(a.plan, a.unit)
    if not plan:
        sys.exit("empty plan for this unit; nothing to apply (copy the inputs instead)")
    agp = read_agp(a.agp)
    probs = problems(plan, agp, fai_lengths(a.contigs))
    if probs:
        sys.exit("plan does not apply:\n  " + "\n  ".join(probs))
    drop = {r["record"] for r in plan}
    for d in ("contigs", "scaffolds_rg", "provenance"):
        os.makedirs(os.path.join(part, d))
    t0 = datetime.datetime.now().isoformat(timespec="seconds")
    c_out = os.path.join(part, "contigs", a.contigs_name)
    s_out = os.path.join(part, "scaffolds_rg", a.scaffolds_name)
    relabel = parse_relabel(a.relabel_linkage)
    template = None
    if a.agp_header:
        with open(a.agp_header) as fh:
            template = fh.read()
        if not template.split("\n", 1)[0].replace("\t", " ").startswith("##agp-version 2.1"):
            sys.exit("--agp-header: first line must be '##agp-version 2.1'")
        if any(l and not l.startswith("#") for l in template.splitlines()):
            sys.exit("--agp-header: every line must be a comment")
    c_info, c_kept, _, c_sha = filter_fasta(a.contigs, c_out, drop, a.threads, a.bgzip, a.samtools)
    s_info, s_kept, _, _ = filter_fasta(a.scaffolds, s_out, drop, a.threads, a.bgzip, a.samtools)
    for k in drop:
        if c_info[k] != s_info[k]:
            sys.exit(f"{k}: scaffold record differs from the contig record ({s_info[k]} vs {c_info[k]})")
    P = lambda name: os.path.join(part, "provenance", name)
    base = lambda p: os.path.basename(p)
    n_agp = filter_lines(a.agp, P(base(a.agp)), lambda l: not l.startswith("#") and l.split("\t", 1)[0] in drop)
    header = None
    if template is not None:
        rendered = template.format(contigs_name=a.contigs_name, contigs_sha256=c_sha, scaffolds_name=a.scaffolds_name,
                                   excluded=len(drop))
        header = [l + "\n" for l in rendered.rstrip("\n").split("\n")]
    n_relabel = 0
    if header is not None or relabel:
        new, n_relabel = rewrite_agp(read_agp(P(base(a.agp))), header, relabel)
        with open(P(base(a.agp)), "w") as fh:
            fh.writelines(new)
    n_map = filter_lines(a.map, P(base(a.map)), lambda l: l.split("\t", 1)[0] in drop)
    n_gap = filter_lines(a.gaps_bed, P(base(a.gaps_bed)), lambda l: l.split("\t", 1)[0] in drop)
    n_chain = filter_chain(a.chain, P(base(a.chain)), drop)
    if (n_agp, n_map, n_chain, n_gap) != (len(drop), len(drop), len(drop), 0):
        sys.exit(f"unexpected removal counts agp={n_agp} map={n_map} chain={n_chain} gaps={n_gap} for {len(drop)} records")
    ex = os.path.join(part, f"{a.prefix}.EXCLUDED_RECORDS.tsv")
    with open(ex, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["unit", "record", "length", "seq_sha256_upper", "reason_class", "evidence",
                    "former_contig_record", "former_scaffold_record"])
        for r in plan:
            L, h = c_info[r["record"]]
            w.writerow([r["unit"], r["record"], L, h, r["reason_class"], r["evidence"],
                        f"{r['record']}:1-{L}", f"{r['record']}:1-{L}"])
    inputs = {k: getattr(a, k) for k in ("plan", "contigs", "scaffolds", "agp", "chain", "map", "gaps_bed")}
    build = {
        "schema": "g2f-release-exclude/1", "tool": "g2f_release_exclude.py", "tool_sha256": SCRIPT_SHA256,
        "started": t0, "finished": datetime.datetime.now().isoformat(timespec="seconds"), "unit": a.unit,
        "inputs": {k: {"path": os.path.abspath(v), "sha256": sha256_file(v)} for k, v in inputs.items()},
        "excluded": len(drop), "excluded_bp": sum(v[0] for v in c_info.values()),
        "contig_records_kept": c_kept, "scaffold_records_kept": s_kept,
        "removed_lines": {"agp": n_agp, "contig_to_scaffold": n_map, "chain": n_chain, "gaps_bed": n_gap},
        "contigs_sha256_uncompressed": c_sha,
        "agp_header": None if template is None else {"template_sha256": hashlib.sha256(template.encode()).hexdigest(),
                                                     "template": template},
        "linkage_relabel": None if not relabel else {"from": relabel[0], "to": relabel[1], "gap_lines": n_relabel},
        "outputs": {os.path.relpath(p, part): sha256_file(p) for p in
                    [c_out, s_out, ex] + [P(base(x)) for x in (a.agp, a.map, a.gaps_bed, a.chain)]},
        "verified": False, "note": "verify with g2f_release_exclude_verify.py before any use",
    }
    with open(os.path.join(part, "BUILD.json"), "w") as fh:
        json.dump(build, fh, indent=1)
    os.rename(part, out)
    print(f"{a.unit}: excluded {len(drop)} records ({build['excluded_bp']} bp); kept {c_kept} contigs, {s_kept} scaffolds -> {out}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check-plan")
    p.add_argument("--plan", required=True)
    p.add_argument("--agp", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--unit")
    p = sub.add_parser("apply")
    for k in ("plan", "unit", "contigs", "scaffolds", "agp", "chain", "map", "gaps-bed", "out", "contigs-name",
              "scaffolds-name", "prefix"):
        p.add_argument(f"--{k}", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--agp-header", help="comment-line template for the output AGP (see above)")
    p.add_argument("--relabel-linkage", help="OLD=NEW linkage-evidence term on gap lines")
    p.add_argument("--bgzip", default=os.environ.get("G2F_BGZIP", "bgzip"))
    p.add_argument("--samtools", default=os.environ.get("G2F_SAMTOOLS", "samtools"))
    a = ap.parse_args(argv)
    return cmd_check_plan(a) if a.cmd == "check-plan" else cmd_apply(a)


if __name__ == "__main__":
    sys.exit(main())
