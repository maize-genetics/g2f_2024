#!/usr/bin/env python3
"""AGP-level scaffold revision from an explicit, reviewed change plan. Contigs are never edited.

Subcommands:
  orient-evidence  for chosen chromosome ends: align the leading (or trailing) components, in scaffold orientation, to a
                   reference (ideally telomere-to-telomere) and report each component's strand relative to the scaffold's
                   overall direction; the maximal reversed run at the tip is proposed as a REVERSE_ORIENTATION plan row
                   (one component) or MANUAL_REVIEW_BLOCK (several: order and orientation would both change)
  check-plan       validate a change plan against an AGP without writing anything
  apply            write the revised AGP, the scaffold FASTA rebuilt from the contigs (bgzip + .fai/.gzi), chain, gap BED,
                   contig-to-scaffold map and a line-level change manifest, into a new directory (built as .partial;
                   BUILD.json records inputs and outputs). Verify it with g2f_agp_revise_verify.py before any use.

Plan TSV columns: action, object, component, old_orientation, new_orientation, new_status, evidence_class, evidence
(an optional `unit` column selects rows with --unit). Actions:
  REVERSE_ORIENTATION          flip the orientation of one component; coordinates and every other line unchanged
  UNPLACE_TERMINAL_EXTENSION   remove trailing components (each with the gap before it) from the end of an object; each
                               becomes its own single-component object (ID and sequence unchanged), inserted in ID order
                               among the unplaced objects if those are in ID order, otherwise appended
Deciding a plan is a human decision (structural edits, docs/workflow_guide.md); this tool only applies and records it.

  g2f_agp_revise.py orient-evidence --agp X.agp --scaffolds X.fa.gz --reference REF.fa --objects-regex REGEX --end start
                                    --out evidence.tsv [--reference-region chr2:1-8000000] [--telomere-integrated DIR]
  g2f_agp_revise.py check-plan --plan plan.tsv --agp X.agp [--unit U]
  g2f_agp_revise.py apply --plan plan.tsv --agp X.agp --contigs X.contigs.fa.gz --out DIR --prefix NAME [--unit U]
                          [--previous-map map.tsv] [--join-table joins.tsv] [--revision-note TEXT] [--threads 8]
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
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")
ACTIONS = ("REVERSE_ORIENTATION", "UNPLACE_TERMINAL_EXTENSION")
with open(os.path.abspath(__file__), "rb") as _fh:
    SCRIPT_SHA256 = hashlib.sha256(_fh.read()).hexdigest()


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def write_tsv(p, rows, keys=None):
    if keys is None:
        keys = []
        for r in rows:
            keys += [k for k in r if k not in keys]
    with open(p, "w") as fh:
        w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def read_agp(path):
    head, body = [], []
    with open(path) as fh:
        for l in fh:
            if l.startswith("#"):
                head.append(l)
            elif l.strip():
                body.append(l.rstrip("\n").split("\t"))
    return head, body


def load_plan(path, unit):
    rows = read_tsv(path)
    if unit is not None:
        if rows and "unit" not in rows[0]:
            sys.exit("--unit given but the plan has no unit column")
        rows = [r for r in rows if r["unit"] == unit]
    for r in rows:
        if r["action"] not in ACTIONS:
            sys.exit(f"unknown action {r['action']!r} (allowed: {ACTIONS})")
    return rows


def check(body, plan):
    """Return a list of problems (empty = the plan applies cleanly)."""
    probs = []
    w = {(x[0], x[5]): x for x in body if x[4] not in ("N", "U")}
    comps = collections.Counter(x[5] for x in body if x[4] not in ("N", "U"))
    seen = set()
    for p in plan:
        k = (p["object"], p["component"])
        if k in seen:
            probs.append(f"duplicate plan row {k}")
        seen.add(k)
        if comps[p["component"]] != 1:
            probs.append(f"{p['component']}: used {comps[p['component']]} times in the AGP (must be exactly once)")
        x = w.get(k)
        if not x:
            probs.append(f"{k}: not in the AGP")
            continue
        if p["action"] == "REVERSE_ORIENTATION":
            if x[8] != p["old_orientation"] or p["new_orientation"] not in "+-" or p["new_orientation"] == x[8]:
                probs.append(f"{k}: orientation {x[8]} in AGP, plan {p['old_orientation']}->{p['new_orientation']}")
        else:
            if x[6] != "1":
                probs.append(f"{k}: component used partially ({x[6]}-{x[7]}); unplacing it would change its sequence")
    # unplaced components must form the tail of their object
    by = collections.defaultdict(list)
    for x in body:
        by[x[0]].append(x)
    un = collections.defaultdict(set)
    for p in plan:
        if p["action"] == "UNPLACE_TERMINAL_EXTENSION":
            un[p["object"]].add(p["component"])
    for o, cs in un.items():
        lines = by.get(o, [])
        idx = [i for i, x in enumerate(lines) if x[4] not in ("N", "U") and x[5] in cs]
        if not idx:
            continue
        first = min(idx)
        tail = lines[first:]
        if any(x[4] not in ("N", "U") and x[5] not in cs for x in tail):
            probs.append(f"{o}: components to unplace are not the terminal run of the object")
        if first == 0:
            probs.append(f"{o}: unplacing every component would leave an empty object")
        elif lines[first - 1][4] not in ("N", "U"):
            probs.append(f"{o}: no gap line before {lines[first][5]}")
    return probs


def apply_plan(head, body, plan, note):
    probs = check(body, plan)
    if probs:
        sys.exit("plan does not apply:\n  " + "\n  ".join(probs))
    order, objs = [], collections.OrderedDict()
    for x in body:
        if x[0] not in objs:
            objs[x[0]] = []
            order.append(x[0])
        objs[x[0]].append(list(x))
    changes = []
    flips = {(p["object"], p["component"]): p for p in plan if p["action"] == "REVERSE_ORIENTATION"}
    for o in order:
        for x in objs[o]:
            if (o, x[5]) in flips:
                old = "\t".join(x)
                x[8] = flips[(o, x[5])]["new_orientation"]
                changes.append(dict(change="orientation", object=o, part_number=x[3], component=x[5],
                                    old_line=old, new_line="\t".join(x)))
    un = {(p["object"], p["component"]) for p in plan if p["action"] == "UNPLACE_TERMINAL_EXTENSION"}
    new_objs = []
    for o in list(order):
        idx = [i for i, x in enumerate(objs[o]) if x[4] not in ("N", "U") and (o, x[5]) in un]
        if not idx:
            continue
        first = min(idx) - 1                       # include the gap before the first removed component
        for x in objs[o][first:]:
            is_gap = x[4] in ("N", "U")
            changes.append(dict(change="removed from object", object=o, part_number=x[3],
                                component=f"gap {x[5]} N" if is_gap else x[5], old_line="\t".join(x), new_line=""))
            if not is_gap:
                L = int(x[7])
                nl = [x[5], "1", str(L), "1", x[4], x[5], "1", str(L), "+"]
                new_objs.append(nl)
                changes.append(dict(change="new unplaced object", object=x[5], part_number="1", component=x[5],
                                    old_line="", new_line="\t".join(nl)))
        objs[o] = objs[o][:first]
    existing_unplaced = [o for o in order if len(objs[o]) == 1 and objs[o][0][5] == o]
    placed = [o for o in order if o not in existing_unplaced]
    for nl in new_objs:
        if nl[0] in objs:
            sys.exit(f"object name {nl[0]} already exists")
        objs[nl[0]] = [nl]
    if existing_unplaced == sorted(existing_unplaced) and order == placed + existing_unplaced:
        out_order = placed + sorted(existing_unplaced + [nl[0] for nl in new_objs])
    else:
        out_order = order + [nl[0] for nl in new_objs]
    lines = list(head) + [f"# REVISION: {note}\n"] if note else list(head)
    for o in out_order:
        for x in objs[o]:
            lines.append("\t".join(x) + "\n")
    return lines, changes, out_order, objs


def fasta_dict(path):
    seqs, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def cmd_check_plan(a):
    head, body = read_agp(a.agp)
    plan = load_plan(a.plan, a.unit)
    probs = check(body, plan)
    print(f"{len(plan)} plan rows; " + ("OK: applies cleanly" if not probs else f"{len(probs)} problems:\n  " + "\n  ".join(probs)))
    return 0 if not probs else 1


def cmd_apply(a):
    out = os.path.abspath(a.out)
    part = out + ".partial"
    for p in (out, part):
        if os.path.exists(p):
            sys.exit(f"refusing to overwrite {p}")
    head, body = read_agp(a.agp)
    plan = load_plan(a.plan, a.unit)
    if not plan:
        sys.exit("empty plan for this unit; nothing to apply (copy the inputs instead)")
    lines, changes, order, objs = apply_plan(head, body, plan, a.revision_note)
    os.makedirs(part)
    P = lambda ext: os.path.join(part, f"{a.prefix}{ext}")
    with open(P(".agp"), "w") as fh:
        fh.writelines(lines)
    seqs = fasta_dict(a.contigs)
    lens = {}
    fa = P(".fa.gz")
    with open(fa, "wb") as fo:
        bg = subprocess.Popen([a.bgzip, "-@", str(a.threads), "-c"], stdin=subprocess.PIPE, stdout=fo)
        for o in order:
            parts = []
            for x in objs[o]:
                if x[4] in ("N", "U"):
                    parts.append("N" * int(x[5]))
                else:
                    seg = seqs[x[5]][int(x[6]) - 1:int(x[7])]
                    parts.append(seg.translate(COMP)[::-1] if x[8] == "-" else seg)
            s = "".join(parts)
            if len(s) != int(objs[o][-1][2]):
                sys.exit(f"{o}: built length {len(s)} != AGP {objs[o][-1][2]}")
            lens[o] = len(s)
            bg.stdin.write(f">{o}\n".encode())
            for i in range(0, len(s), a.line_width):
                bg.stdin.write((s[i:i + a.line_width] + "\n").encode())
        bg.stdin.close()
        if bg.wait() != 0:
            sys.exit("bgzip failed")
    subprocess.run([a.samtools, "faidx", fa], check=True)
    n = 0
    with open(P(".chain"), "w") as ch, open(P(".gaps.bed"), "w") as gb:
        for o in order:
            for x in objs[o]:
                if x[4] in ("N", "U"):
                    gb.write(f"{o}\t{int(x[1]) - 1}\t{x[2]}\n")
                    continue
                n += 1
                cl, cs, ce = len(seqs[x[5]]), int(x[6]), int(x[7])
                S, E, SL = int(x[1]), int(x[2]), lens[o]
                qs, qe = (S - 1, E) if x[8] != "-" else (SL - E, SL - S + 1)
                ch.write(f"chain {ce - cs + 1} {x[5]} {cl} + {cs - 1} {ce} {o} {SL} {'-' if x[8] == '-' else '+'} {qs} {qe} {n}\n{ce - cs + 1}\n\n")
    prev = {r["contig"]: r for r in read_tsv(a.previous_map)} if a.previous_map else {}
    pl = {p["component"]: p for p in plan}
    rows = []
    for o in order:
        k = 0
        for x in objs[o]:
            if x[4] in ("N", "U"):
                continue
            k += 1
            r = dict(prev.get(x[5], {}))
            single = len([y for y in objs[o] if y[4] not in ("N", "U")]) == 1 and o == x[5]
            p = pl.get(x[5])
            status = (p["new_status"] if p and p.get("new_status") else None) or r.get("status") or ("UNPLACED" if single else "PLACED")
            if p and p["action"] == "UNPLACE_TERMINAL_EXTENSION":
                status = p.get("new_status") or "UNPLACED"
            r.update(contig=x[5], contig_length=len(seqs[x[5]]), scaffold=o, scaffold_start=x[1], scaffold_end=x[2],
                     orientation=x[8], order_in_scaffold=k, status=status)
            if p:
                r["revision_note"] = f"{p['action']} ({p.get('evidence_class', '')})"
            rows.append(r)
    base = ["contig", "contig_length", "scaffold", "scaffold_start", "scaffold_end", "orientation", "order_in_scaffold", "status"]
    write_tsv(P(".contig_to_scaffold.tsv"), rows, base + [k for k in (list(next(iter(prev.values()))) if prev else []) if k not in base] + ["revision_note"])
    cm = []
    for c in changes:
        p = pl.get(c["component"], {})
        o_, n_ = (c["old_line"].split("\t") if c["old_line"] else [""] * 9), (c["new_line"].split("\t") if c["new_line"] else [""] * 9)
        cm.append(dict(change=c["change"], object_old=o_[0], object_new=n_[0], part_number_old=o_[3], part_number_new=n_[3],
                       component=c["component"], coords_old=f"{o_[1]}-{o_[2]}" if o_[0] else "", coords_new=f"{n_[1]}-{n_[2]}" if n_[0] else "",
                       orientation_old=o_[8], orientation_new=n_[8], evidence_class=p.get("evidence_class", ""),
                       new_status=p.get("new_status", ""), agp_line_old=c["old_line"], agp_line_new=c["new_line"]))
    write_tsv(P(".change_manifest.tsv"), cm)
    if a.join_table:
        js = read_tsv(a.join_table)
        gone = {p["component"] for p in plan if p["action"] == "UNPLACE_TERMINAL_EXTENSION"}
        flip = {p["component"]: p for p in plan if p["action"] == "REVERSE_ORIENTATION"}
        kept, removed = [], []
        for r in js:
            if r["contig_A"] in gone or r["contig_B"] in gone:
                removed.append(r)
                continue
            for s in ("A", "B"):
                if r[f"contig_{s}"] in flip:
                    if f"{s}_orientation" in r:
                        r[f"{s}_orientation"] = flip[r[f"contig_{s}"]]["new_orientation"]
                    r["evidence_call"] = "ORIENTATION_REVISED_NOT_RETESTED"
                    r["revision_note"] = f"{r[f'contig_{s}']} reversed; the new contig-end adjacency was not re-tested with reads"
            kept.append(r)
        write_tsv(P(".join_table.tsv"), kept)
        write_tsv(P(".join_table_removed.tsv"), removed, list(js[0]) if js else ["contig_A"])
    del seqs
    files = {f: sha256_file(os.path.join(part, f)) for f in sorted(os.listdir(part))}
    rec = {"schema": "g2f-agp-revise/1", "step": "apply", "state": "built (not yet verified)", "prefix": a.prefix,
           "inputs": {"agp": os.path.abspath(a.agp), "agp_sha256": sha256_file(a.agp), "contigs": os.path.abspath(a.contigs),
                      "contigs_sha256": sha256_file(a.contigs), "plan": os.path.abspath(a.plan), "plan_sha256": sha256_file(a.plan),
                      "unit": a.unit, "previous_map": a.previous_map, "join_table": a.join_table},
           "plan_rows": plan, "agp_line_changes": len(changes), "records": len(order), "bp": sum(lens.values()),
           "script": os.path.abspath(__file__), "script_sha256": SCRIPT_SHA256, "python": sys.version.split()[0],
           "built": datetime.datetime.now().astimezone().isoformat(timespec="seconds"), "outputs_sha256": files,
           "next_step": "g2f_agp_revise_verify.py (independent); use the outputs only if it passes"}
    with open(os.path.join(part, "BUILD.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    os.rename(part, out)
    print(f"built {out}: {len(changes)} AGP line changes, {len(order)} records")
    return 0


def cmd_orient_evidence(a):
    if os.path.exists(a.out):
        sys.exit(f"refusing to overwrite {a.out}")
    for f in (a.scaffolds + ".fai",) + ((a.scaffolds + ".gzi",) if a.scaffolds.endswith(".gz") else ()):
        if not os.path.exists(f):
            sys.exit(f"missing index {f}; create it yourself (samtools faidx writes it next to the input)")
    head, body = read_agp(a.agp)
    rx = re.compile(a.objects_regex)
    by = collections.OrderedDict()
    for x in body:
        if rx.search(x[0]):
            by.setdefault(x[0], []).append(x)
    tel = {}
    if a.telomere_integrated:
        for r in read_tsv(os.path.join(a.telomere_integrated, "contig_ends.tsv")):
            tel[(r["contig"], r["side"])] = r["class"]
    tmp = tempfile.mkdtemp(prefix="g2f_orient_", dir=a.tmpdir)
    ref = a.reference
    if a.reference_region:
        ref = os.path.join(tmp, "ref.fa")
        with open(ref, "w") as fh:
            fh.write(subprocess.run([a.samtools, "faidx", a.reference, a.reference_region], capture_output=True, text=True, check=True).stdout)
    rows, plan = [], []
    for o, xs in by.items():
        W = [x for x in xs if x[4] not in ("N", "U")]
        if a.end == "end":
            W = W[::-1]
        L = int(xs[-1][2])
        lead = [c for i, c in enumerate(W) if i < 2 or (int(c[1]) <= a.window if a.end == "start" else L - int(c[2]) + 1 <= a.window)]
        q = os.path.join(tmp, "q.fa")
        with open(q, "w") as fh:
            for c in lead:
                s, e = int(c[1]), int(c[2])
                if a.end == "start":
                    e = min(e, max(a.window, s + 1_000_000))
                else:
                    s = max(s, min(L - a.window + 1, e - 1_000_000))
                seq = subprocess.run([a.samtools, "faidx", a.scaffolds, f"{o}:{s}-{e}"], capture_output=True, text=True, check=True).stdout
                fh.write(f">{c[5]}\n" + "".join(seq.splitlines()[1:]) + "\n")
        paf = subprocess.run([a.minimap2, "-t", str(a.threads), "-x", "asm20", "--secondary=no", ref, q],
                             capture_output=True, text=True, check=True).stdout
        al = collections.defaultdict(lambda: {"+": 0, "-": 0})
        for l in paf.splitlines():
            x = l.split("\t")
            if int(x[11]) >= a.min_mapq:
                al[x[0]][x[4]] += int(x[3]) - int(x[2])
        strand = {}
        for c in lead:
            g = al[c[5]]
            t = g["+"] + g["-"]
            strand[c[5]] = None if t == 0 else ("+" if g["+"] >= 0.7 * t else "-" if g["-"] >= 0.7 * t else "?")
        inner = [c for c in lead[1:]]
        tot = collections.Counter()
        for c in inner:
            for s_ in "+-":
                tot[s_] += al[c[5]][s_]
        overall = "+" if tot["+"] >= tot["-"] else "-"
        block, open_ = [], True
        for c in lead:
            rev = strand[c[5]] is not None and strand[c[5]] not in ("?", overall)
            open_ = open_ and rev
            if open_:
                block.append(c)
            g = al[c[5]]
            rows.append(dict(object=o, end=a.end, component=c[5], orientation=c[8], scaffold_start=c[1], scaffold_end=c[2],
                             component_length=int(c[7]) - int(c[6]) + 1, aligned_plus=g["+"], aligned_minus=g["-"],
                             component_strand_vs_reference=strand[c[5]] or "NO_ALIGNMENT", scaffold_overall_strand=overall,
                             reversed_relative_to_scaffold=rev, in_tip_block=open_ and rev,
                             telomere_contig_L=tel.get((c[5], "L"), ""), telomere_contig_R=tel.get((c[5], "R"), "")))
        if len(block) == 1:
            c = block[0]
            outward_after = ("R" if c[8] == "+" else "L") if a.end == "start" else ("L" if c[8] == "+" else "R")
            tcls = tel.get((c[5], outward_after), "")
            ev = "TELOMERE_AND_SYNTENY" if tcls in ("TEL_KB", "TEL_SHORT") else "SYNTENY_ONLY"
            plan.append(dict(action="REVERSE_ORIENTATION", object=o, component=c[5], old_orientation=c[8],
                             new_orientation="-" if c[8] == "+" else "+", new_status="PLACED", evidence_class=ev,
                             evidence=f"tip component reversed vs reference ({al[c[5]]['-' if overall == '+' else '+']} of "
                                      f"{al[c[5]]['+'] + al[c[5]]['-']} aligned bp on the opposite strand)"
                                      + (f"; contig end that becomes outward carries {tcls}" if tcls else "")))
        elif len(block) > 1:
            plan.append(dict(action="MANUAL_REVIEW_BLOCK", object=o, component=",".join(c[5] for c in block), old_orientation="",
                             new_orientation="", new_status="", evidence_class="SYNTENY",
                             evidence="several leading components reversed: order and orientation would both change; not applied by this tool"))
    shutil.rmtree(tmp)                                   # this step's own temporary files
    write_tsv(a.out, rows)
    if plan:
        write_tsv(a.out + ".proposed_plan.tsv", plan)
    print(f"{len(rows)} components examined; {len(plan)} proposed plan rows"
          + (f" -> {a.out}.proposed_plan.tsv (review before use)" if plan else ""))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check-plan")
    p.add_argument("--plan", required=True)
    p.add_argument("--agp", required=True)
    p.add_argument("--unit")
    p.set_defaults(func=cmd_check_plan)
    p = sub.add_parser("apply")
    for k in ("--plan", "--agp", "--contigs", "--out", "--prefix"):
        p.add_argument(k, required=True)
    p.add_argument("--unit")
    p.add_argument("--previous-map")
    p.add_argument("--join-table")
    p.add_argument("--revision-note", default="AGP-level revision; see the change manifest")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--line-width", type=int, default=80)
    p.add_argument("--bgzip", default=os.environ.get("G2F_BGZIP", "bgzip"))
    p.add_argument("--samtools", default=os.environ.get("G2F_SAMTOOLS", "samtools"))
    p.set_defaults(func=cmd_apply)
    p = sub.add_parser("orient-evidence")
    for k in ("--agp", "--scaffolds", "--reference", "--objects-regex", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--end", choices=("start", "end"), required=True)
    p.add_argument("--reference-region")
    p.add_argument("--telomere-integrated")
    p.add_argument("--window", type=int, default=3_000_000)
    p.add_argument("--min-mapq", type=int, default=5)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--tmpdir")
    p.add_argument("--minimap2", default=os.environ.get("G2F_MINIMAP2", "minimap2"))
    p.add_argument("--samtools", default=os.environ.get("G2F_SAMTOOLS", "samtools"))
    p.set_defaults(func=cmd_orient_evidence)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
