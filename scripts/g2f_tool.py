#!/usr/bin/env python3
"""Command-line helpers for the G2F workflow.

Most subcommands are called by Snakemake rules. `discover` is also meant for
people: it shows which inputs would be used and why, without processing reads.

  python3 scripts/g2f_tool.py discover            # review inputs (read-only)
  python3 scripts/g2f_tool.py discover --json     # machine-readable inventory
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from g2f import discovery, gfa, kmer, provenance, readstats, summary, validate  # noqa: E402

PIPELINE_DIR = os.path.dirname(HERE)


def _kv(items):
    out = {}
    for item in items or []:
        k, _, v = item.partition("=")
        if not v:
            sys.exit(f"expected KEY=VALUE, got {item!r}")
        out[k] = v
    return out


def _emit(result, out, fail_msg):
    provenance.write_json(out, result)
    if not result.get("ok", True):
        sys.stderr.write(f"{fail_msg}:\n  " + "\n  ".join(result.get("errors", [])) + "\n")
        sys.exit(1)


def cmd_discover(a):
    import yaml
    from g2f.config import ConfigError, load_config
    with open(a.config) as fh:
        top = yaml.safe_load(fh) or {}
    top.update(_kv(a.set))
    if a.only:
        top["only_samples"] = a.only
    try:
        cfg = load_config(top, PIPELINE_DIR)
    except ConfigError as exc:
        sys.exit(f"Configuration error: {exc}")
    inputs, issues = discovery.discover(cfg)
    records = [discovery.inventory_records(s, inputs[s], cfg) for s in sorted(inputs)]
    if a.json:
        print(json.dumps({"samples": records, "issues": [vars(i) for i in issues]}, indent=2))
    else:
        for r in records:
            print(f"\n== {r['sample']}  (aliases: {', '.join(r['aliases']) or '-'})")
            print(f"   method: {r['selected_method']}  [{r['selection_reason']}]")
            for i in r["inputs"]:
                ident = i["identity"]
                desc = ", ".join(f"{k}={v}" for k, v in ident.items()
                                 if k in ("sampled_movies", "readtypes", "sampled_runids",
                                          "sampled_len_min", "basecall_models"))
                print(f"   {i['status']:8s} {i['tech']:4s} {i['format'] or '-':5s} "
                      f"{i['size'] / 1e9:8.2f} GB  {i['path']}")
                if desc:
                    print(f"            {desc}")
                for n in i["notes"]:
                    print(f"            note: {n}")
        if not records:
            print(f"No samples found under {cfg['data_dir']} or in external_inputs.")
        if issues:
            print("\nIssues:")
            for i in issues:
                print(i.render())
    errors = [i for i in issues if i.severity == "error"]
    sys.exit(1 if errors else 0)


def cmd_run(a):
    if not a.command:
        sys.exit("run: missing command after --")
    command = " ".join(a.command) if len(a.command) > 1 else a.command[0]
    rc = provenance.run_step(a.prov, a.step, command, sample=a.sample, method=a.method,
                             version_cmds=_kv(a.version), expect=_kv(a.expect),
                             inputs=a.inputs or [], outputs=a.outputs or [],
                             pipeline_dir=PIPELINE_DIR)
    sys.exit(rc)


def cmd_readstats(a):
    res = readstats.aggregate(sys.stdin, a.tech, a.has_qual, a.genome_size, a.ids_out)
    res["source"] = a.source
    provenance.write_json(a.out, res)


def cmd_check_gfa(a):
    res = gfa.compare(a.gfa, a.fasta)
    res.update(gfa=a.gfa, fasta=a.fasta)
    _emit(res, a.out, "GFA/FASTA consistency check failed")


def cmd_kmer_screen(a):
    res = kmer.screen(kmer.read_histogram(a.hist), (a.window_lo, a.window_hi), a.max_ratio,
                      a.min_local_max_ratio, a.smooth_frac, a.k)
    res["histogram"] = a.hist
    _emit(res, a.out, "k-mer screen could not be computed")


def cmd_kmer_assess(a):
    """Combine the k-mer screen with any recorded human review.

    Default: diagnostic only (a flagged sample gets a recorded warning). With
    --block, a flagged sample without an accepted review fails, which stops
    assembly. An invalid/uncomputable screen always fails."""
    import yaml
    with open(a.screen) as fh:
        scr = json.load(fh)
    with open(a.samples) as fh:
        entry = ((yaml.safe_load(fh) or {}).get("samples") or {}).get(a.sample) or {}
    review = entry.get("kmer_review") or None
    errors, warnings = [], []
    if not scr.get("ok"):
        errors.append(f"k-mer screen could not be computed: {scr.get('error')}")
        status = "invalid"
    elif not scr.get("flagged"):
        status = "not_flagged"
    elif review and review.get("status") == "accepted" and review.get("reason"):
        status = "flagged_reviewed"
    else:
        status = "flagged_unreviewed"
        msg = (f"k-mer screen flagged {a.sample} ({scr.get('flag_reason')}). The screen is an "
               "uncalibrated heuristic; inspect " + scr.get("histogram", "the histogram") +
               " and, if appropriate, record kmer_review: {status: accepted, reason, reviewer, "
               "date} in config/samples.yaml.")
        (errors if a.block else warnings).append(msg)
    res = {"sample": a.sample, "status": status, "blocking": a.block, "screen": scr,
           "review": review, "review_id": a.review_id, "warnings": warnings,
           "ok": not errors, "errors": errors}
    for w in warnings:
        sys.stderr.write("WARNING: " + w + "\n")
    _emit(res, a.out, "Assembly blocked pending human review" if a.block else
          "k-mer assessment failed")


def cmd_method_snapshot(a):
    import yaml
    from g2f.config import method_fingerprint
    with open(a.method_file) as fh:
        m = yaml.safe_load(fh)
    provenance.write_json(a.out, {"method_id": m["method_id"], "fingerprint": method_fingerprint(m),
                                  "method": m, "science_version": a.science_version,
                                  "method_file": a.method_file})


def cmd_import_assembly(a):
    """Record identity of a read-only historical assembly and link it into results."""
    from g2f.discovery import fingerprint
    links = dict(zip(("gfa", "fa", "log"), (a.gfa, a.fa, a.log)))
    targets = dict(zip(("gfa", "fa", "log"), (a.link_gfa, a.link_fa, a.link_log)))
    rec = {"sample": a.sample, "assembly_id": a.aid, "assembly_source": "imported_legacy",
           "statement": "Historical assembly registered read-only; NOT produced by this workflow.",
           "register": a.register, "files": {}}
    for k, src in links.items():
        if not os.path.isfile(src):
            sys.exit(f"registered file not found: {src}")
        size = os.path.getsize(src)
        rec["files"][k] = {"path": src, "size": size, "mtime": os.path.getmtime(src),
                           "fingerprint": fingerprint(src, size)}
        dst = targets[k]
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.lexists(dst):
            os.unlink(dst)          # only ever a symlink created by this command
        os.symlink(src, dst)
    provenance.write_json(a.out, rec)


def cmd_validate_hifi(a):
    res = validate.validate_hifi(_kv(a.stats), _kv(a.bam_count), a.inventory, a.genome_size,
                                 readstats.aggregator_id())
    _emit(res, a.out, "HiFi input validation failed")


def cmd_validate_ont(a):
    res = validate.validate_ont(_kv(a.input_stats), _kv(a.retained_stats), _kv(a.input_ids),
                                _kv(a.retained_ids), a.min_len, a.min_retained, a.inventory,
                                a.genome_size, readstats.aggregator_id())
    _emit(res, a.out, "ONT input validation failed")


def cmd_validate_busco(a):
    pin = json.loads(a.pin) if a.pin else None
    if pin is not None:
        pin = dict(pin, name=a.lineage)
    _emit(validate.validate_busco(a.json, a.dataset_cfg, a.expected_version, pin), a.out,
          "BUSCO validation failed")


def cmd_validate_merqury(a):
    _emit(validate.validate_merqury(a.qv, a.completeness, a.asm_name), a.out,
          "Merqury validation failed")


def cmd_validate_asm(a):
    _emit(validate.validate_assembly_stats(a.stats, a.check), a.out,
          "Assembly statistics validation failed")


def cmd_summarize(a):
    L = summary._load
    opt = lambda p: L(p) if p else None   # noqa: E731
    hlog = validate.parse_hifiasm_log(a.hifiasm_log)
    method_rec = opt(a.method_record)
    if method_rec and a.method_fingerprint and method_rec["fingerprint"] != a.method_fingerprint:
        sys.exit(f"{a.sample}/{a.method}: assembly was produced with method fingerprint "
                 f"{method_rec['fingerprint']}, but the current method definition is "
                 f"{a.method_fingerprint}. Method files are immutable: restore it or create a "
                 "new method_id.")
    imported = opt(a.imported)
    rows = summary.build(a.sample, a.method, a.science_version, opt(a.hifi), opt(a.ont),
                         opt(a.kmer), L(a.asm_p), opt(a.asm_a), [L(b) for b in a.busco or []],
                         opt(a.merqury), hlog,
                         source="imported_legacy" if imported else "pipeline",
                         validation_id=a.validation_id,
                         method_fingerprint=(method_rec or {}).get("fingerprint", ""))
    summary.write_tsv(rows, a.out_tsv)
    provenance.write_json(a.out_json, {
        "ok": True, "sample": a.sample, "method_id": a.method,
        "assembly_source": "imported_legacy" if imported else "pipeline",
        "validation_id": a.validation_id, "science_version": a.science_version,
        "method_record": method_rec, "imported": imported, "rows": rows, "hifiasm_log": hlog,
        "validated_records": [p for p in [a.hifi, a.ont, a.kmer, a.asm_p, a.asm_a, a.merqury,
                                          *(a.busco or [])] if p]})


def cmd_combine(a):
    rows = [r for p in a.inputs for r in summary.read_tsv(p)]
    rows.sort(key=lambda r: (r["sample"], r["method_id"], r.get("busco_lineage", "")))
    summary.write_tsv(rows, a.out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("discover", help="discover and validate inputs (read-only)")
    s.add_argument("--config", default=os.path.join(PIPELINE_DIR, "config", "config.yaml"))
    s.add_argument("--set", action="append", help="override top-level config KEY=VALUE")
    s.add_argument("--only", nargs="*", help="restrict to these samples")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_discover)

    s = sub.add_parser("run", help="run a command and record provenance")
    s.add_argument("--prov", required=True)
    s.add_argument("--step", required=True)
    s.add_argument("--sample")
    s.add_argument("--method")
    s.add_argument("--version", action="append", help="TOOL=VERSION_COMMAND")
    s.add_argument("--expect", action="append", help="TOOL=REQUIRED_VERSION_SUBSTRING")
    s.add_argument("--inputs", nargs="*")
    s.add_argument("--outputs", nargs="*")
    s.add_argument("command", nargs=argparse.REMAINDER)
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("identity", help="write the validation-identity record")
    s.add_argument("--out", required=True)
    s.add_argument("--json", required=True)
    s.set_defaults(func=lambda a: provenance.write_json(a.out, json.loads(a.json)))

    s = sub.add_parser("readstats", help="aggregate `seqkit fx2tab -n -l [-q]` from stdin")
    s.add_argument("--tech", choices=["hifi", "ont"], required=True)
    s.add_argument("--has-qual", action="store_true")
    s.add_argument("--genome-size", type=int)
    s.add_argument("--ids-out")
    s.add_argument("--source")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_readstats)

    s = sub.add_parser("check-gfa")
    s.add_argument("--gfa", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_check_gfa)

    s = sub.add_parser("kmer-screen")
    s.add_argument("--hist", required=True)
    s.add_argument("--window-lo", type=float, default=0.35)
    s.add_argument("--window-hi", type=float, default=0.65)
    s.add_argument("--max-ratio", type=float, default=0.25)
    s.add_argument("--min-local-max-ratio", type=float, default=0.05)
    s.add_argument("--smooth-frac", type=float, default=0.1)
    s.add_argument("--k", type=int, default=21)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_kmer_screen)

    s = sub.add_parser("kmer-assess")
    s.add_argument("--screen", required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--samples", required=True)
    s.add_argument("--review-id", default="")
    s.add_argument("--block", action="store_true")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_kmer_assess)

    s = sub.add_parser("method-snapshot")
    s.add_argument("--method-file", required=True)
    s.add_argument("--science-version", type=int, required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_method_snapshot)

    s = sub.add_parser("import-assembly")
    for k in ("sample", "aid", "gfa", "fa", "log", "link-gfa", "link-fa", "link-log", "out"):
        s.add_argument(f"--{k}", required=True)
    s.add_argument("--register")
    s.set_defaults(func=cmd_import_assembly)

    s = sub.add_parser("validate-hifi")
    s.add_argument("--stats", nargs="+", required=True, help="INPUT_ID=stats.json")
    s.add_argument("--bam-count", nargs="*", help="INPUT_ID=count.txt")
    s.add_argument("--inventory", required=True)
    s.add_argument("--genome-size", type=int, required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_hifi)

    s = sub.add_parser("validate-ont")
    s.add_argument("--input-stats", nargs="+", required=True, help="INPUT_ID=stats.json")
    s.add_argument("--retained-stats", nargs="+", required=True)
    s.add_argument("--input-ids", nargs="+", required=True)
    s.add_argument("--retained-ids", nargs="+", required=True)
    s.add_argument("--min-len", type=int, required=True)
    s.add_argument("--min-retained", type=int, default=1)
    s.add_argument("--inventory", required=True)
    s.add_argument("--genome-size", type=int, required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_ont)

    s = sub.add_parser("validate-busco")
    s.add_argument("--json", required=True)
    s.add_argument("--dataset-cfg")
    s.add_argument("--lineage")
    s.add_argument("--expected-version")
    s.add_argument("--pin", help="JSON object of pinned dataset fields")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_busco)

    s = sub.add_parser("validate-merqury")
    s.add_argument("--qv", required=True)
    s.add_argument("--completeness", required=True)
    s.add_argument("--asm-name", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_merqury)

    s = sub.add_parser("validate-asm-stats")
    s.add_argument("--stats", required=True)
    s.add_argument("--check", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_asm)

    s = sub.add_parser("summarize")
    for k in ("sample", "method", "asm-p", "hifiasm-log", "out-json", "out-tsv", "validation-id"):
        s.add_argument(f"--{k}", required=True)
    s.add_argument("--science-version", type=int, required=True)
    for k in ("hifi", "kmer", "asm-a", "ont", "merqury", "method-record", "method-fingerprint",
              "imported"):
        s.add_argument(f"--{k}")
    s.add_argument("--busco", nargs="*")
    s.set_defaults(func=cmd_summarize)

    s = sub.add_parser("combine")
    s.add_argument("--inputs", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_combine)

    a = p.parse_args(argv)
    if a.cmd == "run" and a.command and a.command[0] == "--":
        a.command = a.command[1:]
    a.func(a)


if __name__ == "__main__":
    main()
