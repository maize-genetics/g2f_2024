"""Validation of read sets and QC outputs.

Each function returns a dict with "ok", "errors" and the verified metrics.
The CLI writes the dict and exits non-zero when ok is False, so Snakemake
never creates the downstream ".validated.json" marker for failed checks.
"""
import gzip
import json
import math
import os
import re


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def _yield(stats):
    keys = ("n_reads", "total_bases", "read_n50", "mean_len", "max_len",
            "mean_read_qual", "approx_coverage")
    return {k: stats.get(k) for k in keys if k in stats}


def _check_aggregator(iid, stats, expected, errors):
    got = stats.get("aggregator_id")
    if expected and got != expected:
        errors.append(f"{iid}: read statistics were produced by different read-statistics code "
                      f"({got} vs current {expected}); regenerate them, e.g. "
                      "snakemake --forcerun hifi_read_stats ont_input_stats ont_filtered_stats")


def validate_hifi(stats_paths, bam_count_paths, inventory_path, genome_size, aggregator=None):
    """stats_paths/bam_count_paths: {input_id: path}."""
    inv = _load(inventory_path)
    acks = {i["input_id"]: i.get("notes", []) for i in inv["inputs"] if i["tech"] == "hifi"}
    errors, per_input, movie_owner = [], {}, {}
    total_reads = total_bases = 0
    for iid in sorted(stats_paths):
        s = _load(stats_paths[iid])
        _check_aggregator(iid, s, aggregator, errors)
        per_input[iid] = _yield(s)
        if s["n_reads"] <= 0:
            errors.append(f"{iid}: no reads")
        kinds = s.get("name_kinds", {})
        if kinds.get("subread"):
            errors.append(f"{iid}: {kinds['subread']} subread-style names (not HiFi)")
        if kinds.get("ccs/fwd") or kinds.get("ccs/rev"):
            errors.append(f"{iid}: by-strand CCS reads present")
        if s.get("unparsed_names") and "ack:unverified_hifi" not in acks.get(iid, []):
            errors.append(f"{iid}: {s['unparsed_names']} read names are not PacBio CCS names")
        if s.get("duplicate_read_ids"):
            errors.append(f"{iid}: {s['duplicate_read_ids']} duplicated movie/ZMW read names")
        for movie in s.get("census", {}):
            if movie in movie_owner:
                errors.append(f"movie {movie} occurs in both {movie_owner[movie]} and {iid}: "
                              "inputs overlap (duplicate or derived data); exclude one in samples.yaml")
            movie_owner.setdefault(movie, iid)
        if iid in bam_count_paths:
            with open(bam_count_paths[iid]) as fh:
                n_bam = int(fh.read().split()[0])
            per_input[iid]["bam_primary_records"] = n_bam
            if n_bam != s["n_reads"]:
                errors.append(f"{iid}: converted FASTQ has {s['n_reads']} reads but BAM has "
                              f"{n_bam} primary records (conversion incomplete or silently failed)")
        total_reads += s["n_reads"]
        total_bases += s["total_bases"]
    return {
        "ok": not errors, "errors": errors, "tech": "hifi",
        "per_input": per_input,
        "movies": {m: movie_owner[m] for m in sorted(movie_owner)},
        "total_reads": total_reads, "total_bases": total_bases,
        "approx_coverage": round(total_bases / genome_size, 2),
        "genome_size_bp": genome_size,
    }


def _read_ids(path):
    """Read the (read_id, length) table written by `readstats --ids-out`."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.strip():
                rid, _, n = line.rstrip("\n").partition("\t")
                out.append((rid, int(n) if n else None))
    return out


def _n50(lengths):
    srt = sorted(lengths, reverse=True)
    half, acc = sum(srt) / 2, 0
    for n in srt:
        acc += n
        if acc >= half:
            return n
    return 0


def validate_ont(input_paths, retained_paths, input_id_paths, retained_id_paths, min_len,
                 min_retained, inventory_path, genome_size, aggregator=None):
    """Validate the ONT reads passed to hifiasm --ul.

    input_*: statistics/IDs of each ONT file as supplied. retained_*: of the
    reads actually used (the length-filtered file, or the supplied file itself
    for an acknowledged pre-filtered input). Yields of pre-filtered inputs are
    reported as "supplied", never as raw sequencing yield."""
    inv = _load(inventory_path)
    notes = {i["input_id"]: i.get("notes", []) for i in inv["inputs"] if i["tech"] == "ont"}
    errors, warnings, per_input, owner = [], [], {}, {}
    sup_reads = sup_bases = raw_reads = raw_bases = ret_reads = ret_bases = 0
    any_prefiltered, retained_lengths = False, []
    for iid in sorted(input_paths):
        r, f = _load(input_paths[iid]), _load(retained_paths[iid])
        _check_aggregator(iid, r, aggregator, errors)
        _check_aggregator(iid, f, aggregator, errors)
        pre = "ack:prefiltered_ont" in notes.get(iid, [])
        any_prefiltered |= pre
        per_input[iid] = {"yield_type": "prefiltered_supplied" if pre else "raw",
                          "supplied": _yield(r), "retained": _yield(f),
                          "refiltered": not pre}
        if r["n_reads"] <= 0:
            errors.append(f"{iid}: ONT file contains no reads")
        if f["n_reads"] == 0 and "ack:no_qualifying_reads" not in notes.get(iid, []):
            errors.append(f"{iid}: no reads >= {min_len} bp; human review required "
                          "(the workflow does not silently drop ONT inputs)")
        if f["n_reads"] and f["min_len"] < min_len:
            errors.append(f"{iid}: {'acknowledged pre-filtered file' if pre else 'filtered reads'} "
                          f"contain reads shorter than {min_len} bp (shortest {f['min_len']}); "
                          "it does not satisfy the method's ONT length threshold")
        if r.get("duplicate_read_ids"):
            errors.append(f"{iid}: {r['duplicate_read_ids']} duplicated read IDs within file")
        for rid, _ in _read_ids(input_id_paths[iid]):
            prev = owner.get(rid)
            if prev is not None and prev != iid:
                errors.append(f"read ID {rid} occurs in both {prev} and {iid}: ONT inputs overlap")
                break
            owner[rid] = iid
        retained_lengths += [n for _, n in _read_ids(retained_id_paths[iid]) if n is not None]
        sup_reads += r["n_reads"]; sup_bases += r["total_bases"]
        if not pre:
            raw_reads += r["n_reads"]; raw_bases += r["total_bases"]
        ret_reads += f["n_reads"]; ret_bases += f["total_bases"]
    if len(retained_lengths) != ret_reads or sum(retained_lengths) != ret_bases:
        errors.append("retained read-length table is inconsistent with retained statistics")
    if ret_reads < min_retained:
        errors.append(f"only {ret_reads} ONT reads retained (< {min_retained}); human review required")
    if any_prefiltered:
        warnings.append("one or more ONT inputs were supplied already length-filtered; their "
                        "yield is the supplied (derived) yield, and the original sequencing yield "
                        "is unknown")
    raw_known = not any_prefiltered
    return {
        "ok": not errors, "errors": errors[:50], "warnings": warnings, "tech": "ont",
        "min_read_len": min_len, "per_input": per_input,
        "input_type": ("raw" if raw_known else
                       "prefiltered" if all(p["yield_type"] != "raw" for p in per_input.values())
                       else "mixed"),
        "supplied_reads": sup_reads, "supplied_bases": sup_bases,
        "raw_reads": raw_reads if raw_known else None,
        "raw_bases": raw_bases if raw_known else None,
        "raw_approx_coverage": round(raw_bases / genome_size, 3) if raw_known else None,
        "retained_reads": ret_reads, "retained_bases": ret_bases,
        "retained_read_n50": _n50(retained_lengths) if retained_lengths else 0,
        "retained_approx_coverage": round(ret_bases / genome_size, 3),
        "genome_size_bp": genome_size,
    }


def validate_busco(json_path, lineage_cfg_path, expected_version=None, pin=None):
    """Check a BUSCO short-summary JSON for completeness, internal consistency,
    the BUSCO version, and the pinned lineage dataset identity."""
    errors = []
    d = _load(json_path)
    res = d.get("results") or {}
    params = d.get("parameters") or {}
    need = ["Complete percentage", "Single copy percentage", "Multi copy percentage",
            "Fragmented percentage", "Missing percentage", "n_markers", "Complete BUSCOs",
            "Single copy BUSCOs", "Multi copy BUSCOs", "Fragmented BUSCOs", "Missing BUSCOs"]
    for k in need:
        if k not in res:
            errors.append(f"BUSCO JSON missing results['{k}']")
    if not errors:
        counts = sum(res[k] for k in ("Single copy BUSCOs", "Multi copy BUSCOs",
                                      "Fragmented BUSCOs", "Missing BUSCOs"))
        if counts != res["n_markers"]:
            errors.append(f"BUSCO categories sum to {counts}, expected n_markers {res['n_markers']}")
        for k in need[:5]:
            if not 0 <= res[k] <= 100:
                errors.append(f"{k} out of range: {res[k]}")
    dataset = {}
    if lineage_cfg_path and os.path.exists(lineage_cfg_path):
        with open(lineage_cfg_path) as fh:
            dataset = dict(l.strip().split("=", 1) for l in fh if "=" in l)
    elif lineage_cfg_path:
        errors.append(f"lineage dataset.cfg not found: {lineage_cfg_path}")
    version = (d.get("versions") or {}).get("busco")
    if expected_version and version != str(expected_version):
        errors.append(f"BUSCO version {version} != required {expected_version}")
    lineage = (d.get("lineage_dataset") or {}).get("name")
    for key, want in (pin or {}).items():
        if want is None or key == "archive_md5":
            continue
        if key == "datasets_version":
            got = params.get("datasets_version")
        elif key == "name":
            got = lineage
        else:
            got = dataset.get(key)
        if str(got) != str(want):
            errors.append(f"lineage {key} = {got!r}, but the pinned dataset requires {want!r}")
    if pin is not None and dataset.get("name") and lineage and dataset["name"] != lineage:
        errors.append(f"BUSCO used lineage {lineage}, but dataset.cfg is for {dataset['name']}")
    if pin is not None and dataset.get("number_of_BUSCOs") and res.get("n_markers") is not None \
            and int(dataset["number_of_BUSCOs"]) != res["n_markers"]:
        errors.append("n_markers does not match the lineage dataset.cfg number_of_BUSCOs")
    return {
        "ok": not errors, "errors": errors,
        "lineage": lineage,
        "datasets_version": params.get("datasets_version"),
        "lineage_dataset_cfg": dataset,
        "pinned": pin,
        "versions": d.get("versions"),
        "mode": params.get("mode"),
        "gene_predictor": params.get("gene_predictor"),
        "one_line_summary": res.get("one_line_summary"),
        "complete_pct": res.get("Complete percentage"),
        "single_pct": res.get("Single copy percentage"),
        "duplicated_pct": res.get("Multi copy percentage"),
        "fragmented_pct": res.get("Fragmented percentage"),
        "missing_pct": res.get("Missing percentage"),
        "n_markers": res.get("n_markers"),
        "internal_stop_codon_count": res.get("internal_stop_codon_count"),
        "internal_stop_codon_percent": res.get("internal_stop_codon_percent"),
    }


def validate_merqury(qv_path, completeness_path, asm_name):
    """Parse Merqury <prefix>.qv and <prefix>.completeness.stats for one assembly."""
    errors, qv_row, comp_row = [], None, None
    with open(qv_path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    for r in rows:
        if r[0] == asm_name:
            qv_row = r
    if qv_row is None and len(rows) == 1:
        qv_row = rows[0]
    with open(completeness_path) as fh:
        crow = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    for r in crow:
        if r[0] == asm_name:
            comp_row = r
    if comp_row is None and len(crow) == 1:
        comp_row = crow[0]
    out = {"asm_name": asm_name}
    try:
        out.update(kmers_asm_only=int(qv_row[1]), kmers_total=int(qv_row[2]),
                   qv=float(qv_row[3]), error_rate=float(qv_row[4]))
        if not (0 < out["qv"] < 100 or math.isinf(out["qv"])):
            errors.append(f"QV out of range: {out['qv']}")
    except (TypeError, IndexError, ValueError):
        errors.append(f"cannot parse QV row for {asm_name} in {qv_path}")
    try:
        out.update(solid_kmers_in_asm=int(comp_row[2]), solid_kmers_total=int(comp_row[3]),
                   kmer_completeness_pct=float(comp_row[4]))
        if not 0 <= out["kmer_completeness_pct"] <= 100:
            errors.append(f"completeness out of range: {out['kmer_completeness_pct']}")
    except (TypeError, IndexError, ValueError):
        errors.append(f"cannot parse completeness row for {asm_name} in {completeness_path}")
    out.update(ok=not errors, errors=errors)
    return out


def parse_seqkit_stats(path):
    """Parse `seqkit stats -a` output (one row per file), tabular (-T) or aligned."""
    with open(path) as fh:
        first, second = fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
    if "\t" in first:
        header, row = first.split("\t"), second.split("\t")
    else:
        header, row = first.split(), second.split()
    d = dict(zip(header, row))

    def num(k, f=int):
        v = d.get(k)
        return f(v.replace(",", "")) if v not in (None, "") else None
    return {
        "n_contigs": num("num_seqs"), "total_bp": num("sum_len"), "min_len": num("min_len"),
        "max_len": num("max_len"), "n50": num("N50"), "l50": num("N50_num"),
        "gc_pct": num("GC(%)", float), "gap_bp": num("sum_gap"), "n_bp": num("sum_n"),
    }


def validate_assembly_stats(stats_path, gfa_check_path):
    st = parse_seqkit_stats(stats_path)
    chk = _load(gfa_check_path)
    errors = []
    if not chk.get("ok"):
        errors.append("GFA/FASTA consistency check did not pass")
    if st["n_contigs"] != chk.get("n_fasta_records"):
        errors.append(f"seqkit contig count {st['n_contigs']} != {chk.get('n_fasta_records')}")
    if st["total_bp"] != chk.get("fasta_total_bp"):
        errors.append(f"seqkit total {st['total_bp']} != {chk.get('fasta_total_bp')}")
    if st["n_contigs"] and st["gc_pct"] is not None and not 0 < st["gc_pct"] < 100:
        errors.append(f"GC% out of range: {st['gc_pct']}")
    if st["n_contigs"] == 0:
        st["empty"] = True     # e.g. hifiasm -l0 a_ctg with no alternate sequence
    return {"ok": not errors, "errors": errors, **st}


def parse_hifiasm_log(path):
    """Extract version, command and coverage estimates reported by hifiasm."""
    out = {"peaks": [], "purge_dups": []}
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("[M::main] Version:"):
                out["version"] = line.split(":", 3)[-1].strip()
            elif line.startswith("[M::main] CMD:"):
                out["command"] = line.split("CMD:", 1)[1].strip()
            elif line.startswith("[M::main] Real time:"):
                m = re.search(r"Real time: ([\d.]+) sec; CPU: ([\d.]+) sec; Peak RSS: ([\d.]+) GB", line)
                if m:
                    out.update(real_s=float(m.group(1)), cpu_s=float(m.group(2)),
                               peak_rss_gb=float(m.group(3)))
            elif "peak_hom:" in line:
                m = re.search(r"\[M::(\w+)\] peak_hom: (-?\d+); peak_het: (-?\d+)", line)
                if m:
                    out["peaks"].append({"stage": m.group(1), "peak_hom": int(m.group(2)),
                                         "peak_het": int(m.group(3))})
            elif line.startswith("[M::purge_dups]"):
                out["purge_dups"].append(line.split("]", 1)[1].strip())
    return out
