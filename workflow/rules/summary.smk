# Per-assembly validated summaries and the project tables. Only validated
# JSON records feed these tables.
#
# Names are content-addressed: <VID> is the validation identity (science
# settings, method definitions, validation code) and <SID> a hash of the
# sample's samples.yaml entry (human decisions such as kmer_review). Under
# Snakemake's mtime-only rerun policy this is what makes a changed decision
# or validation logic produce a new, current summary (cheap jobs only),
# while summaries made under older settings keep their old names and can
# never be mistaken for current ones.
import hashlib as _hl


def sample_decision_id(s):
    entry = CFG["samples"].get(s) or {}
    return _hl.sha256(json.dumps(entry, sort_keys=True, default=str).encode()).hexdigest()[:6]


def _set_id(samples):
    return _hl.sha256(",".join(f"{s}:{sample_decision_id(s)}" for s in sorted(samples))
                      .encode()).hexdigest()[:6]


PROJECT_SUMMARY = f"{RES}/summary/assembly_summary.{VID}.{_set_id(SAMPLES)}.tsv"
IMPORTED_SUMMARY = (f"{RES}/summary/imported_assemblies_summary.{VID}."
                    f"{_set_id({s for s, _ in IMPORTED})}.tsv")

localrules: sample_summary, project_summary, imported_summary


def validated_marker(s, akey):
    return f"{RES}/{s}/summary/{s}.{akey}.validated.{VID}.{sample_decision_id(s)}.json"


def summary_tsv(s, akey):
    return f"{RES}/{s}/summary/{s}.{akey}.summary.{VID}.{sample_decision_id(s)}.tsv"


def reads_qc_targets(s):
    out = [inventory_path(s), hifi_validated(s), kmer_assessment_path(s)]
    if uses_ont(METHOD[s]):
        out.append(ont_validated(s, method(METHOD[s])["ont"]["min_read_len"]))
    return out


def summary_inputs(wc):
    s, akey = wc.sample, wc.akey
    if wc.sid != sample_decision_id(s):
        raise WorkflowError(f"summary for {s} requested with decision id {wc.sid}, but the current "
                            f"samples.yaml entry has id {sample_decision_id(s)}")
    d = {"stamp": STAMP, "asm_p": f"{qc_dir(s, akey)}/stats/p_ctg.validated.json",
         "small": f"{qc_dir(s, akey)}/stats/p_ctg.small_contigs.tsv",
         "hifiasm_log": assembly_log(s, akey)}
    if is_imported(akey):
        rec = imported(s, akey)
        d["imported"] = import_record_path(s, rec["aid"])
        qcs = rec["qc"]
        if "busco" in qcs:
            d["busco"] = [f"{qc_dir(s, akey)}/busco/{l}.validated.json" for l in busco_lineages()]
        if "merqury" in qcs:
            d["merqury"] = f"{qc_dir(s, akey)}/merqury/merqury.validated.json"
        return d
    m = akey
    d.update(hifi=hifi_validated(s), kmer=kmer_assessment_path(s),
             asm_a=f"{qc_dir(s, m)}/stats/a_ctg.validated.json",
             method_record=method_record_path(s, m),
             busco=[f"{qc_dir(s, m)}/busco/{l}.validated.json" for l in busco_lineages()])
    if SCI["qc"]["merqury"].get("enabled", True):
        d["merqury"] = f"{qc_dir(s, m)}/merqury/merqury.validated.json"
    if uses_ont(m):
        d["ont"] = ont_validated(s, method(m)["ont"]["min_read_len"])
    return d


def _opt(input, key, flag):
    return f"--{flag} {q(input[key])}" if key in input.keys() else ""


rule sample_summary:
    input:
        unpack(summary_inputs),
    output:
        tsv=f"{RES}/{{sample}}/summary/{{sample}}.{{akey}}.summary.{VID}.{{sid}}.tsv",
        marker=f"{RES}/{{sample}}/summary/{{sample}}.{{akey}}.validated.{VID}.{{sid}}.json",
    wildcard_constraints:
        sid=r"[0-9a-f]{6}",
    params:
        opts=lambda wc, input: " ".join(_opt(input, k, f) for k, f in (
            ("hifi", "hifi"), ("ont", "ont"), ("kmer", "kmer"), ("asm_a", "asm-a"),
            ("merqury", "merqury"), ("method_record", "method-record"), ("imported", "imported"))),
        busco=lambda wc, input: ("--busco " + " ".join(q(b) for b in input.busco)
                                 if "busco" in input.keys() else ""),
        fp=lambda wc: "" if is_imported(wc.akey) else f"--method-fingerprint {METHOD_FP[wc.akey]}",
    shell:
        "{TOOL} summarize --sample {wildcards.sample} --method {wildcards.akey} "
        "--science-version {SCI[science_version]} --validation-id {VID} "
        "--asm-p {input.asm_p:q} {params.opts} {params.busco} {params.fp} "
        "--hifiasm-log {input.hifiasm_log:q} --out-tsv {output.tsv:q} --out-json {output.marker:q}"


rule project_summary:
    """Assemblies produced by this workflow only (one per sample)."""
    input:
        lambda wc: [summary_tsv(s, METHOD[s]) for s in SAMPLES],
    output:
        PROJECT_SUMMARY,
    shell:
        "{TOOL} combine --inputs {input:q} --out {output:q}"


def imported_targets():
    return [validated_marker(s, k) for (s, k), r in sorted(IMPORTED.items()) if r["qc"]]


rule imported_summary:
    """Registered historical assemblies only; never mixed with pipeline output."""
    input:
        lambda wc: [summary_tsv(s, k) for (s, k), r in sorted(IMPORTED.items()) if r["qc"]],
    output:
        IMPORTED_SUMMARY,
    shell:
        "{TOOL} combine --inputs {input:q} --out {output:q}"
