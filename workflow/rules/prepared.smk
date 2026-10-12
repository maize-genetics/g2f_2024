# Verification of the prepared hifiasm inputs (gate before hifiasm).
#
# Order per genotype: raw inputs validated (discovery, checksums) -> HiFi/ONT
# preprocessing and read QC -> each prepared file checked in one full pass
# (gzip integrity, MD5, record count) -> verify_prepared_inputs reconciles
# runs, read counts, ONT length filtering and the manifest -> hifiasm.
# Raw-copy deletion eligibility is reported separately and never acted on:
#   python3 scripts/g2f_tool.py raw-release   (see docs/full_scale_run_plan.md)


def prepared_check_path(s, tech, iid):
    if tech == "hifi":
        return f"{RES}/{s}/reads/hifi/checks/{iid}.json"
    return f"{ont_dir(s, method(METHOD[s])['ont']['min_read_len'])}/checks/{iid}.json"


def prepared_checks(s):
    """(tech, input_id, check path) for every prepared file that is not the raw file itself."""
    out = [("hifi", i.input_id, prepared_check_path(s, "hifi", i.input_id))
           for i in used(s, "hifi") if i.format == "bam"]
    if uses_ont(METHOD[s]):
        out += [("ont", i.input_id, prepared_check_path(s, "ont", i.input_id))
                for i in used(s, "ont") if not ont_used_as_supplied(s, i.input_id)]
    return out


def prepared_verified(s):
    return f"{RES}/{s}/inputs/prepared_inputs.verified.json"


def hifiasm_input_files(s):
    m = METHOD[s]
    files = all_hifi_reads(s)
    if uses_ont(m):
        files += all_ont_ul_reads(s, method(m)["ont"]["min_read_len"])
    return files


rule check_prepared_hifi:
    """One full pass over a converted HiFi FASTQ: gzip integrity, MD5, records."""
    input:
        f"{RES}/{{sample}}/reads/hifi/converted/{{iid}}.fastq.gz",
    output:
        f"{RES}/{{sample}}/reads/hifi/checks/{{iid}}.json",
    resources:
        mem_mb=2000,
    shell:
        "{TOOL} check-gz --path {input:q} --out {output:q}"


rule check_prepared_ont:
    """One full pass over a filtered ONT FASTQ: gzip integrity, MD5, records."""
    input:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/{{iid}}.fastq.gz",
    output:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/checks/{{iid}}.json",
    resources:
        mem_mb=2000,
    shell:
        "{TOOL} check-gz --path {input:q} --out {output:q}"


def _verify_inputs(wc):
    s = wc.sample
    d = {"checks": [p for _, _, p in prepared_checks(s)],
         "hifi_ok": hifi_validated(s),
         "md5": all_checksums(s),
         "inventory": ancient(inventory_path(s))}
    if uses_ont(METHOD[s]):
        d["ont_ok"] = ont_validated(s, method(METHOD[s])["ont"]["min_read_len"])
    return d


rule verify_prepared_inputs:
    """Gate: every prepared hifiasm input verified and every run accounted for."""
    input:
        unpack(_verify_inputs),
    output:
        f"{RES}/{{sample}}/inputs/prepared_inputs.verified.json",
    params:
        checks=lambda wc: " ".join(q(f"{t}:{i}={p}") for t, i, p in prepared_checks(wc.sample)),
        md5=lambda wc: " ".join(q(f"{i.tech}:{i.input_id}={checksum_path(wc.sample, i.tech, i.input_id)}")
                                for i in used(wc.sample, "hifi") + used(wc.sample, "ont")),
        ont=lambda wc: (f"--ont-validated {q(ont_validated(wc.sample, method(METHOD[wc.sample])['ont']['min_read_len']))} "
                        f"--min-len {method(METHOD[wc.sample])['ont']['min_read_len']}"
                        if uses_ont(METHOD[wc.sample]) else ""),
        manifest=lambda wc: f"--manifest {q(CFG['raw_manifest'])}" if CFG.get("raw_manifest") else "",
        hifiasm=lambda wc: " ".join(q(p) for p in hifiasm_input_files(wc.sample)),
    shell:
        "{TOOL} verify-prepared --sample {wildcards.sample} --inventory {input.inventory:q} "
        "--hifi-validated {input.hifi_ok:q} {params.ont} --checks {params.checks} "
        "--raw-md5 {params.md5} {params.manifest} --hifiasm-inputs {params.hifiasm} "
        "--out {output:q}"
