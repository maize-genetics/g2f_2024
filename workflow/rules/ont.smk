# ONT preprocessing: QC of each supplied file, ultra-long length filtering,
# and the read list passed to hifiasm --ul.
#
# Inputs are basecalled FASTQ(.gz), BAM (primary records only), or an
# uncompressed tar of MinKNOW FASTQ (pass members streamed in place, fail
# members excluded; scripts/g2f/archive.py). Raw signal (POD5/FAST5) and ONT
# FASTA are rejected at discovery.
#
# No unnecessary copies: hifiasm takes a comma-separated --ul file list, so
# filtered files are never merged, and an acknowledged pre-filtered FASTQ is
# used as supplied (after verifying every read meets the method threshold)
# instead of being re-filtered into a duplicate.


def ont_stream(s, iid, threads):
    """Shell fragment that writes the ONT reads of one input as FASTQ to stdout."""
    i = get_input(s, "ont", iid)
    if discovery.is_released(i):
        return released_guard(s, "ont", iid, "")
    if i.format == "bam":
        return f"{tool('samtools')} fastq -F 0x900 -@ {threads} {q(i.path)} 2>/dev/null"
    if i.format == "tar":
        # Pass members streamed from the archive in place (g2f.archive); never extracted.
        return f"{TOOL} tar-fastq {q(i.path)} | {tool('seqkit')} seq -j {threads} -"
    return f"{tool('seqkit')} seq -j {threads} {q(i.path)}"


def ont_used_as_supplied(s, iid):
    i = get_input(s, "ont", iid)
    return discovery.is_prefiltered(i) and i.format == "fastq"


def ont_input_stats(s, iid):
    return f"{RES}/{s}/reads/ont/stats/{iid}.input.json"


def ont_input_ids(s, iid):
    return f"{RES}/{s}/reads/ont/stats/{iid}.input_ids.tsv.gz"


def ont_dir(s, minlen):
    return f"{RES}/{s}/reads/ont/ul{minlen}"


def ont_ul_reads(s, iid, minlen):
    if ont_used_as_supplied(s, iid):
        return get_input(s, "ont", iid).path
    return f"{ont_dir(s, minlen)}/{iid}.fastq.gz"


def ont_retained_stats(s, iid, minlen):
    if ont_used_as_supplied(s, iid):
        return ont_input_stats(s, iid)
    return f"{ont_dir(s, minlen)}/stats/{iid}.json"


def ont_retained_ids(s, iid, minlen):
    if ont_used_as_supplied(s, iid):
        return ont_input_ids(s, iid)
    return f"{ont_dir(s, minlen)}/stats/{iid}.ids.tsv.gz"


def all_ont_ul_reads(s, minlen):
    return [ont_ul_reads(s, i.input_id, minlen) for i in used(s, "ont")]


def ont_validated(s, minlen):
    return f"{ont_dir(s, minlen)}/ont_reads.validated.json"


def ont_tools(s):
    return ["seqkit"] + (["samtools"] if any(i.format == "bam" for i in used(s, "ont")) else [])


rule ont_input_stats:
    """QC of each ONT file as supplied (raw, or acknowledged pre-filtered)."""
    input:
        reads=lambda wc: raw_input(wc.sample, "ont", wc.iid),
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
    output:
        stats=f"{RES}/{{sample}}/reads/ont/stats/{{iid}}.input.json",
        ids=f"{RES}/{{sample}}/reads/ont/stats/{{iid}}.input_ids.tsv.gz",
    log:
        prov=f"{RES}/{{sample}}/provenance/ont_input_stats/{{iid}}.json",
    threads: R("read_stats")["threads"]
    resources:
        mem_mb=R("read_stats")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "ont_input_stats",
            f"{ont_stream(wc.sample, wc.iid, threads)} | {tool('seqkit')} fx2tab -n -l -q "
            f"| {TOOL} readstats --tech ont --has-qual --ids-out {q(output.ids)} "
            f"--source {q(input.reads)} --out {q(output.stats)}",
            wc.sample, prov_path(wc.sample, "ont_input_stats", wc.iid),
            tools=ont_tools(wc.sample), inputs=[input.reads], outputs=[output.stats, output.ids]),
    shell:
        "{params.cmd}"


rule ont_filter:
    """Keep reads >= minlen bp (method parameter ont.min_read_len)."""
    input:
        reads=lambda wc: raw_input(wc.sample, "ont", wc.iid),
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
    output:
        reads=f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/{{iid}}.fastq.gz",
    log:
        prov=f"{RES}/{{sample}}/provenance/ont_filter/ul{{minlen}}.{{iid}}.json",
    threads: R("ont_filter")["threads"]
    resources:
        mem_mb=R("ont_filter")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "ont_filter",
            f"{ont_stream(wc.sample, wc.iid, threads)} | {tool('seqkit')} seq -m {wc.minlen} "
            f"-j {threads} -o {q(output.reads)}",
            wc.sample, prov_path(wc.sample, "ont_filter", f"ul{wc.minlen}.{wc.iid}"),
            tools=ont_tools(wc.sample), inputs=[input.reads], outputs=[output.reads]),
    shell:
        "{params.cmd}"


rule ont_filtered_stats:
    input:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/{{iid}}.fastq.gz",
    output:
        stats=f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/stats/{{iid}}.json",
        ids=f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/stats/{{iid}}.ids.tsv.gz",
    threads: R("read_stats")["threads"]
    resources:
        mem_mb=R("read_stats")["mem_mb"],
    params:
        seqkit=tool("seqkit"),
    shell:
        "{params.seqkit} fx2tab -n -l -q -j {threads} {input:q} | {TOOL} readstats "
        "--tech ont --has-qual --ids-out {output.ids:q} --source {input:q} --out {output.stats:q}"


def _pairs(s, fn):
    return " ".join(q(f"{i.input_id}={fn(i.input_id)}") for i in used(s, "ont"))


rule validate_ont_reads:
    """Gate: every ONT input usable, yields reported, no duplicated reads."""
    input:
        sup=lambda wc: [ont_input_stats(wc.sample, i.input_id) for i in used(wc.sample, "ont")],
        sup_ids=lambda wc: [ont_input_ids(wc.sample, i.input_id) for i in used(wc.sample, "ont")],
        ret=lambda wc: [ont_retained_stats(wc.sample, i.input_id, wc.minlen)
                        for i in used(wc.sample, "ont")],
        ret_ids=lambda wc: [ont_retained_ids(wc.sample, i.input_id, wc.minlen)
                            for i in used(wc.sample, "ont")],
        md5=lambda wc: [checksum_path(wc.sample, "ont", i.input_id) for i in used(wc.sample, "ont")],
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/ont_reads.validated.json",
    params:
        sup=lambda wc: _pairs(wc.sample, lambda i: ont_input_stats(wc.sample, i)),
        sup_ids=lambda wc: _pairs(wc.sample, lambda i: ont_input_ids(wc.sample, i)),
        ret=lambda wc: _pairs(wc.sample, lambda i: ont_retained_stats(wc.sample, i, wc.minlen)),
        ret_ids=lambda wc: _pairs(wc.sample, lambda i: ont_retained_ids(wc.sample, i, wc.minlen)),
        min_retained=lambda wc: min_retained_for(wc.minlen),
    shell:
        "{TOOL} validate-ont --input-stats {params.sup} --retained-stats {params.ret} "
        "--input-ids {params.sup_ids} --retained-ids {params.ret_ids} "
        "--min-len {wildcards.minlen} --min-retained {params.min_retained} "
        "--inventory {input.inventory:q} --genome-size {GENOME_SIZE} --out {output:q}"


def min_retained_for(minlen):
    vals = {(m.get("ont") or {}).get("min_retained_reads", 1)
            for m in CFG["methods"].values()
            if m.get("ont") and str(m["ont"]["min_read_len"]) == str(minlen)}
    return max(vals) if vals else 1
