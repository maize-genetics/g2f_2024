# ONT preprocessing: raw read QC, ultra-long length filtering, merge for --ul.
#
# Inputs are basecalled FASTQ(.gz) or unaligned/aligned BAM (primary records
# only); ONT FASTA is rejected at discovery because qualities are required. Raw signal (POD5/FAST5) is rejected at discovery.


def ont_stream(s, iid, threads):
    """Shell fragment that writes the ONT reads of one input as FASTQ to stdout."""
    i = get_input(s, "ont", iid)
    if i.format == "bam":
        return f"{tool('samtools')} fastq -F 0x900 -@ {threads} {q(i.path)} 2>/dev/null"
    return f"{tool('seqkit')} seq -j {threads} {q(i.path)}"


def ont_dir(s, minlen):
    return f"{RES}/{s}/reads/ont/ul{minlen}"


def ont_merged(s, minlen):
    return f"{ont_dir(s, minlen)}/{s}.ont_ul{minlen}.fastq.gz"


def ont_validated(s, minlen):
    return f"{ont_dir(s, minlen)}/ont_reads.validated.json"


def ont_tools(s):
    return ["seqkit"] + (["samtools"] if any(i.format == "bam" for i in used(s, "ont")) else [])


rule ont_raw_stats:
    input:
        reads=lambda wc: get_input(wc.sample, "ont", wc.iid).path,
        inventory=f"{RES}/{{sample}}/inputs/inventory.json",
    output:
        stats=f"{RES}/{{sample}}/reads/ont/stats/{{iid}}.raw.json",
        ids=f"{RES}/{{sample}}/reads/ont/stats/{{iid}}.raw_ids.txt.gz",
    log:
        prov=f"{RES}/{{sample}}/provenance/ont_raw_stats/{{iid}}.json",
    threads: R("read_stats")["threads"]
    resources:
        mem_mb=R("read_stats")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "ont_raw_stats",
            f"{ont_stream(wc.sample, wc.iid, threads)} | {tool('seqkit')} fx2tab -n -l -q "
            f"| {TOOL} readstats --tech ont --has-qual --genome-size {GENOME_SIZE} --ids-out {q(output.ids)} "
            f"--source {q(input.reads)} --out {q(output.stats)}",
            wc.sample, prov_path(wc.sample, "ont_raw_stats", wc.iid),
            tools=ont_tools(wc.sample), inputs=[input.reads], outputs=[output.stats]),
    shell:
        "{params.cmd}"


rule ont_filter:
    """Keep reads >= minlen bp (method parameter ont.min_read_len)."""
    input:
        reads=lambda wc: get_input(wc.sample, "ont", wc.iid).path,
        inventory=f"{RES}/{{sample}}/inputs/inventory.json",
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
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/stats/{{iid}}.json",
    threads: R("read_stats")["threads"]
    resources:
        mem_mb=R("read_stats")["mem_mb"],
    params:
        seqkit=tool("seqkit"),
    shell:
        "{params.seqkit} fx2tab -n -l -q -j {threads} {input:q} | {TOOL} readstats "
        "--tech ont --has-qual --genome-size {GENOME_SIZE} --source {input:q} --out {output:q}"


rule ont_merge:
    """Concatenate filtered ONT inputs (gzip members) into the single --ul file."""
    input:
        lambda wc: [f"{ont_dir(wc.sample, wc.minlen)}/{i.input_id}.fastq.gz"
                    for i in used(wc.sample, "ont")],
    output:
        reads=f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/{{sample}}.ont_ul{{minlen}}.fastq.gz",
    log:
        prov=f"{RES}/{{sample}}/provenance/ont_merge/ul{{minlen}}.json",
    params:
        cmd=lambda wc, input, output: wrap(
            "ont_merge", f"cat {' '.join(q(p) for p in input)} > {q(output.reads)}",
            wc.sample, prov_path(wc.sample, "ont_merge", f"ul{wc.minlen}"),
            inputs=list(input), outputs=[output.reads]),
    shell:
        "{params.cmd}"


use rule ont_filtered_stats as ont_merged_stats with:
    input:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/{{sample}}.ont_ul{{minlen}}.fastq.gz",
    output:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/stats/merged.json",


rule validate_ont_reads:
    """Gate: every ONT input usable, yields reported, no duplicated reads."""
    input:
        raw=lambda wc: [f"{RES}/{wc.sample}/reads/ont/stats/{i.input_id}.raw.json"
                        for i in used(wc.sample, "ont")],
        ids=lambda wc: [f"{RES}/{wc.sample}/reads/ont/stats/{i.input_id}.raw_ids.txt.gz"
                        for i in used(wc.sample, "ont")],
        filt=lambda wc: [f"{ont_dir(wc.sample, wc.minlen)}/stats/{i.input_id}.json"
                         for i in used(wc.sample, "ont")],
        merged=f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/stats/merged.json",
        md5=lambda wc: [checksum_path(wc.sample, "ont", i.input_id) for i in used(wc.sample, "ont")],
        inventory=f"{RES}/{{sample}}/inputs/inventory.json",
    output:
        f"{RES}/{{sample}}/reads/ont/ul{{minlen}}/ont_reads.validated.json",
    params:
        raw=lambda wc: " ".join(q(f"{i.input_id}={RES}/{wc.sample}/reads/ont/stats/{i.input_id}.raw.json")
                                for i in used(wc.sample, "ont")),
        ids=lambda wc: " ".join(q(f"{i.input_id}={RES}/{wc.sample}/reads/ont/stats/{i.input_id}.raw_ids.txt.gz")
                                for i in used(wc.sample, "ont")),
        filt=lambda wc: " ".join(q(f"{i.input_id}={ont_dir(wc.sample, wc.minlen)}/stats/{i.input_id}.json")
                                 for i in used(wc.sample, "ont")),
        min_retained=lambda wc: min_retained_for(wc.minlen),
    shell:
        "{TOOL} validate-ont --raw {params.raw} --filtered {params.filt} --ids {params.ids} "
        "--merged {input.merged:q} --min-len {wildcards.minlen} "
        "--min-retained {params.min_retained} --inventory {input.inventory:q} "
        "--genome-size {GENOME_SIZE} --out {output:q}"


def min_retained_for(minlen):
    vals = {(m.get("ont") or {}).get("min_retained_reads", 1)
            for m in CFG["methods"].values()
            if m.get("ont") and str(m["ont"]["min_read_len"]) == str(minlen)}
    return max(vals) if vals else 1
