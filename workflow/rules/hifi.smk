# HiFi preprocessing: BAM -> FASTQ (pbtk), read QC, full-file identity checks.
#
# FASTQ/FASTA inputs are used in place (never copied). BAM inputs are
# converted once with pbtk bam2fastq. hifiasm accepts several read files, so
# no concatenated copy of the HiFi reads is made.


def hifi_reads(s, iid):
    i = get_input(s, "hifi", iid)
    if i.format == "bam":
        return f"{RES}/{s}/reads/hifi/converted/{iid}.fastq.gz"
    return i.path


def hifi_has_qual(s, iid):
    return get_input(s, "hifi", iid).format in ("bam", "fastq")


def all_hifi_reads(s):
    return [hifi_reads(s, i.input_id) for i in used(s, "hifi")]


def hifi_validated(s):
    return f"{RES}/{s}/reads/hifi/hifi_reads.validated.json"


def bam2fastq_cmd(wc, input, output, threads):
    """pbtk needs a .pbi index next to the BAM. To keep data/ untouched, work on
    a symlink in results/: reuse an existing .pbi, otherwise build one there.
    bam2fastq can exit 0 with an empty FASTQ when the index is unusable, so
    validate_hifi_reads compares read counts with the BAM record count."""
    if discovery.is_released(get_input(wc.sample, "hifi", wc.iid)):
        return released_guard(wc.sample, "hifi", wc.iid, "")
    bam = input.bam
    work = output[0][:-len(".fastq.gz")] + ".pbtk"
    link = f"{work}/in.bam"
    index = (f"ln -sfn {q(bam + '.pbi')} {q(link + '.pbi')}" if os.path.exists(bam + ".pbi")
             else f"{tool('pbindex')} -j {threads} {q(link)}")
    cmd = (f"mkdir -p {q(work)} && ln -sfn {q(bam)} {q(link)} && {index} && "
           f"{tool('bam2fastq')} -j {threads} -o {q(output[0][:-len('.fastq.gz')])} {q(link)}")
    return released_guard(wc.sample, "hifi", wc.iid, wrap(
        "hifi_bam_to_fastq", cmd, wc.sample, prov_path(wc.sample, "hifi_bam_to_fastq", wc.iid),
        tools=["bam2fastq", "pbindex"], inputs=[bam], outputs=[output[0]]))


rule hifi_bam_to_fastq:
    input:
        bam=lambda wc: raw_input(wc.sample, "hifi", wc.iid),
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
    output:
        # Prepared hifiasm input. Kept by default; temporary when
        # retain_prepared_hifi is false (regenerable from the BAM while it exists).
        ((lambda p: p) if RETAIN_HIFI else temp)(
            f"{RES}/{{sample}}/reads/hifi/converted/{{iid}}.fastq.gz"),
    log:
        prov=f"{RES}/{{sample}}/provenance/hifi_bam_to_fastq/{{iid}}.json",
    threads: R("hifi_bam_to_fastq")["threads"]
    resources:
        mem_mb=R("hifi_bam_to_fastq")["mem_mb"],
    params:
        cmd=bam2fastq_cmd,
    shell:
        "{params.cmd}"


rule bam_record_count:
    """Primary-record count of a HiFi BAM, to verify the conversion is complete."""
    input:
        bam=lambda wc: raw_input(wc.sample, "hifi", wc.iid),
    output:
        f"{RES}/{{sample}}/reads/hifi/converted/{{iid}}.bam_primary_records.txt",
    log:
        prov=f"{RES}/{{sample}}/provenance/bam_record_count/{{iid}}.json",
    threads: R("bam_record_count")["threads"]
    resources:
        mem_mb=R("bam_record_count")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: released_guard(wc.sample, "hifi", wc.iid, wrap(
            "bam_record_count",
            f"{tool('samtools')} view -c -F 0x900 -@ {threads} {q(input.bam)} > {q(output[0])}",
            wc.sample, prov_path(wc.sample, "bam_record_count", wc.iid),
            tools=["samtools"], inputs=[input.bam], outputs=[output[0]])),
    shell:
        "{params.cmd}"


rule hifi_read_stats:
    """One streaming pass: yield, length/quality distributions, movie census."""
    input:
        reads=lambda wc: hifi_reads(wc.sample, wc.iid),
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
    output:
        f"{RES}/{{sample}}/reads/hifi/stats/{{iid}}.json",
    log:
        prov=f"{RES}/{{sample}}/provenance/hifi_read_stats/{{iid}}.json",
    threads: R("read_stats")["threads"]
    resources:
        mem_mb=R("read_stats")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "hifi_read_stats",
            f"{tool('seqkit')} fx2tab -n -l {'-q ' if hifi_has_qual(wc.sample, wc.iid) else ''}"
            f"-j {threads} {q(input.reads)} | {TOOL} readstats --tech hifi "
            f"{'--has-qual ' if hifi_has_qual(wc.sample, wc.iid) else ''}"
            f"--source {q(input.reads)} --out {q(output[0])}",
            wc.sample, prov_path(wc.sample, "hifi_read_stats", wc.iid),
            tools=["seqkit"], inputs=[input.reads], outputs=[output[0]]),
    shell:
        "{params.cmd}"


rule validate_hifi_reads:
    """Gate: complete conversion, CCS read names, no movie shared between inputs."""
    input:
        stats=lambda wc: [f"{RES}/{wc.sample}/reads/hifi/stats/{i.input_id}.json"
                          for i in used(wc.sample, "hifi")],
        counts=lambda wc: [f"{RES}/{wc.sample}/reads/hifi/converted/{i.input_id}.bam_primary_records.txt"
                           for i in used(wc.sample, "hifi") if i.format == "bam"],
        md5=lambda wc: [checksum_path(wc.sample, "hifi", i.input_id) for i in used(wc.sample, "hifi")],
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/reads/hifi/hifi_reads.validated.json",
    params:
        stats=lambda wc: " ".join(
            q(f"{i.input_id}={RES}/{wc.sample}/reads/hifi/stats/{i.input_id}.json")
            for i in used(wc.sample, "hifi")),
        counts=lambda wc: " ".join(
            q(f"{i.input_id}={RES}/{wc.sample}/reads/hifi/converted/{i.input_id}.bam_primary_records.txt")
            for i in used(wc.sample, "hifi") if i.format == "bam"),
    shell:
        "{TOOL} validate-hifi --stats {params.stats} --bam-count {params.counts} "
        "--inventory {input.inventory:q} --genome-size {GENOME_SIZE} --out {output:q}"
