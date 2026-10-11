# Assembly QC: contiguity statistics, BUSCO, Merqury.


def qc_dir(s, m):
    return f"{RES}/{s}/qc/{m}"


def busco_lineages():
    return list(SCI["qc"]["busco"]["lineages"])


def busco_lineage_path(lineage):
    base = CFG["local"].get("busco_lineage_dir")
    if not base:
        raise WorkflowError("busco_lineage_dir is not set in config/local/local.yaml")
    return os.path.join(os.path.expanduser(base), lineage)


def busco_json(s, m, lineage):
    return f"{qc_dir(s, m)}/busco/{lineage}/short_summary.specific.{lineage}.{lineage}.json"


rule assembly_stats:
    input:
        fa=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.{{ctg}}.fa",
    output:
        stats=f"{RES}/{{sample}}/qc/{{method}}/stats/{{ctg}}.seqkit_stats.tsv",
        contigs=f"{RES}/{{sample}}/qc/{{method}}/stats/{{ctg}}.contigs.tsv",
    log:
        prov=f"{RES}/{{sample}}/provenance/assembly_stats/{{method}}.{{ctg}}.json",
    threads: R("assembly_stats")["threads"]
    resources:
        mem_mb=R("assembly_stats")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "assembly_stats",
            f"{tool('seqkit')} stats -a -T -j {threads} {q(input.fa)} > {q(output.stats)} && "
            f"{tool('seqkit')} fx2tab -n -i -l -g -B N {q(input.fa)} > {q(output.contigs)}",
            wc.sample, prov_path(wc.sample, "assembly_stats", f"{wc.method}.{wc.ctg}"),
            tools=["seqkit"], method_id=wc.method, inputs=[input.fa],
            outputs=[output.stats, output.contigs]),
    shell:
        "{params.cmd}"


rule small_contigs:
    """Contigs <= qc.small_contig_report_max_len, for organelle/contamination follow-up."""
    input:
        f"{RES}/{{sample}}/qc/{{method}}/stats/p_ctg.contigs.tsv",
    output:
        f"{RES}/{{sample}}/qc/{{method}}/stats/p_ctg.small_contigs.tsv",
    params:
        maxlen=SCI["qc"]["small_contig_report_max_len"],
    shell:
        "(printf 'contig\\tlength\\tgc_pct\\tn_pct\\n'; "
        "awk -F'\\t' -v m={params.maxlen} '$2 <= m' {input:q} | sort -t$'\\t' -k2,2n) > {output:q}"


rule validate_assembly_stats:
    input:
        stats=f"{RES}/{{sample}}/qc/{{method}}/stats/{{ctg}}.seqkit_stats.tsv",
        check=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.{{ctg}}.gfa_check.json",
    output:
        f"{RES}/{{sample}}/qc/{{method}}/stats/{{ctg}}.validated.json",
    shell:
        "{TOOL} validate-asm-stats --stats {input.stats:q} --check {input.check:q} --out {output:q}"


rule busco:
    input:
        fa=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.p_ctg.fa",
        check=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.p_ctg.gfa_check.json",
    output:
        f"{RES}/{{sample}}/qc/{{method}}/busco/{{lineage}}/short_summary.specific.{{lineage}}.{{lineage}}.json",
    log:
        prov=f"{RES}/{{sample}}/provenance/busco/{{method}}.{{lineage}}.json",
        stdout=f"{RES}/{{sample}}/qc/{{method}}/busco/{{lineage}}.busco.log",
    threads: R("busco")["threads"]
    resources:
        mem_mb=R("busco")["mem_mb"],
    params:
        cmd=lambda wc, input, threads: wrap(
            "busco",
            f"{tool('busco')} -i {q(input.fa)} -m {SCI['qc']['busco']['mode']} "
            f"-l {q(busco_lineage_path(wc.lineage))} -o {wc.lineage} "
            f"--out_path {q(qc_dir(wc.sample, wc.method) + '/busco')} -c {threads} --offline "
            f"--{SCI['qc']['busco']['gene_predictor']} --skip_bbtools -f "
            f"> {q(qc_dir(wc.sample, wc.method) + '/busco/' + wc.lineage + '.busco.log')} 2>&1",
            wc.sample, prov_path(wc.sample, "busco", f"{wc.method}.{wc.lineage}"),
            tools=["busco"], method_id=wc.method, inputs=[input.fa],
            outputs=[busco_json(wc.sample, wc.method, wc.lineage)]),
    shell:
        "mkdir -p $(dirname $(dirname {output:q})) && {params.cmd}"


rule validate_busco:
    input:
        f"{RES}/{{sample}}/qc/{{method}}/busco/{{lineage}}/short_summary.specific.{{lineage}}.{{lineage}}.json",
    output:
        f"{RES}/{{sample}}/qc/{{method}}/busco/{{lineage}}.validated.json",
    params:
        cfg=lambda wc: os.path.join(busco_lineage_path(wc.lineage), "dataset.cfg"),
    shell:
        "{TOOL} validate-busco --json {input:q} --dataset-cfg {params.cfg:q} --out {output:q}"


rule merqury:
    input:
        meryl=f"{RES}/{{sample}}/kmer/hifi.k{K}.meryl",
        fa=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.p_ctg.fa",
        check=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.p_ctg.gfa_check.json",
    output:
        qv=f"{RES}/{{sample}}/qc/{{method}}/merqury/{{sample}}.qv",
        comp=f"{RES}/{{sample}}/qc/{{method}}/merqury/{{sample}}.completeness.stats",
    log:
        prov=f"{RES}/{{sample}}/provenance/merqury/{{method}}.json",
        stdout=f"{RES}/{{sample}}/qc/{{method}}/merqury/merqury.log",
    threads: R("merqury")["threads"]
    resources:
        mem_mb=R("merqury")["mem_mb"],
    params:
        cmd=lambda wc, input, threads: wrap(
            "merqury",
            f"cd {q(qc_dir(wc.sample, wc.method) + '/merqury')} && "
            f"OMP_NUM_THREADS={threads} SLURM_CPUS_PER_TASK={threads} "
            f"{tool('merqury')} {q(input.meryl)} {q(input.fa)} {wc.sample} > merqury.log 2>&1",
            wc.sample, prov_path(wc.sample, "merqury", wc.method), tools=["merqury"],
            method_id=wc.method, inputs=[input.meryl, input.fa]),
    shell:
        "mkdir -p $(dirname {output.qv:q}) && {params.cmd}"


rule validate_merqury:
    input:
        qv=f"{RES}/{{sample}}/qc/{{method}}/merqury/{{sample}}.qv",
        comp=f"{RES}/{{sample}}/qc/{{method}}/merqury/{{sample}}.completeness.stats",
    output:
        f"{RES}/{{sample}}/qc/{{method}}/merqury/merqury.validated.json",
    params:
        asm=lambda wc: f"{wc.sample}.{wc.method}.p_ctg",
    shell:
        "{TOOL} validate-merqury --qv {input.qv:q} --completeness {input.comp:q} "
        "--asm-name {params.asm} --out {output:q}"
