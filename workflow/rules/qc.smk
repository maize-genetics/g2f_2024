# Assembly QC: GFA/FASTA consistency, contiguity statistics, BUSCO, Merqury.
#
# {akey} is either a pipeline method_id or "imported-<id>" for a registered
# historical assembly. Pipeline assemblies have p_ctg and a_ctg; imported
# ones are checked for p_ctg only.


def is_imported(akey):
    return akey.startswith("imported-")


def qc_dir(s, akey):
    return f"{RES}/{s}/qc/{akey}"


def assembly_fa(s, akey, ctg):
    if is_imported(akey):
        if ctg != "p_ctg":
            raise WorkflowError(f"{akey}: only p_ctg is registered for imported assemblies")
        return imported_fasta(s, imported(s, akey)["aid"])
    return asm_fasta(s, akey, ctg)


def assembly_gfa(s, akey, ctg):
    if is_imported(akey):
        return imported_gfa(s, imported(s, akey)["aid"])
    return asm_gfa(s, akey, ctg)


def assembly_log(s, akey):
    if is_imported(akey):
        return imported_log(s, imported(s, akey)["aid"])
    return hifiasm_log(s, akey)


def busco_lineages():
    return list(SCI["qc"]["busco"]["lineages"])


def busco_pin(lineage):
    return dict(SCI["qc"]["busco"]["lineages"][lineage] or {})


def busco_lineage_path(lineage):
    base = CFG["local"].get("busco_lineage_dir")
    if not base:
        raise WorkflowError("busco_lineage_dir is not set in config/local/local.yaml")
    return os.path.join(os.path.expanduser(base), lineage)


def busco_json(s, akey, lineage):
    return f"{qc_dir(s, akey)}/busco/{lineage}/short_summary.specific.{lineage}.{lineage}.json"


def gfa_check(s, akey, ctg):
    return f"{qc_dir(s, akey)}/stats/{ctg}.gfa_check.json"


rule check_gfa_fasta:
    """Same names, order, lengths and sequence hashes in GFA and FASTA; ACGTN only."""
    input:
        gfa=lambda wc: assembly_gfa(wc.sample, wc.akey, wc.ctg),
        fa=lambda wc: assembly_fa(wc.sample, wc.akey, wc.ctg),
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.gfa_check.json",
    resources:
        mem_mb=R("validate_fasta")["mem_mb"],
    shell:
        "{TOOL} check-gfa --gfa {input.gfa:q} --fasta {input.fa:q} --out {output:q}"


rule assembly_stats:
    input:
        fa=lambda wc: assembly_fa(wc.sample, wc.akey, wc.ctg),
    output:
        stats=f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.seqkit_stats.tsv",
        contigs=f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.contigs.tsv",
    log:
        prov=f"{RES}/{{sample}}/provenance/assembly_stats/{{akey}}.{{ctg}}.json",
    threads: R("assembly_stats")["threads"]
    resources:
        mem_mb=R("assembly_stats")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads: wrap(
            "assembly_stats",
            f"{tool('seqkit')} stats -a -T -j {threads} {q(input.fa)} > {q(output.stats)} && "
            f"{tool('seqkit')} fx2tab -n -i -l -g -B N {q(input.fa)} > {q(output.contigs)}",
            wc.sample, prov_path(wc.sample, "assembly_stats", f"{wc.akey}.{wc.ctg}"),
            tools=["seqkit"], method_id=wc.akey, inputs=[input.fa],
            outputs=[output.stats, output.contigs]),
    shell:
        "{params.cmd}"


rule small_contigs:
    """Contigs <= qc.small_contig_report_max_len, for organelle/contamination follow-up."""
    input:
        contigs=f"{RES}/{{sample}}/qc/{{akey}}/stats/p_ctg.contigs.tsv",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/stats/p_ctg.small_contigs.tsv",
    params:
        maxlen=SCI["qc"]["small_contig_report_max_len"],
    shell:
        "(printf 'contig\\tlength\\tgc_pct\\tn_pct\\n'; "
        "awk -F'\\t' -v m={params.maxlen} '$2 <= m' {input.contigs:q} | sort -t$'\\t' -k2,2n) "
        "> {output:q}"


rule validate_assembly_stats:
    input:
        stats=f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.seqkit_stats.tsv",
        check=f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.gfa_check.json",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/stats/{{ctg}}.validated.json",
    shell:
        "{TOOL} validate-asm-stats --stats {input.stats:q} --check {input.check:q} --out {output:q}"


rule busco:
    input:
        fa=lambda wc: assembly_fa(wc.sample, wc.akey, "p_ctg"),
        check=lambda wc: ancient(gfa_check(wc.sample, wc.akey, "p_ctg")),
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/busco/{{lineage}}/short_summary.specific.{{lineage}}.{{lineage}}.json",
    log:
        prov=f"{RES}/{{sample}}/provenance/busco/{{akey}}.{{lineage}}.json",
        stdout=f"{RES}/{{sample}}/qc/{{akey}}/busco/{{lineage}}.busco.log",
    threads: R("busco")["threads"]
    resources:
        mem_mb=R("busco")["mem_mb"],
    params:
        cmd=lambda wc, input, threads: wrap(
            "busco",
            f"{tool('busco')} -i {q(input.fa)} -m {SCI['qc']['busco']['mode']} "
            f"-l {q(busco_lineage_path(wc.lineage))} -o {wc.lineage} "
            f"--out_path {q(qc_dir(wc.sample, wc.akey) + '/busco')} -c {threads} --offline "
            f"--{SCI['qc']['busco']['gene_predictor']} --skip_bbtools --opt-out-run-stats -f "
            f"> {q(qc_dir(wc.sample, wc.akey) + '/busco/' + wc.lineage + '.busco.log')} 2>&1",
            wc.sample, prov_path(wc.sample, "busco", f"{wc.akey}.{wc.lineage}"),
            tools=["busco"], method_id=wc.akey, inputs=[input.fa],
            outputs=[busco_json(wc.sample, wc.akey, wc.lineage)]),
    shell:
        "mkdir -p $(dirname $(dirname {output:q})) && {params.cmd}"


rule validate_busco:
    input:
        json=f"{RES}/{{sample}}/qc/{{akey}}/busco/{{lineage}}/short_summary.specific.{{lineage}}.{{lineage}}.json",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/busco/{{lineage}}.validated.json",
    params:
        cfg=lambda wc: os.path.join(busco_lineage_path(wc.lineage), "dataset.cfg"),
        pin=lambda wc: q(json.dumps(busco_pin(wc.lineage))),
        version=SCI["qc"]["busco"].get("expected_version", ""),
    shell:
        "{TOOL} validate-busco --json {input.json:q} --dataset-cfg {params.cfg:q} "
        "--lineage {wildcards.lineage} --expected-version {params.version} --pin {params.pin} "
        "--out {output:q}"


rule merqury:
    input:
        meryl=f"{RES}/{{sample}}/kmer/hifi.k{K}.meryl",
        fa=lambda wc: assembly_fa(wc.sample, wc.akey, "p_ctg"),
        check=lambda wc: ancient(gfa_check(wc.sample, wc.akey, "p_ctg")),
    output:
        qv=f"{RES}/{{sample}}/qc/{{akey}}/merqury/{{sample}}.qv",
        comp=f"{RES}/{{sample}}/qc/{{akey}}/merqury/{{sample}}.completeness.stats",
    log:
        prov=f"{RES}/{{sample}}/provenance/merqury/{{akey}}.json",
        stdout=f"{RES}/{{sample}}/qc/{{akey}}/merqury/merqury.log",
    threads: R("merqury")["threads"]
    resources:
        mem_mb=R("merqury")["mem_mb"],
    params:
        cmd=lambda wc, input, threads: wrap(
            "merqury",
            f"cd {q(qc_dir(wc.sample, wc.akey) + '/merqury')} && "
            f"OMP_NUM_THREADS={threads} "
            f"{tool('merqury')} {q(input.meryl)} {q(input.fa)} {wc.sample} > merqury.log 2>&1",
            wc.sample, prov_path(wc.sample, "merqury", wc.akey), tools=["merqury", "meryl"],
            method_id=wc.akey, inputs=[input.meryl, input.fa]),
    shell:
        "mkdir -p $(dirname {output.qv:q}) && {params.cmd}"


rule validate_merqury:
    input:
        qv=f"{RES}/{{sample}}/qc/{{akey}}/merqury/{{sample}}.qv",
        comp=f"{RES}/{{sample}}/qc/{{akey}}/merqury/{{sample}}.completeness.stats",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/qc/{{akey}}/merqury/merqury.validated.json",
    params:
        asm=lambda wc: os.path.basename(assembly_fa(wc.sample, wc.akey, "p_ctg"))[:-len(".fa")],
    shell:
        "{TOOL} validate-merqury --qv {input.qv:q} --completeness {input.comp:q} "
        "--asm-name {params.asm} --out {output:q}"


for _lin in busco_lineages():
    _base = CFG["local"].get("busco_lineage_dir")
    if _base and not os.path.isfile(os.path.join(os.path.expanduser(_base), _lin, "dataset.cfg")):
        sys.stderr.write(
            f"WARNING: BUSCO lineage {_lin} is not installed under {_base}; BUSCO jobs will fail "
            "until it is downloaded (see README, 'One-time setup').\n")
