# Read-only registration of historical assemblies (config/external_assemblies.yaml)
# for QC validation and benchmarking. Nothing is written to the original
# directory: files are linked into results/<SAMPLE>/imported/<id>/, so any
# tool that writes next to its input (index files, Merqury links) writes in
# results/. Their QC is reported as assembly_source "imported_legacy" and
# never mixed with assemblies produced by this workflow.

localrules: import_assembly


def imported(s, akey):
    return IMPORTED[(s, akey)]


def imported_dir(s, aid):
    return f"{RES}/{s}/imported/{aid}"


def imported_fasta(s, aid):
    return f"{imported_dir(s, aid)}/{s}.{aid}.p_ctg.fa"


def imported_gfa(s, aid):
    return f"{imported_dir(s, aid)}/{s}.{aid}.p_ctg.gfa"


def imported_log(s, aid):
    return f"{imported_dir(s, aid)}/{s}.{aid}.hifiasm.stderr.log"


rule import_assembly:
    """Record identity (size, mtime, fingerprint) and link the historical files."""
    input:
        gfa=lambda wc: IMPORTED[(wc.sample, f"imported-{wc.aid}")]["p_ctg_gfa"],
        fa=lambda wc: IMPORTED[(wc.sample, f"imported-{wc.aid}")]["p_ctg_fa"],
        log=lambda wc: IMPORTED[(wc.sample, f"imported-{wc.aid}")]["hifiasm_log"],
    output:
        record=f"{RES}/{{sample}}/imported/{{aid}}/import.json",
        gfa=f"{RES}/{{sample}}/imported/{{aid}}/{{sample}}.{{aid}}.p_ctg.gfa",
        fa=f"{RES}/{{sample}}/imported/{{aid}}/{{sample}}.{{aid}}.p_ctg.fa",
        log=f"{RES}/{{sample}}/imported/{{aid}}/{{sample}}.{{aid}}.hifiasm.stderr.log",
    params:
        register=lambda wc: (f"--register {q(IMPORTED[(wc.sample, f'imported-{wc.aid}')]['register'])}"
                             if IMPORTED[(wc.sample, f"imported-{wc.aid}")].get("register") else ""),
    shell:
        "{TOOL} import-assembly --sample {wildcards.sample} --aid {wildcards.aid} "
        "--gfa {input.gfa:q} --fa {input.fa:q} --log {input.log:q} "
        "--link-gfa {output.gfa:q} --link-fa {output.fa:q} --link-log {output.log:q} "
        "{params.register} --out {output.record:q}"
