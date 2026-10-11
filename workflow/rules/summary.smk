# Per-sample validated summary and the project-wide table. Only validated
# JSON records feed these tables; the marker file is written last.

PROJECT_SUMMARY = f"{RES}/summary/assembly_summary.tsv"

localrules: sample_summary, project_summary


def validated_marker(s, m):
    return f"{RES}/{s}/summary/{s}.{m}.validated.json"


def reads_qc_targets(s):
    out = [inventory_path(s), hifi_validated(s), kmer_screen_path(s)]
    if uses_ont(METHOD[s]):
        out.append(ont_validated(s, method(METHOD[s])["ont"]["min_read_len"]))
    return out


def summary_inputs(wc):
    s, m = wc.sample, wc.method
    d = {
        "hifi": hifi_validated(s),
        "kmer": kmer_gate_path(s),
        "asm_p": f"{qc_dir(s, m)}/stats/p_ctg.validated.json",
        "asm_a": f"{qc_dir(s, m)}/stats/a_ctg.validated.json",
        "small": f"{qc_dir(s, m)}/stats/p_ctg.small_contigs.tsv",
        "busco": [f"{qc_dir(s, m)}/busco/{l}.validated.json" for l in busco_lineages()],
        "hifiasm_log": hifiasm_log(s, m),
        "gfa": asm_gfa(s, m, "p_ctg"),
    }
    if SCI["qc"]["merqury"].get("enabled", True):
        d["merqury"] = f"{qc_dir(s, m)}/merqury/merqury.validated.json"
    if uses_ont(m):
        d["ont"] = ont_validated(s, method(m)["ont"]["min_read_len"])
    return d


rule sample_summary:
    input:
        unpack(summary_inputs),
    output:
        tsv=f"{RES}/{{sample}}/summary/{{sample}}.{{method}}.summary.tsv",
        marker=f"{RES}/{{sample}}/summary/{{sample}}.{{method}}.validated.json",
    params:
        ont=lambda wc, input: f"--ont {q(input.ont)}" if "ont" in input.keys() else "",
        merqury=lambda wc, input: f"--merqury {q(input.merqury)}" if "merqury" in input.keys() else "",
        busco=lambda wc, input: " ".join(q(b) for b in input.busco),
        kmer_screen=lambda wc: kmer_screen_path(wc.sample),
    shell:
        "{TOOL} summarize --sample {wildcards.sample} --method {wildcards.method} "
        "--science-version {SCI[science_version]} --hifi {input.hifi:q} {params.ont} "
        "--kmer {params.kmer_screen:q} --asm-p {input.asm_p:q} --asm-a {input.asm_a:q} "
        "--busco {params.busco} {params.merqury} --hifiasm-log {input.hifiasm_log:q} "
        "--out-tsv {output.tsv:q} --out-json {output.marker:q}"


rule project_summary:
    input:
        lambda wc: [f"{RES}/{s}/summary/{s}.{METHOD[s]}.summary.tsv" for s in SAMPLES],
    output:
        PROJECT_SUMMARY,
    shell:
        "{TOOL} combine --inputs {input:q} --out {output:q}"
