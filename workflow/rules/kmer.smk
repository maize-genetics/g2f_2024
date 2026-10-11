# HiFi read k-mer database (reused by Merqury), histogram, and a screening
# gate that halts assembly for human review if the spectrum suggests
# residual heterozygosity.


def meryl_db(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.meryl"


def kmer_screen_path(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.screen.json"


def kmer_gate_path(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.gate.json"


localrules: kmer_screen, kmer_gate


rule meryl_count:
    input:
        reads=lambda wc: all_hifi_reads(wc.sample),
        hifi_ok=lambda wc: hifi_validated(wc.sample),
    output:
        directory(f"{RES}/{{sample}}/kmer/hifi.k{K}.meryl"),
    log:
        prov=f"{RES}/{{sample}}/provenance/meryl_count/k{K}.json",
    threads: R("meryl_count")["threads"]
    resources:
        mem_mb=R("meryl_count")["mem_mb"],
    params:
        cmd=lambda wc, input, output, threads, resources: wrap(
            "meryl_count",
            f"{tool('meryl')} count k={K} threads={threads} "
            f"memory={max(1, resources.mem_mb // 1024)} output {q(output[0])} "
            + " ".join(q(r) for r in input.reads),
            wc.sample, prov_path(wc.sample, "meryl_count", f"k{K}"),
            tools=["meryl"], inputs=list(input.reads), outputs=[output[0]]),
    shell:
        "{params.cmd}"


rule meryl_histogram:
    input:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.meryl",
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.hist",
    params:
        meryl=tool("meryl"),
    shell:
        "{params.meryl} histogram {input:q} > {output:q}"


rule kmer_screen:
    input:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.hist",
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.screen.json",
    params:
        lo=SCI["kmer"]["screen"]["half_peak_window"][0],
        hi=SCI["kmer"]["screen"]["half_peak_window"][1],
        ratio=SCI["kmer"]["screen"]["max_half_peak_height_ratio"],
    shell:
        "{TOOL} kmer-screen --hist {input:q} --window-lo {params.lo} --window-hi {params.hi} "
        "--max-ratio {params.ratio} --out {output:q}"


rule kmer_gate:
    """Fails (blocking assembly) if the screen is flagged and no review is recorded.

    samples.yaml is read via params, not declared as an input, so editing it
    never triggers re-assembly of other samples."""
    input:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.screen.json",
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.gate.json",
    params:
        samples=SAMPLES_YAML,
        block="--block" if SCI["kmer"]["screen"].get("block_assembly_on_flag", True) else "",
    shell:
        "{TOOL} kmer-gate --screen {input:q} --sample {wildcards.sample} "
        "--samples {params.samples:q} {params.block} --out {output:q}"
