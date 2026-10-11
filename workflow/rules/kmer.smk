# HiFi read k-mer database (reused by Merqury), histogram, and a diagnostic
# screen for possible residual heterozygosity.
#
# The screen is an uncalibrated heuristic. By default it only records a
# warning (kmer_screen.block_assembly_on_flag: false). A human review in
# samples.yaml (kmer_review) is recorded in the assessment, whose file name
# carries a hash of the review so that editing the review regenerates the
# assessment and summaries, but never the assembly.
import hashlib as _hashlib


def meryl_db(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.meryl"


def kmer_screen_path(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.screen.json"


def review_id(s):
    rev = (CFG["samples"].get(s) or {}).get("kmer_review")
    if not rev:
        return "none"
    return _hashlib.sha256(json.dumps(rev, sort_keys=True, default=str).encode()).hexdigest()[:8]


def kmer_assessment_path(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.assessment.{review_id(s)}.json"


def kmer_gate_path(s):
    return f"{RES}/{s}/kmer/hifi.k{K}.gate.json"


KMER_BLOCKING = bool(SCI["kmer"]["screen"].get("block_assembly_on_flag", False))

localrules: kmer_screen, kmer_assess, kmer_gate


rule meryl_count:
    input:
        reads=lambda wc: all_hifi_reads(wc.sample),
        hifi_ok=lambda wc: ancient(hifi_validated(wc.sample)),
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
            f"OMP_NUM_THREADS={threads} {tool('meryl')} count k={K} threads={threads} "
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
        "OMP_NUM_THREADS=1 {params.meryl} histogram {input:q} > {output:q}"


rule kmer_screen:
    input:
        hist=f"{RES}/{{sample}}/kmer/hifi.k{K}.hist",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.screen.json",
    params:
        lo=SCI["kmer"]["screen"]["half_peak_window"][0],
        hi=SCI["kmer"]["screen"]["half_peak_window"][1],
        ratio=SCI["kmer"]["screen"]["max_half_peak_height_ratio"],
        lmr=SCI["kmer"]["screen"].get("min_local_max_ratio", 0.05),
        smooth=SCI["kmer"]["screen"].get("smooth_frac", 0.1),
    shell:
        "{TOOL} kmer-screen --hist {input.hist:q} --window-lo {params.lo} --window-hi {params.hi} "
        "--max-ratio {params.ratio} --min-local-max-ratio {params.lmr} "
        "--smooth-frac {params.smooth} --k {K} --out {output:q}"


rule kmer_assess:
    """Screen + recorded human review. Fails only for an invalid screen, or (if
    blocking is configured) for a flagged sample without an accepted review."""
    input:
        screen=f"{RES}/{{sample}}/kmer/hifi.k{K}.screen.json",
        stamp=STAMP,
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.assessment.{{reviewid}}.json",
    params:
        samples=SAMPLES_YAML,
        block="--block" if KMER_BLOCKING else "",
    shell:
        "{TOOL} kmer-assess --screen {input.screen:q} --sample {wildcards.sample} "
        "--samples {params.samples:q} --review-id {wildcards.reviewid} {params.block} "
        "--out {output:q}"


rule kmer_gate:
    """Only used when blocking is configured: fixed-name copy of a passing
    assessment, so later review edits never re-trigger assembly."""
    input:
        lambda wc: ancient(kmer_assessment_path(wc.sample)),
    output:
        f"{RES}/{{sample}}/kmer/hifi.k{K}.gate.json",
    shell:
        "cp {input:q} {output:q}"
