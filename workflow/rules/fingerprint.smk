# Optional genotype fingerprinting (k-mer sample-identity QC; docs/fingerprinting.md).
#
#   snakemake --profile profiles/biohpc --cores 32 fingerprint
#
# Independent of assembly: it reads the original input files directly (BAM,
# FASTQ/FASTA, ONT tar), never converts or copies reads, and nothing in the
# assembly DAG depends on it. Every ONT read set is compared with every HiFi
# genotype in the run; HiFi labels define the reference panel, ONT labels
# (sample directories) are the expectations being tested.
#
# Outputs are keyed by a hash of the settings and code that produce them, so a
# change never overwrites an earlier result set:
#   results/fingerprint/sketches-<SID>/<sample>/<tech>/<input_id>.{sig.gz,json}
#   results/fingerprint/report-<RID>/  (assignments.tsv, similarity_matrix.tsv, ...)
import hashlib as _fphash

FPCFG = CFG["fingerprint"]


def _fp_id(obj, code):
    h = _fphash.sha256(json.dumps(obj, sort_keys=True, default=str).encode())
    for f in code:
        with open(os.path.join(PIPELINE_DIR, "scripts", "g2f", f), "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:10]


_FP_SKETCH_KEYS = ("ksizes", "scaled", "hifi_max_bases_per_input", "ont_max_bases_per_input",
                   "read_fraction", "sampling_seed", "min_read_len")
FP_SID = _fp_id({k: FPCFG.get(k) for k in _FP_SKETCH_KEYS}, ["fpsketch.py", "archive.py"])
FP_SAMPLES = SAMPLES
FP_RID = _fp_id({"sketches": FP_SID, "samples": FP_SAMPLES, "analysis": FPCFG.get("analysis"),
                 "primary_k": FPCFG.get("primary_k")}, ["fingerprint.py", "fpanalysis.py"])
FP_DIR = f"{RES}/fingerprint"
FP_REPORT = f"{FP_DIR}/report-{FP_RID}"


def fp_sketch_json(s, tech, iid):
    return f"{FP_DIR}/sketches-{FP_SID}/{s}/{tech}/{iid}.json"


def fp_all_sketches():
    return [fp_sketch_json(s, i.tech, i.input_id) for s in FP_SAMPLES
            for i in used(s, "hifi") + used(s, "ont")]


def fp_sketch_cmd(wc, input, output, threads):
    i = get_input(wc.sample, wc.tech, wc.iid)
    cap = FPCFG.get(f"{wc.tech}_max_bases_per_input", 0)
    cmd = (f"{tool('fingerprint_python')} {q(os.path.join(PIPELINE_DIR, 'scripts', 'g2f_fingerprint.py'))} "
           f"sketch --input {q(i.path)} --format {i.format} --tech {wc.tech} "
           f"--label {q(wc.sample)} --run-id {q(wc.iid)} "
           f"--ksizes {','.join(str(k) for k in FPCFG['ksizes'])} --scaled {FPCFG['scaled']} "
           f"--max-bases {cap} --read-fraction {FPCFG.get('read_fraction', 1.0)} "
           f"--min-read-len {FPCFG.get('min_read_len', 0)} --seed {FPCFG.get('sampling_seed', 1)} "
           f"--samtools {q(tool('samtools'))} --threads {threads} "
           f"--out-sig {q(output.sig)} --out-json {q(output.json)}")
    return released_guard(wc.sample, wc.tech, wc.iid, wrap("fingerprint_sketch", cmd, wc.sample,
                f"{FP_DIR}/sketches-{FP_SID}/{wc.sample}/{wc.tech}/{wc.iid}.prov.json",
                tools=["fingerprint_python"], inputs=[i.path],
                outputs=[output.sig, output.json]))


rule fingerprint_sketch:
    """Streaming FracMinHash sketch of one original read file (bounded prefix)."""
    input:
        reads=lambda wc: raw_input(wc.sample, wc.tech, wc.iid),
    output:
        sig=f"{FP_DIR}/sketches-{FP_SID}/{{sample}}/{{tech}}/{{iid}}.sig.gz",
        json=f"{FP_DIR}/sketches-{FP_SID}/{{sample}}/{{tech}}/{{iid}}.json",
    wildcard_constraints:
        tech="hifi|ont",
    threads: R("fingerprint_sketch")["threads"]
    resources:
        mem_mb=R("fingerprint_sketch")["mem_mb"],
    params:
        cmd=fp_sketch_cmd,
    shell:
        "{params.cmd}"


rule fingerprint_analyze:
    """Compare every ONT read set with every HiFi genotype; calibrate on HiFi controls."""
    input:
        fp_all_sketches(),
    output:
        assignments=f"{FP_REPORT}/assignments.tsv",
        matrix=f"{FP_REPORT}/similarity_matrix.tsv",
        summary=f"{FP_REPORT}/qc_summary.md",
        prov=f"{FP_REPORT}/parameters_and_provenance.json",
    log:
        f"{FP_REPORT}/analysis.log",
    threads: 1
    resources:
        mem_mb=R("fingerprint_analyze")["mem_mb"],
    params:
        ks=",".join(str(k) for k in FPCFG.get("ksizes", [21])),
        primary=FPCFG.get("primary_k", 21),
        analysis=lambda wc: q(json.dumps(FPCFG.get("analysis") or {})),
        prov=lambda wc: q(json.dumps({"report_id": FP_RID, "sketch_id": FP_SID,
                                      "samples": FP_SAMPLES, "settings": FPCFG,
                                      "settings_file": config.get("fingerprint")})),
        outdir=FP_REPORT,
        script=os.path.join(PIPELINE_DIR, "scripts", "g2f_fingerprint.py"),
        listing=f"{FP_REPORT}/sketch_records.txt",
    shell:
        "mkdir -p {params.outdir:q} && printf '%s\\n' {input:q} > {params.listing:q} && "
        "/usr/bin/time -v {PY:q} {params.script:q} analyze "
        "--records-from {params.listing:q} --outdir {params.outdir:q} --ksizes {params.ks} "
        "--primary-k {params.primary} --params {params.analysis} --provenance {params.prov} "
        "> {log:q} 2>&1"
