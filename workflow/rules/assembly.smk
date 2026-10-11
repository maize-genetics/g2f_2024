# hifiasm assembly (method selected per sample) and GFA -> FASTA conversion.
#
# Assembly outputs are write-protected after success. Snakemake will refuse
# to overwrite them; a human must move them aside deliberately.


def asm_dir(s, m):
    return f"{RES}/{s}/assembly/{m}"


def asm_gfa(s, m, ctg):
    return f"{asm_dir(s, m)}/{s}.asm.bp.{ctg}.gfa"


def asm_fasta(s, m, ctg):
    return f"{asm_dir(s, m)}/{s}.{m}.{ctg}.fa"


def hifiasm_log(s, m):
    return f"{asm_dir(s, m)}/{s}.asm.stderr.log"


def hifiasm_inputs(wc):
    s, m = wc.sample, wc.method
    d = {
        "hifi": all_hifi_reads(s),
        "hifi_ok": hifi_validated(s),
        "kmer_gate": kmer_gate_path(s),
        "md5": all_checksums(s),
        "inventory": inventory_path(s),
    }
    if uses_ont(m):
        if not used(s, "ont"):
            raise WorkflowError(f"method {m} needs ONT reads, but sample {s} has no usable ONT inputs")
        L = method(m)["ont"]["min_read_len"]
        d["ont"] = ont_merged(s, L)
        d["ont_ok"] = ont_validated(s, L)
    return d


def hifiasm_threads(wc):
    return R("hifiasm_hifi_ont" if uses_ont(wc.method) else "hifiasm_hifi")["threads"]


def hifiasm_mem(wc):
    return R("hifiasm_hifi_ont" if uses_ont(wc.method) else "hifiasm_hifi")["mem_mb"]


def hifiasm_cmd(wc, input, threads):
    s, m = wc.sample, wc.method
    meth = method(m)
    argv = [tool("hifiasm"), "-o", f"{s}.asm", "-t", str(threads)] + [str(a) for a in meth["args"]]
    if uses_ont(m):
        argv += ["--ul", q(input.ont)]
    argv += [q(p) for p in input.hifi]
    run = " ".join(argv)
    d = asm_dir(s, m)
    return wrap(
        "hifiasm",
        f"cd {q(d)} && "
        f"if ls {s}.asm.*.bin >/dev/null 2>&1; then echo 'Refusing to start: intermediate "
        f"hifiasm .bin files from an earlier attempt exist in {d}. Inspect them and move them "
        f"aside before re-running (hifiasm would silently reuse them).' >&2; exit 2; fi && "
        f"{run} > {s}.asm.stdout.log 2> {s}.asm.stderr.log",
        s, prov_path(s, "hifiasm", m), tools=["hifiasm"], method_id=m,
        expect={"hifiasm": meth["expected_version"]},
        inputs=list(input.hifi) + ([input.ont] if uses_ont(m) else []),
        outputs=[asm_gfa(s, m, "p_ctg"), asm_gfa(s, m, "a_ctg")])


rule hifiasm:
    input:
        unpack(hifiasm_inputs),
    output:
        p=protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.p_ctg.gfa"),
        a=protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.a_ctg.gfa"),
    log:
        prov=f"{RES}/{{sample}}/provenance/hifiasm/{{method}}.json",
        stderr=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.stderr.log",
    threads: hifiasm_threads
    resources:
        mem_mb=hifiasm_mem,
    params:
        cmd=hifiasm_cmd,
    shell:
        "mkdir -p $(dirname {output.p:q}) && {params.cmd}"


rule gfa_to_fasta:
    """hifiasm-documented conversion of GFA segment lines to FASTA."""
    input:
        f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.{{ctg}}.gfa",
    output:
        protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.{{ctg}}.fa"),
    log:
        prov=f"{RES}/{{sample}}/provenance/gfa_to_fasta/{{method}}.{{ctg}}.json",
    resources:
        mem_mb=R("gfa_to_fasta")["mem_mb"],
    params:
        cmd=lambda wc, input, output: wrap(
            "gfa_to_fasta",
            "awk '/^S/{print \">\"$2; print $3}' " + f"{q(input[0])} > {q(output[0])}",
            wc.sample, prov_path(wc.sample, "gfa_to_fasta", f"{wc.method}.{wc.ctg}"),
            method_id=wc.method, inputs=[input[0]], outputs=[output[0]]),
    shell:
        "{params.cmd}"


rule check_gfa_fasta:
    """Same names, order, lengths and sequence hashes in GFA and FASTA; ACGTN only."""
    input:
        gfa=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.{{ctg}}.gfa",
        fa=f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.{{ctg}}.fa",
    output:
        f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.{{ctg}}.gfa_check.json",
    resources:
        mem_mb=R("validate_fasta")["mem_mb"],
    shell:
        "{TOOL} check-gfa --gfa {input.gfa:q} --fasta {input.fa:q} --out {output:q}"
