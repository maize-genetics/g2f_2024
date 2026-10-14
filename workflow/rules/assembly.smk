# hifiasm assembly (method selected per sample) and GFA -> FASTA conversion.
#
# Assembly outputs are write-protected after success; Snakemake refuses to
# overwrite them. A snapshot of the method definition (with its fingerprint)
# is written next to the assembly; the workflow refuses to continue if the
# method file is later edited in place.
#
# Validation gates are declared ancient(): re-validation under new code or
# science settings never re-triggers an expensive assembly.


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
        "hifi_ok": ancient(hifi_validated(s)),
        "prepared_ok": ancient(prepared_verified(s)),
        "md5": [ancient(p) for p in all_checksums(s)],
        "inventory": ancient(inventory_path(s)),
        "method_file": ancient(method(m)["_file"]),
    }
    if KMER_BLOCKING:
        d["kmer_gate"] = ancient(kmer_gate_path(s))
    if uses_ont(m):
        if not used(s, "ont"):
            raise WorkflowError(f"method {m} needs ONT reads, but sample {s} has no usable ONT inputs")
        L = method(m)["ont"]["min_read_len"]
        d["ont"] = all_ont_ul_reads(s, L)
        d["ont_ok"] = ancient(ont_validated(s, L))
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
        argv += ["--ul", q(",".join(input.ont))]
    argv += [q(p) for p in input.hifi]
    run = " ".join(argv)
    d = asm_dir(s, m)
    # Only after hifiasm exited 0 (within the same job, so never while it runs),
    # and only if the primary-contig GFA has sequence lines and the
    # alternate-contig GFA exists: list, then remove hifiasm's reusable
    # intermediate .bin files (ec/ovlp/uidx; tens of GB). Nothing downstream
    # reads them. A failed or interrupted run keeps them for diagnosis (and the
    # next attempt refuses to start until a person moves them aside).
    tidy = ""
    if CLEANUP.get("remove_hifiasm_bins_after_success"):
        rec = f"{s}.asm.removed_bin_files.txt"
        tidy = (f" && if grep -q -m1 '^S' {s}.asm.bp.p_ctg.gfa && test -f {s}.asm.bp.a_ctg.gfa; then "
                f"{{ date '+removed %F %T after successful hifiasm:'; ls -l {s}.asm.*.bin; }} > {rec} "
                f"2>/dev/null || true; rm -f {s}.asm.*.bin; cat {rec} >&2; "
                f"else echo 'hifiasm outputs not verified: .bin files kept' >&2; fi")
    snap = (f"{TOOL} method-snapshot --method-file {q(meth['_file'])} "
            f"--science-version {SCI['science_version']} --out {q(method_record_path(s, m))}")
    # Reusing .bin files from an interrupted attempt (hifiasm's checkpoints) is a
    # human decision: a person writes {s}.asm.reuse_bins.approved (who, why). It is
    # consumed here, so a later attempt needs a fresh approval, and the files
    # reused are recorded in {s}.asm.reused_bin_files.txt.
    ok, used = f"{s}.asm.reuse_bins.approved", f"{s}.asm.reused_bin_files.txt"
    guard = (f"if ls {s}.asm.*.bin >/dev/null 2>&1; then if test -s {ok}; then "
             f"{{ echo 'reusing hifiasm .bin files from an earlier attempt; approval:'; cat {ok}; "
             f"date '+started %F %T with:'; ls -l {s}.asm.*.bin; }} > {used} && rm -f {ok} && "
             f"cat {used} >&2; else echo 'Refusing to start: intermediate "
             f"hifiasm .bin files from an earlier attempt exist in {d}. Inspect them and move them "
             f"aside before re-running (hifiasm would silently reuse them), or approve their reuse "
             f"in {ok}.' >&2; exit 2; fi; fi")
    return wrap(
        "hifiasm",
        f"cd {q(d)} && {guard} && "
        f"{snap} && {run} > {s}.asm.stdout.log 2> {s}.asm.stderr.log{tidy}",
        s, prov_path(s, "hifiasm", m), tools=["hifiasm"], method_id=m,
        expect={"hifiasm": meth["expected_version"]},
        inputs=list(input.hifi) + (list(input.ont) if uses_ont(m) else []),
        outputs=[asm_gfa(s, m, "p_ctg"), asm_gfa(s, m, "a_ctg"), method_record_path(s, m)])


rule hifiasm:
    input:
        unpack(hifiasm_inputs),
    output:
        p=protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.p_ctg.gfa"),
        a=protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.asm.bp.a_ctg.gfa"),
        method=protected(f"{RES}/{{sample}}/assembly/{{method}}/{{sample}}.{{method}}.method.json"),
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
