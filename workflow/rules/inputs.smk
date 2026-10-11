# Input inventory, validation-identity stamp, and full-file checksums.
#
# Checksums and read statistics are computed once and reused: they depend
# only on the (unchanged) input files, never on the validation identity.

localrules: inventory, identity_stamp


rule identity_stamp:
    """Record what the current validation identity consists of (auditable)."""
    output:
        STAMP,
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as fh:
            json.dump({"validation_id": VID, "components": VID_PARTS,
                       "science_version": SCI["science_version"]}, fh, indent=2)
            fh.write("\n")


rule inventory:
    """Record the discovered inputs and selected method (identity of the input set)."""
    output:
        f"{RES}/{{sample}}/inputs/inventory.json",
    run:
        rec = dict(INV[wildcards.sample])
        rec["science_version"] = SCI["science_version"]
        rec["method_file"] = os.path.relpath(method(rec["selected_method"])["_file"], PIPELINE_DIR)
        rec["method_fingerprint"] = METHOD_FP[rec["selected_method"]]
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as fh:
            json.dump(rec, fh, indent=2)
            fh.write("\n")


def checksum_path(s, tech, iid):
    return f"{RES}/{s}/inputs/md5/{tech}/{iid}.md5"


def all_checksums(s):
    return [checksum_path(s, i.tech, i.input_id) for i in used(s, "hifi") + used(s, "ont")]


rule checksum:
    """Full-file MD5 of an original input (read-only); computed once, recorded for provenance/SRA."""
    input:
        reads=lambda wc: get_input(wc.sample, wc.tech, wc.iid).path,
        inventory=ancient(f"{RES}/{{sample}}/inputs/inventory.json"),
    output:
        f"{RES}/{{sample}}/inputs/md5/{{tech}}/{{iid}}.md5",
    log:
        prov=f"{RES}/{{sample}}/provenance/checksum/{{tech}}.{{iid}}.json",
    wildcard_constraints:
        tech="hifi|ont",
    threads: R("checksum")["threads"]
    resources:
        mem_mb=R("checksum")["mem_mb"],
    params:
        cmd=lambda wc, input, output: wrap(
            "checksum", f"{tool('md5sum')} {q(input.reads)} > {q(output[0])}", wc.sample,
            prov_path(wc.sample, "checksum", f"{wc.tech}.{wc.iid}"), tools=["md5sum"],
            inputs=[input.reads], outputs=[output[0]]),
    shell:
        "{params.cmd}"
