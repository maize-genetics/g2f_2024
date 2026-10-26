#!/usr/bin/env python3
"""Pre-handoff check that a FASTA is acceptable to GeneCAD's input preparation (no GPU, no prediction).

Re-implements, from a reading of plantcad/genecad (record the commit you read with --genecad-commit), what
GeneCAD does with FASTA identifiers and bases before any model runs. GeneCAD code itself is NOT executed.
  - predict.sh `list_sequences`: ID = first whitespace token after '>'; a record is used if a non-empty
    sequence line follows;
  - scripts/extract_fasta.py: Bio.SeqIO record.id; the ID becomes a directory name
    (<out>/<ID>/sequences_<ID>.zarr), a Zarr group "<species_id>/<ID>" (xarray.to_zarr, zarr_format=2) and
    the GFF seqid; upper-case A-Z are unmasked.
The Zarr round trip is made with the same libraries on the first --subset-bp bases of every record.

Requires: biopython, numpy, xarray, zarr (GeneCAD's pins: xarray 2025.1.2, zarr 3.0.x, Python 3.12).

  g2f_genecad_input_check.py --fasta GENOME.fa --species-id SAMPLE_A --workdir /path/to/scratch \
      --genecad-commit <sha> --json out.json
Exit status 0 only if every check passes.
"""
import argparse
import json
import os
import re
import sys


def discovery(path):
    """predict.sh list_sequences: [(has_bases, id)] per header."""
    out, cur, seen = [], None, 0
    for ln in open(path):
        if ln.startswith(">"):
            if cur is not None:
                out.append((seen, cur))
            f = ln[1:].split()
            cur, seen = (f[0] if f else ""), 0
        elif not seen and ln.split():
            seen = 1
    if cur is not None:
        out.append((seen, cur))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--species-id", required=True)
    p.add_argument("--workdir", required=True, help="scratch directory for the Zarr round trip")
    p.add_argument("--genecad-commit", default="unknown")
    p.add_argument("--subset-bp", type=int, default=20_000)
    p.add_argument("--json")
    a = p.parse_args()

    import numpy as np
    import xarray as xr
    from Bio import SeqIO

    fai_path = a.fasta + ".fai"
    if not os.path.exists(fai_path):
        sys.exit(f"missing index {fai_path} (samtools faidx)")
    fai = [x.split("\t")[:2] for x in open(fai_path)]
    ids, lens = [x[0] for x in fai], {x[0]: int(x[1]) for x in fai}
    r = {}
    disc = discovery(a.fasta)
    r["predict_sh_discovery_ids_equal_fai"] = [i for _, i in disc] == ids
    r["predict_sh_no_empty_sequences"] = all(s == 1 for s, _ in disc)
    bio_ids, len_ok, alpha_bad, subset = [], True, {}, []
    for rec in SeqIO.parse(a.fasta, "fasta"):
        s = bytes(rec.seq)
        bio_ids.append(rec.id)
        len_ok &= len(s) == lens.get(rec.id, -1)
        extra = set(s) - set(b"ACGT")
        if extra:
            alpha_bad[rec.id] = sorted(chr(c) for c in extra)
        subset.append((rec.id, s[:a.subset_bp]))
    r["biopython_ids_equal_fai"] = bio_ids == ids
    r["biopython_lengths_equal_fai"] = len_ok
    r["alphabet_ACGT_uppercase_only"] = not alpha_bad
    r["ids_unique"] = len(set(ids)) == len(ids)
    r["ids_case_insensitive_unique"] = len({i.lower() for i in ids}) == len(ids)
    r["ids_charset_A-Za-z0-9_.-"] = all(re.fullmatch(r"[A-Za-z0-9_.\-]+", i) for i in ids)
    ok = True
    for cid, s in subset:
        d = os.path.join(a.workdir, a.species_id, cid)
        os.makedirs(d, exist_ok=True)
        arr = np.frombuffer(s, dtype="S1")
        ds = xr.Dataset({"sequence_tokens": (["sequence"], arr)}, coords={"sequence": np.arange(len(arr))},
                        attrs={"species_id": a.species_id, "chromosome_id": cid})
        path = os.path.join(d, f"sequences_{cid}.zarr")
        ds.to_zarr(path, group=f"{a.species_id}/{cid}", zarr_format=2, consolidated=True, mode="w")
        got = xr.open_datatree(path, engine="zarr", consolidated=False)[f"{a.species_id}/{cid}"].to_dataset()
        ok &= (got.attrs["chromosome_id"] == cid and got.attrs["species_id"] == a.species_id
               and b"".join(got["sequence_tokens"].values.tolist()) == s)
    r["zarr_group_and_dir_roundtrip_all_ids"] = bool(ok)
    res = {"fasta": os.path.abspath(a.fasta), "species_id": a.species_id, "records": len(ids),
           "genecad_repo": "https://github.com/plantcad/genecad", "genecad_commit_read": a.genecad_commit,
           "method": "re-implementation of GeneCAD input handling; GeneCAD code not executed",
           "checks": r, "passed": all(r.values())}
    if alpha_bad:
        res["non_ACGT_examples"] = dict(list(alpha_bad.items())[:5])
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=1)
    print(json.dumps(res))
    sys.exit(0 if res["passed"] else 1)


if __name__ == "__main__":
    main()
