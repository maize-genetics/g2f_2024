# Rerunning individual QC components

Checked against the implementation on 2026-10-11 using dry runs on the real
PHBW8 results, without running any QC. All commands run from `pipeline/`
after `module load python/3.12.7`.

## Summary

| Question | Answer |
|---|---|
| Can BUSCO run with another lineage without redoing read QC or hifiasm? | **Yes**, for any assembly whose FASTA already exists |
| Are the lineage and BUSCO version recorded? | **Yes**. The output path contains the lineage. The validated record holds the lineage, the `dataset.cfg` contents, the pins and the BUSCO version. The summary tables have one row per lineage (`busco_lineage`, `busco_datasets_version`, `busco_dataset`) |
| Do previous results stay available? | **Yes, when a lineage is added.** Replacing a lineage keeps its files on disk, but newer summaries no longer include it. Rerunning the *same* lineage overwrites it (see limitations) |
| Is there a "QC only" target? | **No.** `imported_qc` exists for historical assemblies. For workflow assemblies, request the QC files directly and add the `--allowed-rules` guard shown below |

## Adding a BUSCO lineage

1. **Install the dataset** next to `poales_odb12.2` (needs network
   access). Replace `<lineage>` with the name of a lineage BUSCO provides:

   ```bash
   apptainer exec --bind /home/$USER/busco_downloads /programs/busco-6.1.0/busco.sif \
     busco --download_path /home/$USER/busco_downloads --download <lineage>
   cat /home/$USER/busco_downloads/lineages/<lineage>/dataset.cfg
   ```

2. **Add the lineage under `qc.busco.lineages` in `config/science.yaml`.**
   Keep `poales_odb12.2`. Copy the pins from the `dataset.cfg` you just
   inspected; never type them from memory. The validator checks every pin
   against `dataset.cfg` and the BUSCO result.

   ```yaml
       lineages:
         poales_odb12.2:
           ...                      # unchanged
         <lineage>:
           datasets_version: ...    # from dataset.cfg / the BUSCO download
           creation_date: "..."
           OrthoDB_version: "..."
           dataset_version: "..."
           number_of_BUSCOs: "..."
   ```

   Also record the decision in `docs/scientific_decisions.md`. Editing
   `science.yaml` changes the validation identity, so new summaries get new
   file names and the old ones are kept.

3. **Dry-run and check the job list:**

   ```bash
   snakemake --profile profiles/biohpc --dry-run imported_qc   # historical assemblies
   snakemake --profile profiles/biohpc --dry-run               # workflow assemblies
   ```

   **What re-runs:**
   * one `busco` job per assembly for the new lineage;
   * cheap jobs that depend on the validation identity: `identity_stamp`,
     `check_gfa_fasta` and `small_contigs` (each re-reads the 2.3 GB
     assembly; minutes), the validations and the summaries.

   **What does not re-run:** read QC, meryl, hifiasm, Merqury, and BUSCO for
   lineages that already have results. Existing BUSCO results are only
   re-validated, because the heavy jobs depend on their validation inputs
   through `ancient()`.

4. **Run it** (inside `screen`). Use the same target you dry-ran:
   `imported_qc` for historical assemblies, or no target for workflow
   assemblies. If any workflow assembly is not finished yet, the run with no
   target also builds it; to avoid that, use the targeted form below.

   ```bash
   snakemake --profile profiles/biohpc --cores 32 imported_qc
   ```

## Targeting one QC output, with a guard

Requesting a QC file directly runs only what that file needs. **If the
assembly does not exist yet, "what it needs" includes read QC and hifiasm.**
`--allowed-rules` turns that case into an error rather than an assembly:

```bash
R=$(realpath ../results)     # targets must use the resolved path (/local/workdir/...)
A=hifiasm-0.25.0_hifi-ont-ul50k_l0_v1     # or imported-legacy_hifi_ont_l0
snakemake --profile profiles/biohpc --dry-run \
  --allowed-rules identity_stamp check_gfa_fasta busco validate_busco \
  -- "$R/PHBW8/qc/$A/busco/<lineage>.validated.json"
```

* The guard was verified on the real results. For the historical assembly it
  plans 4 jobs (`identity_stamp`, `check_gfa_fasta`, `busco`,
  `validate_busco`). For the not-yet-built workflow assembly it exits with
  `Missing input files for rule check_gfa_fasta` and nothing runs.
* Drop `--dry-run` (and add `--cores 32`) to run it.
* Targets must use the resolved path: `$PWD/../results/...` from
  `/workdir/...` does not match the workflow's output paths.

Equivalent targets for the other QC components:

| Component | Target under `$R/<LINE>/qc/<method_id or imported-id>/` | Rules to allow |
|---|---|---|
| BUSCO | `busco/<lineage>.validated.json` | `identity_stamp check_gfa_fasta busco validate_busco` |
| Merqury | `merqury/merqury.validated.json` | `identity_stamp check_gfa_fasta merqury validate_merqury` (also needs `kmer/hifi.k21.meryl` from `reads_qc`) |
| Contiguity | `stats/p_ctg.validated.json` | `identity_stamp check_gfa_fasta assembly_stats validate_assembly_stats` |

## Regenerating an existing result

Snakemake does not re-run a job whose outputs exist. To recompute one on
purpose, first **move the old output directory aside** so that it is kept,
then force the rule:

```bash
mv "$R/PHBW8/qc/$A/busco/poales_odb12.2" "$R/PHBW8/qc/$A/busco/poales_odb12.2.prev-$(date +%F)"
snakemake --profile profiles/biohpc --cores 32 \
  --allowed-rules identity_stamp check_gfa_fasta busco validate_busco \
  --forcerun busco -- "$R/PHBW8/qc/$A/busco/poales_odb12.2.validated.json"
```

`--forcerun` on its own writes over the previous result in place.

## Current limitations and recommended minimal changes

None of these are fixed in this release; they are recorded here so that no
result is lost by accident.

1. **The BUSCO version is not part of the output path.** Results go to
   `busco/<lineage>/`.
   * Upgrading BUSCO (`expected_version` in `science.yaml`, plus the
     container in `config/tools/biohpc.yaml`) does not re-run BUSCO.
     Instead, the old result fails validation with
     `BUSCO version … != required …`, and `--forcerun` would overwrite it.
   * *Minimal change:* write to `busco/busco-<version>/<lineage>/`. That is
     a path change for existing results, so it should be added together
     with a one-off migration (moving the current `poales_odb12.2` results
     into the versioned path) rather than silently orphaning them.
2. **No QC-only target for workflow assemblies.** Running with no target
   (`all`) builds any missing assembly first.
   * *Minimal change:* add a `qc_only` target rule that requests QC outputs
     only for assemblies whose validated FASTA already exists. Until then,
     use the `--allowed-rules` guard above.
3. **Lineages are project-wide.** Every assembly is checked against every
   listed lineage. This is intentional (all lines are compared on the same
   datasets). A line-specific lineage would need a per-sample override,
   which does not exist.
4. **Editing `science.yaml` re-runs the cheap QC jobs that read the
   assembly** (GFA/FASTA check, small-contig listing): a few minutes per
   assembly, never hours.
