# Scientific decisions

The scientific configuration lives in `config/science.yaml` and
`config/methods/`. Its current version is `science_version: 2`.

Every change to that configuration must be recorded here and must bump
`science_version`. A method file must not be edited after it has been used
for a production assembly; create a new `method_id` instead.

## D1. Assembler baseline: hifiasm 0.25.0 with `-l0` (2026-10-11)

* **Methods:**
  * `hifiasm-0.25.0_hifi_l0_v1`: HiFi only.
  * `hifiasm-0.25.0_hifi-ont-ul50k_l0_v1`: HiFi plus ONT ≥ 50 kb via `--ul`.
* **Why:** the lines are highly inbred, and `-l0` (no purging) is the
  documented hifiasm setting for homozygous samples. PHBW8 precedent:
  legacy runs `hifi_only/l0` and `ont_hifi`.
* **No manual coverage overrides.** `--hom-cov` is not used. The
  discrepancy in hifiasm's coverage estimates for PHBW8 is unresolved (see
  the audit). A sensitivity test happens only if independent evidence
  justifies it.
* **Runtime check:** each run verifies that the installed hifiasm reports
  version `0.25.0-r726` and refuses to start otherwise.

## D2. Automatic method selection; ONT is always used when present (2026-10-11)

* **Selection rule:**
  * Usable ONT reads present → HiFi + ONT method.
  * No ONT inputs at all → HiFi-only method.
* **No silent fallback.** If ONT files are present but unusable (bad format,
  raw POD5/FAST5, no reads ≥ the cutoff, duplicates), the workflow stops for
  human review. It never falls back to HiFi-only on its own.
* **Human exclusion.** A human may exclude ONT files in `samples.yaml`, with
  a reason. The sample then uses HiFi-only, and the selection reason records
  that exclusion.
* **Per-sample override.** `method_override` in `samples.yaml` requires an
  `override_reason`.

## D3. ONT ultra-long filter: ≥ 50,000 bp (2026-10-11)

* **Starting value:** the PHBW8 precedent. This is not an immutable rule.
* **Implementation:** `seqkit seq -m 50000`, which keeps reads ≥ 50,000 bp.
  The legacy filter may have been "> 50 kb"; the exact historical command is
  unknown.
* **Reporting:** raw and retained yields, read N50 and approximate coverage
  are reported for every input and for each sample.
* **No minimum coverage.** No minimum ONT coverage threshold is imposed;
  adding one needs a documented justification.
* **Standard for future datasets:** ≥ 50,000 bp, defined by the method's
  `ont.min_read_len`.
* **Pre-filtered inputs (added 2026-10-11, science v2).** An acknowledged
  pre-filtered ONT FASTQ is used as supplied. It is not re-filtered or
  copied. The workflow first checks across the whole file that every read
  is ≥ the method threshold, and fails if any read is shorter. Its yield is
  reported as `supplied` (`ont_input_type: prefiltered`); raw yield fields
  stay empty.

## D3a. PHBW8 ONT input: approved pre-filtered dataset (2026-10-11)

* **Decision.** The project lead approved `PHBW8_ONT_gt50kb.fastq.gz` as a
  derived, already length-filtered ONT input, recorded in `samples.yaml`
  under `acknowledge`.
* **Provenance observed** in the first 200 read headers:
  * run `a6c999f5ed670f6a87f642dcc2be2382abd9c375`;
  * flow cell `PAQ81223`;
  * basecall model `dna_r10.4.1_e8.2_400bps_hac@v4.2.0`;
  * first read started 2024-05-20.
* **Unknown:** the raw basecalled files, the unfiltered yield, the filtering
  command, and the length-boundary convention.
* **Boundary convention.** The historical name suggests "> 50 kb", while
  current methods use "≥ 50,000 bp". Any read set that is > 50,000 is also
  ≥ 50,000, so the difference doesn't justify reprocessing or discarding the
  dataset. The distinction is recorded.

## D4. Duplicate and derived inputs (2026-10-11)

* **Two-stage checks.**
  * *At discovery* (cheap): checks file fingerprints, PacBio movies and read
    groups, and ONT read IDs in the first records of each file.
  * *Before assembly* (full files): checks that no PacBio movie appears in
    two HiFi inputs, that no ONT read ID appears in two inputs, that BAM
    record counts equal the converted read counts, and records MD5
    checksums.
* **BAM is the authoritative source.** A FASTQ/FASTA that shares movies with
  a BAM is presumed to be derived from it. The workflow does **not**
  auto-exclude it, though, because a sample of records cannot prove the file
  holds no extra data. It stops and asks for an explicit `exclude_inputs`
  entry.
* **Pre-filtered ONT files.** Files in which every sampled read is ≥ 50 kb
  need an explicit acknowledgement, because their yield cannot be reported
  as raw yield.

## D5. Genome size for coverage reporting: 2.3 Gb (2026-10-11)

* **Use:** only for reporting approximate coverage. It is never passed to
  any tool.
* **Basis:** the PHBW8 assembly sizes (2.27–2.32 Gb).

## D6. Read k-mer screen: diagnostic, not blocking (2026-10-11, science v2)

* **What it does:** meryl k = 21 on the HiFi reads (the same database Merqury
  uses). The spectrum is smoothed with a moving average over about 10% of
  the main-peak coverage. A sample is **flagged** if either:
  * the smoothed height within 0.35–0.65× of the main peak exceeds 25% of
    the main-peak height; or
  * a smoothed local maximum reaching ≥ 5% of the main peak lies in that
    window.
* **Interpretation (idealised).** Assume Poisson coverage, a diploid genome
  and per-base heterozygosity h ≪ 1/k.
  * Heterozygous k-mers number about 2·k·h·G.
  * Their peak sits at half coverage and is √2 narrower.
  * The height ratio is therefore about 2√2·k·h, i.e. about 59·h for k = 21.
  * A ratio of 0.25 corresponds to h ≈ 0.4%. Each screen reports this
    "idealised heterozygosity equivalent".
  * Real HiFi spectra are over-dispersed, and maize repeats and duplications
    add mass at other multiplicities. **The thresholds are not calibrated**,
    and lower residual heterozygosity (< 0.4%) would generally not be
    flagged by the ratio rule.
* **Evidence for treating it as diagnostic only:**
  * The unsmoothed v1 screen flagged a perfectly homozygous synthetic
    dataset because of sampling noise.
  * hifiasm printed `peak_het: 26` for the same data.
* **Policy.**
  * A flag is a recorded warning (`kmer_assessment: flagged_unreviewed`). It
    does not block assembly (`block_assembly_on_flag: false`).
  * An explicit, auditable review can be recorded in `samples.yaml`
    (`kmer_review`: status, reason, reviewer, date); the sample then shows
    `flagged_reviewed`.
  * Blocking can be turned back on in `science.yaml`.
  * An invalid or uncomputable spectrum is always a hard failure.
* **Not used as evidence:** the PHBW8 37× vs 72× discrepancy in hifiasm's
  logs.
* **Calibration plan:** use the PHBW8 spectrum and further lines, plus a
  model-based tool (e.g. GenomeScope 2, not installed) if needed.

## D7. Assembly QC standards (2026-10-11)

* **Contiguity and content:** seqkit size, contig count, N50, L50, GC and gap
  bases for both p_ctg and a_ctg.
* **Primary BUSCO lineage: `poales_odb12.2`** (decided 2026-10-11, science
  v2), run with BUSCO 6.1.0 and miniprot.
  * **Pinned identity:** `datasets_version: odb12.2` (from the BUSCO JSON,
    which BUSCO derives from the `dataset.cfg` name), `dataset.cfg`
    `creation_date: 2026-05-13`, and the BUSCO version 6.1.0.
  * **Enforcement:** every BUSCO validation checks these; a mismatch fails.
  * **Installed and verified (2026-10-11):** OrthoDB 12.2,
    `dataset_version 01`, 6,188 BUSCOs, archive MD5
    `c2cf3a2db7a23f9b8c0d3d9af54c5132` (checked by BUSCO at download). Every
    field is pinned.
  * The separately installed `poales_odb12` (OrthoDB 12.1, created
    2026-05-22) runs as `odb12` and is not this dataset.
  * **Other lineages:** the historical `liliopsida_odb10` results are kept
    and labelled separately and are not re-run. Results from different
    lineages or releases are not interchangeable.
* **Merqury:** QV and k-mer completeness.
* **Small contigs:** contigs ≤ 1 Mb are listed for organelle and
  contamination follow-up. No screening tool is integrated yet; FCS-GX is
  installed at `/programs/fcs-0.5.5`.
* **Selection principle:** no single metric (N50, BUSCO C) decides between
  assemblies.

## D8. Result identity and rerun policy (2026-10-11)

Snakemake runs with `--rerun-triggers mtime`, so editing code or parameters
never silently re-runs expensive jobs. To keep that from leaving outdated
results labelled as valid, the workflow does the following.

* **Methods are immutable.**
  * Each assembly writes `<LINE>.<method_id>.method.json`, which holds the
    method definition and its fingerprint.
  * If a used method file is edited, the workflow refuses to run.
  * Changed parameters need a new `method_id`, which writes a new results
    directory.
  * The required tool version is part of the method and is checked when the
    job runs.
* **Validation identity (`VID`).** This is a hash of `science.yaml`, every
  method definition, and the validation and summary code.
  * Cheap validation jobs depend on a stamp file named with the VID.
  * Markers and summary tables carry the VID in their file names, so a
    change re-runs only cheap jobs and writes summaries under a new name.
  * Old summaries keep their old names and are never presented as current.
  * Heavy jobs depend on validation results only through `ancient()`, so
    re-validation never re-runs assembly, k-mer counting, BUSCO or Merqury.
* **Sample decision identity (`SID`).** This is a hash of a sample's
  `samples.yaml` entry and is part of its summary names. Editing a review or
  a note refreshes that sample's summary only.
* **Read-statistics code identity.** Every statistics file records the
  identity of the code that produced it. Validation refuses statistics made
  by different code instead of silently re-scanning large files.
* **Inputs.** If a sample's inputs or acknowledgements change after results
  exist, the workflow refuses to run.
* **Unaffected.** Resource and tool-path settings are part of none of these
  identities.
* **Failed runs.**
  * Snakemake removes the outputs of a failed job.
  * The provenance record (with `status: failed`) is kept.
  * hifiasm refuses to start while `.bin` files from a failed attempt exist.
* **Reused, not recomputed:** checksums and full-file statistics.

## D9. Historical assemblies as read-only QC benchmarks (2026-10-11)

* **Registration.** The PHBW8 legacy runs are registered in
  `config/external_assemblies.yaml` and are QC'd as
  `method_id: imported-<id>`, with `assembly_source: imported_legacy`.
* **Separation.** Their results go into a separate table and never into the
  pipeline's assembly table.
* **Protection.**
  * Files are linked into `results/`, never written in place.
  * Their size, mtime and fingerprint are recorded at import.
  * If a registered file later changes, the workflow refuses to run.
* **Scope.** Initially only the HiFi + ONT `-l0` run has `stats` and `busco`
  enabled. Merqury needs the full read k-mer database and is enabled
  separately once authorised.

## D10. ONT tar deliveries are read in place, pass reads only (2026-10-12)

* **Context.** The facility delivered every non-PHBW8 ONT run as one
  uncompressed tar of MinKNOW `fastq_pass/*.fastq.gz` batches (18 tars,
  105 GB). Extracting them would duplicate ~105 GB on a nearly full disk.
* **Decision.** `.tar` is an accepted ONT input format.
  * Members are streamed in natural batch order, straight into the existing
    ONT QC and ≥ 50 kb filter; nothing is extracted.
  * `fastq_fail` members are excluded: MinKNOW's own quality filter, the
    usual definition of the usable read set. None of this project's tars
    contain fail reads.
  * Any other member type stops the workflow.
* **Unchanged:** method selection (D2), the ≥ 50 kb filter (D3), duplicate
  detection (read IDs; D4), and every validation gate.
* **Consequence.** The yields reported for these runs are those of the pass
  reads as delivered.

## D11. Optional k-mer genotype fingerprinting (2026-10-12)

* **Purpose.** Sequence-based evidence that each ONT read set belongs to its
  labelled genotype (ONT labels were metadata-only).
* **Method.** FracMinHash (sourmash 4.9.3, k = 21 primary, k = 31 check,
  scaled = 2000). It uses pairwise "distinct allele" tests against every
  HiFi genotype, normalized by recovery of panel-core k-mers, and a
  threshold calibrated on HiFi leave-one-run-out and split-half controls
  (`docs/fingerprinting.md`).
* **Alternatives considered.**
  * Mash screen / bottom-k MinHash: no abundance-filtered set algebra across
    a panel.
  * Exact k-mer counting of every read set with KMC/meryl: tens of GB per
    read set and full passes over all data.
  * Alignment to a reference genome followed by SNP genotype concordance:
    needs a reference and alignments (hours per run), and uses only
    reference-shared loci.
* **Scope.** QC only. It never relabels or moves files and is not an input
  to assembly. Settings are outside `science.yaml` and outside the
  validation identity.

## D12. Retention, cleanup and raw-copy release (approved 2026-10-12)

* **Priorities:** (1) assemblies, QC, logs, provenance; (2) prepared hifiasm
  inputs (converted HiFi FASTQ, filtered ONT FASTQ); (3) regenerable
  intermediates.
* **hifiasm `.bin` files** are removed after hifiasm exits 0 and its GFA
  outputs are verified; the removal is logged. Failed runs keep them.
* **Prepared HiFi FASTQ:** kept, unless the full-scale launcher's disk gate
  projects that keeping it would endanger later genotypes. Then it is
  deleted after use (`retain_prepared_hifi=false`, recorded per genotype).
* **Local raw BAM/tar copies** are never deleted by the pipeline. Eligibility
  for manual deletion requires:
  * verified prepared inputs;
  * a retained, unchanged prepared file;
  * a recorded raw MD5;
  * a confirmed remote copy.

  After deletion, the workflow uses the recorded inventory for those inputs;
  any other missing input still stops it.
* **Why:** ~1.2 TB free against ~1.4 TB without cleanup
  (`docs/full_scale_run_plan.md`). No result changes.
* **Also decided:**
  * HiFi+ONT `-l0` for all nine genotypes;
  * GEMS-0219 `10474838_9` excluded;
  * GEMS-0061 and GEMS-0113 ONT assignments accepted provisionally
    (fingerprinting ambiguous; cause undetermined; QC to be interpreted
    with this in mind).

## Open questions for review

* What does hifiasm's coverage-estimate discrepancy for PHBW8 mean (audit,
  section 2)?
* Why do 127 complete BUSCOs contain internal stop codons?
* How should the k-mer screen thresholds be calibrated once PHBW8 has a real
  histogram?
