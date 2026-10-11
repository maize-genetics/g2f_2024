# Scientific decisions

The scientific configuration lives in `config/science.yaml` and
`config/methods/`. Its current version is `science_version: 1`.

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

## D6. Read k-mer screen (2026-10-11)

* **What:** meryl k = 21 on the HiFi reads. The same database is reused by
  Merqury.
* **Heuristic:** a sample is flagged if the histogram has a local maximum,
  or a height above 25% of the main peak, within 0.35–0.65× of the main peak
  coverage.
* **Effect:** a flagged sample does not assemble until `samples.yaml`
  records a `kmer_review`.
* **Status:** a screening heuristic for routing samples to review, not an
  estimate of heterozygosity. It has not been calibrated on real data yet;
  calibrate it on PHBW8.

## D7. Assembly QC standards (2026-10-11)

* **Contiguity and content:** seqkit size, contig count, N50, L50, GC and gap
  bases for both p_ctg and a_ctg.
* **Primary BUSCO lineage:** `poales_odb12` (OrthoDB 12.1, created
  2026-05-22, 6,282 markers), run with miniprot.
  * Available offline in `~/busco_downloads`.
  * BUSCO's lineage list also has `poales_odb12.2` (dated 2026-05-13). Which
    of the two is current needs checking before production.
  * The historical `liliopsida_odb10` results are kept for comparison only.
    Results from different lineages or OrthoDB versions are not
    interchangeable.
* **Merqury:** QV and k-mer completeness.
* **Small contigs:** contigs ≤ 1 Mb are listed for organelle and
  contamination follow-up. No screening tool is integrated yet; FCS-GX is
  installed at `/programs/fcs-0.5.5`.
* **Selection principle:** no single metric (N50, BUSCO C) decides between
  assemblies.

## Open questions for review

* What does hifiasm's coverage-estimate discrepancy for PHBW8 mean (audit,
  section 2)?
* Why do 127 complete BUSCOs contain internal stop codons?
* Use `poales_odb12` or `poales_odb12.2`?
* How should the k-mer screen thresholds be calibrated once PHBW8 has a real
  histogram?
