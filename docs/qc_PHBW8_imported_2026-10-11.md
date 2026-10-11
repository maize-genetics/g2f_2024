# PHBW8 historical HiFi+ONT `-l0` assembly: first real QC run (2026-10-11)

## What was run

* **Command:** `snakemake --profile profiles/biohpc --cores 32 imported_qc`,
  run on cbsuxm01 at 01:23–01:35 (12.5 min wall time).
* **Assembly:** `imported-legacy_hifi_ont_l0`, i.e. `hifiasm_v0.25.0/ont_hifi/`
  registered read-only.
  * It was **not produced by this workflow**; see
    `docs/legacy/PHBW8_legacy_runs.yaml` for its provenance limits.
  * Sizes and mtimes of the three registered files are identical before and
    after the run, and no file in `hifiasm_v0.25.0/` is newer than the
    project start.
* **Validation identity:** `efa133413e` (science v2). The pipeline commit was
  `6043a0b`, plus the uncommitted working-tree changes that were committed
  immediately afterwards (provenance shows `dirty: true`; see the note at the
  end).
* **Outputs:** `results/PHBW8/qc/imported-legacy_hifi_ont_l0/` and
  `results/summary/imported_assemblies_summary.efa133413e.484c84.tsv`.

## Results (all validation checks passed)

### GFA vs FASTA consistency

* 490 GFA segments and 490 FASTA records, both totalling 2,266,382,187 bp.
* Names, order, lengths and sequence hashes are identical, and the FASTA
  contains only ACGTN.
* This verifies the output of the unrecorded historical GFA-to-FASTA
  conversion.

### Contiguity (seqkit 2.13.0)

| Metric | Value |
|---|---|
| Contigs | 490 |
| Total | 2,266,382,187 bp |
| N50 / L50 | 121,136,088 bp / 6 |
| Max / min | 345,841,590 / 18,047 bp |
| GC | 46.71% |
| Gaps (N) | 0 bp |

These reproduce the audit values exactly.

### BUSCO 6.1.0, `poales_odb12.2`

* **Dataset:** OrthoDB 12.2, created 2026-05-13, n = 6,188. Run offline with
  miniprot 0.18-r281; 728 s wall time at 32 threads.

| Category | Count | % |
|---|---:|---:|
| Complete | 6,123 | 98.9 |
| Single-copy | 4,575 | 73.9 |
| Duplicated | 1,548 | 25.0 |
| Fragmented | 3 | 0.0 (0.05) |
| Missing | 62 | 1.0 |

* **Internal stop codons:** 143 of 6,123 complete matches (2.3%) contain
  internal stop codons in the miniprot predictions.

**This is a different dataset from the historical results.** The historical
`liliopsida_odb10` result for the same assembly was C 99.1 [S 82.2, D 17.0],
F 0.6, M 0.2, n = 3,236, with 127 internal stops (4.0%). The two differ in
marker set, number of markers, OrthoDB release and completeness criteria, so
**the differences, including the higher duplicated fraction, are not changes
in assembly quality** and must not be read as such. As before, duplication
is not interpreted as haplotig retention: maize has an ancient genome
duplication, and a Poales-level marker set is expected to capture more
lineage-specific duplicates.

### hifiasm log (historical run)

* Version 0.25.0-r726.
* Peak RSS 133.7 GB.
* 12,416 UL reads passed.

## Observations needing follow-up (not conclusions)

* **Small contigs.** 445 of the 490 contigs are ≤ 1 Mb and together hold
  44.0 Mb:

  | Size class | Contigs |
  |---|---:|
  | < 100 kb | 314 |
  | 100–500 kb | 129 |
  | 0.5–1 Mb | 2 |

  The remaining 45 contigs hold about 2.22 Gb.
* **Small-contig GC is bimodal.** 174 contigs fall at 36–40% GC and 114 at
  58–60% GC (listed in `stats/p_ctg.small_contigs.tsv`).
  * *Hypotheses only:* plastid-derived sequence (maize chloroplast GC is
    ~38%) and high-GC tandem repeats or rDNA.
  * *Proposed check:* align the small contigs to the maize chloroplast and
    mitochondrial genomes and to rDNA / known repeat consensus sequences
    (light compute). Organelle and contamination screening (e.g. FCS-GX) is
    still planned.
* **Internal stop codons.** These remain unresolved (2.3% with odb12.2, 4.0%
  with odb10). Merqury QV, from the read-based QC stage, will show whether
  base errors contribute.
* **Not recorded for this run:** BUSCO peak memory. The provenance recorder
  now records the peak RSS of the largest child process and CPU times for
  future jobs.

## Note on provenance

The run's provenance records `pipeline_git.dirty: true`. The working tree
held the `poales_odb12.2` pin and the hifiasm UL-count parser, both committed
right after the run. The validation identity `efa133413e`, which covers
`science.yaml`, all method definitions and the validation code, was
recomputed from the committed tree and matches, so the committed code is the
code that produced these results.
