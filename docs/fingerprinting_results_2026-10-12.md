# Genotype fingerprinting results (2026-10-12)

Method: `docs/fingerprinting.md`. Full run:
`results/fingerprint/report-742a05a563/` (sketches `sketches-90f6312d72`, settings
`config/fingerprint.yaml`, code commit ff68d36). All 9 HiFi genotypes form
the panel; all 19 ONT read sets (18 tars + the PHBW8 ≥ 50 kb file) are the
queries. No files were moved or relabelled.

## Verdict

| genotype | ONT read sets | result | strength |
|---|---|---|---|
| 3AZA1 | 1 | **supported** | has 0.998, lacks 0.050 |
| FBLL | 2 | **supported** (both runs, all halves) | has 0.96-0.99 |
| GEMS-0219 | 3 | **supported** for `_9b`, `_9c`; `_9` **insufficient data** (52 reads) | has 0.88-0.91 |
| GEMS-0227 | 3 | **supported** (all runs, all halves) | has 0.91-0.97 |
| PHB47 | 3 | **supported** (all runs, all halves) | has 0.92-0.98 |
| PHBW8 | 1 (pre-filtered) | **supported** | has 0.94 |
| PHP02 | 1 | **supported** | has 0.997 |
| GEMS-0061 | 3 | **ambiguous**: best match GEMS-0061 in every run and half, but reads carry other alleles above the control null | has 0.68-0.80, lacks 0.15-0.20 vs *every* other line |
| GEMS-0113 | 2 | **ambiguous**: best match GEMS-0113 in every run and half; elevated FBLL-distinct alleles | has 0.78-0.84; lacks 0.08-0.11, but 0.18 vs FBLL |

**No ONT dataset is potentially mismatched.** Each run's best match is its
labelled genotype, by a wide margin. "has" is the lowest pairwise score of the
expected genotype's distinguishing alleles (~1 = all present). "lacks" is the
highest score for another genotype's distinguishing alleles (~0 = absent).

## Calibration (k = 21; k = 31 gives identical classifications)

* Positive controls: 176 pairwise scores from 13 leave-one-run-out (6 genotypes with
  ≥ 2 HiFi runs: independent SMRT cells) and 9 split-half controls (technical).
  Min **0.923**, median 0.997. All 22 controls are classified `supported`.
* Negative controls: max **0.114**, median 0.019.
* Threshold **0.518**. The positive and negative distributions are separated
  by 0.81: the method discriminates all nine genotypes, including the most
  closely related (GEMS-0113 shares 81 % of its solid k-mers with FBLL).
* ONT null (label-free): median 0.042, 95th percentile 0.150, max 0.171. The
  upper tail comes entirely from GEMS-0061 and GEMS-0113 (below).

## GEMS-0061 and GEMS-0113: what the evidence says

1. **Not a mislabel and not contamination by a panel line.** For GEMS-0061, the
   foreign-allele signal is uniform across all eight other lines (0.15-0.20).
   Contamination by a specific line would give a high score for that line
   only. For GEMS-0113 it is mostly FBLL (0.18 vs ≤ 0.11), its closest
   relative in the panel. All runs (independent flow cells) and both
   replicate halves of every run agree, so this is a property of the ONT DNA
   source, not of one run.
2. **The HiFi samples themselves are not homozygous.** In the HiFi k-mer
   spectra, the mass at half the single-copy depth relative to the main peak
   (abundance 4-6 vs 8-13, all at 25 Gb sampled) is:

   | line | half / full |
   |---|---|
   | GEMS-0061 | **0.444** |
   | GEMS-0113 | **0.341** |
   | PHP02 | 0.207 |
   | PHBW8 | 0.203 |

   Both GEMS lines also have the largest numbers of distinct sampled k-mers
   (present 488 k and 469 k vs 448-450 k), and the largest split-half
   negative-control scores (0.10, 0.08). This indicates residual
   heterozygosity or segregation (or pooled, non-identical plants) in the
   HiFi DNA.
3. **Most likely explanation:** the ONT and HiFi libraries were made from
   genetically non-identical individuals of a line that is not fully inbred.
   At loci where the line segregates, the ONT plant carries alleles absent
   from the HiFi sample, and vice versa. This explains the uniform has ≈ 0.8
   / lacks ≈ 0.17. These numbers are not uniquely diagnostic: an ONT sample
   from a related but different seed source would look similar. The data
   cannot distinguish the two without pedigree/seed-lot information.

## Decision (project lead, 2026-10-12)

The GEMS-0061 and GEMS-0113 ONT assignments are **accepted provisionally**,
and both lines are assembled HiFi+ONT `-l0` like the others. The cause of the
alternative alleles is **not determined**: residual heterozygosity,
biological differences between the sampled plants, admixture, or other
causes remain possible. It will be discussed with collaborators. No
definitive explanation should be assigned without further evidence. The
options listed below were considered; option 1's HiFi-only comparison was
**not** adopted.

## Implications for assembly (as reviewed before the decision)

* **Label assignments:** all supported or best-matching the label. Nothing
  suggests swapping files.
* **GEMS-0061 and GEMS-0113 HiFi+ONT assembly:** `--ul` reads from a
  non-identical plant can join haplotypes that differ from the HiFi
  sample's. Options:
  1. assemble as planned (HiFi+ONT) and compare with a HiFi-only assembly of
     the same line (`method_override: hifiasm-0.25.0_hifi_l0_v1` in
     `samples.yaml` creates a separate result);
  2. assemble HiFi-only for these two lines;
  3. defer them until seed-lot/plant provenance of the ONT and HiFi DNA is
     checked.
* **Heterozygosity and `-l0`:** D1 chose `-l0` (no purging) for inbred maize.
  These two HiFi samples are measurably heterozygous, so `-l0` primary
  assemblies may retain duplicated haplotigs. The pipeline's k-mer screen
  will report this (diagnostic, D6). Interpret BUSCO duplication for these
  two lines accordingly. Do not change the method silently; record a
  decision.
* **GEMS-0219 `_9`:** insufficient data, as expected (52 reads). Exclude it
  (`docs/full_scale_run_plan.md`).

## Pilot

Bounded pilot (4 genotypes, 8 Gb HiFi / 4 Gb ONT per file;
`config/fingerprint_pilot.yaml`). The original report
(`results/fingerprint/report-91311470e3/`) is **invalid** and kept only as a
record. At ~3.5x, PHBW8's single HiFi file had no resolvable single-copy
peak, so the solid band locked onto a repeat peak, which corrupted the panel
core. The fix anchors the peak to the expected depth, raises the HiFi cap to
25 Gb per file, and refuses scores when the core is too small. Re-analysis
of the same sketches (`results/fingerprint/pilot-reanalysis-ff68d36/`)
classified every ONT run with data as `supported` and the 52-read run as
insufficient. That justified the full run.

## Resources (measured)

| step | value |
|---|---|
| sketch jobs | 36 (17 HiFi, 19 ONT), 3 threads each, ≤ 48 cores |
| bases sketched | 431 Gb (HiFi ≤ 25 Gb per file, ONT ≤ 10 Gb per run) |
| sketch CPU time | 30,735 s (8.5 h); 15.5 Mbp/s per job (HiFi BAM), ~12 Mbp/s (ONT tar) |
| wall time | 54 min (whole Snakemake run) |
| analysis | 6 min 43 s, peak RSS 2.2 GB, single thread |
| disk | all fingerprint sketches and reports (pilot + full) 616 MB; no read copies |

## Addendum: supplementary diagnostics (2026-10-12)

The text above is preserved as written. `docs/fingerprinting_diagnostics_2026-10-12.md`
re-examines GEMS-0061, GEMS-0113 and GEMS-0219 run 9 with read-only diagnostics
(ONT-matched null, k = 21/31 values, cross-run replication, per-read marker
co-occurrence, allele dosage, in-silico admixture). In brief: all nine
genotype assignments are supported. No evidence of contamination by another
panel genotype was found in any ONT set (indicative detection limit about
1–3 %; lower-level, non-panel or non-maize contamination is not excluded).
For GEMS-0061 and GEMS-0113, the foreign alleles are real sequence at
near-genomic dosage, which indicates genetic discordance between the HiFi
and ONT sources and is not explained by low-level admixture. GEMS-0113's
FBLL signal is part panel composition and part FBLL-like haplotype segments;
its dosage argues against FBLL admixture. The biological provenance is not
established. Classifications and thresholds are unchanged.
