# Post-assembly QC plan (2026-10-12)

Planned checks after the nine-genotype full run
(`logs/assembly_20261012_023729`) has finished. Nothing here starts before
then. Each item states whether it uses existing workflow outputs (cheap,
read-only) or needs new computation (approval first).

## 1. Input completeness (all nine genotypes)

Read-only, from existing logs and manifests:

* HiFi: hifiasm `# bases` = sum of per-file bases in
  `results/<LINE>/inputs/prepared_inputs.verified.json`. Already confirmed for
  FBLL and GEMS-0227; still to do for the other seven.
* ONT: hifiasm `# passed UL reads` (from
  `results/<LINE>/assembly/*/<LINE>.asm.stderr.log`) = sum of retained ONT
  reads over all prepared ONT files. Flag any discrepancy.

## 2. GEMS-0061 and GEMS-0113 vs the other seven lines

Background: fingerprinting supports all nine genotype assignments, but
GEMS-0061 and GEMS-0113 show genetic discordance between their HiFi and ONT
sources (`docs/fingerprinting_diagnostics_2026-10-12.md`). The project lead
decided on 2026-10-12 to keep the HiFi+ONT assemblies of both lines and to
check for consequences here. The question is whether these two lines are
outliers relative to the other seven; seven lines are a small reference
group, so differences are reported descriptively, not as tests.

| check | source | cost |
|---|---|---|
| Contiguity (contig N50/L50, contig count, total and small-contig length, p_ctg vs a_ctg size) | `results/<LINE>/qc/<akey>/stats/*.validated.json`, `p_ctg.small_contigs.tsv` | read-only |
| BUSCO completeness and duplication (`poales_odb12.2`) | `results/<LINE>/qc/<akey>/busco/*.validated.json` | read-only |
| k-mer completeness and QV (HiFi read k-mers vs assembly) | Merqury outputs, meryl k-mer screen | read-only |
| UL reads integrated by hifiasm (passed / used UL reads as a fraction of filtered ONT reads and bases) | hifiasm stderr logs (§1) | read-only |
| HiFi/ONT alignment consistency (e.g. mapping rate, identity and clipping of ONT vs HiFi reads to p_ctg; regions supported by one technology only) | new alignments (minimap2) | **needs approval**: substantial CPU and disk; script and resource estimate to be prepared first |

Interpretation notes:

* BUSCO duplication alone does not indicate haplotig retention in maize
  (ancient genome duplication); read it together with a_ctg size and k-mer
  multiplicity.
* Merqury completeness is measured against HiFi k-mers, so it cannot by
  itself show ONT-derived discordance; the UL integration fraction and the
  alignment check are the more direct tests.
* An outlier result would justify, in this order, a closer look at the
  affected regions and then reconsidering a HiFi-only comparison assembly
  for the affected line(s). Neither is planned now.

## 3. Not planned unless the QC above justifies it

* HiFi-only assemblies of GEMS-0061 or GEMS-0113.
* Read-level spike-in or other mixing experiments.
* Changes to the fingerprinting module or its thresholds (proposals in §7
  of the diagnostics report, for a later versioned change).
