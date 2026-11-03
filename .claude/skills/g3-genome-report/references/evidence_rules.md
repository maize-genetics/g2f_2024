# Evidence rules

## Status labels

**Evidence status** (each row of the provenance table)
- `VERIFIED`: a machine record or recomputation was read and agrees.
- `RECORDED`: a project record states it, with its code and outputs retained, but it was not recomputed here.
- `RECONSTRUCTED`: described after the fact; no complete machine record.
- `OPEN`: not done, not measured, or awaiting a decision.

**Text status** (each manuscript section)
- `READY`: verified text ready for integration.
- `PROVISIONAL`: drafted, but needs a named check or decision.
- `COLLABORATOR`: content pending from a contributor.
- `MISSING`: an analysis, decision or fact does not exist yet.

**Resource class:** final resource; companion resource; exploratory or rejected; unresolved. Do not move an item between classes in prose; it takes a recorded decision.

## Reporting traps in genome papers

- **File existence is not success.** Cite a validation marker or a recomputation.
- **Reference-guided is not validated.** Order and orientation taken from a reference cannot then be cited as evidence of collinearity with that reference. Report contig contiguity first, joins with independent support separately, and untested joins as untested.
- **Chromosome-scale is not chromosome-level validated.** Without Hi-C, optical maps or genetic maps from the same material, say "chromosome-assigned" or "reference-guided".
- **Scaffold N50 is not contiguity.** Report contig N50, and a supported-join N50 if available.
- **Gaps of fixed size are not distances.**
- **Partial phasing is not haplotype resolution.** State the scale over which phase holds and whether labels connect across contigs.
- **BUSCO duplication is not automatically contamination or retained alleles.** Compare with a reference of the same species and lineage version.
- **Read-derived QV is not independent of the assembly reads.** Name the k-mer source.
- **Identity tests have a scope.** A panel of the study's own lines cannot detect two plants of one accession, or confirm a line against external material.
- **Length thresholds matter.** Report the coverage retained after filtering, not only input yield. Use "ultra-long" only if the data justify it.
- **Tallies from different methods differ.** Say which method a count comes from, and explain disagreements.
- **Superlatives** need a comparison and a source, or they go.
- **A pilot or ongoing analysis** stays out of the released workflow unless a recorded decision adopts it.
- **A QC summary can lag a corrected file.** Check which version of the sequence each metric was measured on.

## Provenance table columns

`item_id`, `manuscript_section`, `claim_or_parameter`, `status`, `scope`, `input_evidence`, `command_or_config`, `software_version`, `date`, `git_commit`, `notes`

One row per reportable method or result. Cite exact paths. Write `unresolved` rather than leaving a cell ambiguous.
