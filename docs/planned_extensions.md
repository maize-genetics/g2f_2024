# Planned extensions (not implemented)

**None of the functionality on this page exists yet.** The current
deliverable is the standardized assembly and QC workflow. The items below
are planned and will each be added only after they have been validated on
real data (see `AGENTS.md`).

| Extension | Purpose | Status |
|---|---|---|
| EDTA | Repeat and transposable-element annotation | Planned |
| GeneCAD | Gene annotation | Planned |
| PlantCAD / PlantCAD2 | Regulatory (cis-regulatory element) annotation | Planned |
| Comparative and pangenome analyses | Comparing the nine G2F lines with each other and with reference maize genomes | Planned |
| hifiasm with ropebwt3 | Possible alternative assembly approach | Candidate only: to be benchmarked separately first; not added to the production pipeline unless validated |

## Design for annotation modules

EDTA (repeats), GeneCAD (genes), PlantCAD/PlantCAD2 (CREs) and comparative
or pangenome analyses will be added as further rule files
(`workflow/rules/annotation_*.smk`) in this same workflow. They will follow
the conventions already used for assembly.

**Inputs.** A module consumes only *validated* assemblies:

* `asm_fasta(sample, method_id, "p_ctg")`, gated on that assembly's current
  `validated_marker()`;
* or imported assemblies, explicitly and labelled as such.

**Methods.** Each annotation tool gets its own immutable method files, e.g.
`config/methods/edta-<ver>_<params>_v1.yaml`, with an expected version that
is checked when the job runs.

**Outputs.** They go to
`results/<LINE>/annotation/<assembly method_id>/<annotation method_id>/`,
carry provenance records, and are summarized from validated records only.

**Shared metadata.** Sample IDs, aliases, inventories and decisions come
from the same `samples.yaml` and discovery code.

## New assembly methods

A new assembler or assembler mode, such as hifiasm with ropebwt3, would be
added as a new immutable method file in `config/methods/` with a new
`method_id`. Its results would then go to their own
`results/<LINE>/assembly/<method_id>/` directory, alongside the existing
results and never replacing them. Before entering production it would be
benchmarked against the validated `-l0` baselines on PHBW8, and the
decision would be recorded in `docs/scientific_decisions.md`.
