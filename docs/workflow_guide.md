# Workflow guide: from reads to an annotation-ready assembly release

This guide describes the whole assembly process the pipeline supports, end to end. It covers what each stage is for,
which entry point runs it, what it must produce before the next stage may start, and which decisions belong to a human.

- [README.md](../README.md) covers installation and the Snakemake part (stages 1–4).
- The machine-readable index of stages is [`config/workflow_stages.yaml`](../config/workflow_stages.yaml).
- The record of a particular project is kept outside the repository, under `../project_records/`. That record holds
  its commands, results and decisions (see [AGENTS.md](../AGENTS.md)).

**Status words used below**

| status | meaning |
|---|---|
| **workflow** | a Snakemake rule; resumable; validated outputs; provenance JSON per step |
| **tool** | a tracked CLI in `scripts/` with tests, run by hand or by a launcher |
| **project script** | logic that so far exists only in a project's records. It was executed, but it is not packaged here. Promote it before relying on it for a new project. |
| **manual** | a human decision or review; no code decides it |
| **planned** | not implemented |

## Species profile: what is not generic

The code is species-agnostic. These values are species- or project-specific and live in configuration, not in code.

| parameter | where it is set | maize default (example) |
|---|---|---|
| Genome size for coverage | `config/science.yaml` `genome_size_bp` | 2.3 Gb |
| BUSCO lineage and pinned dataset | `config/science.yaml` `qc.busco.lineages` | `poales_odb12.2` |
| Expected BUSCO duplication of a correct haploid assembly | human interpretation (see stage 5) | ≈ 25 % (odb12.2): maize is an ancient polyploid |
| ONT ultra-long length cut-off | `config/methods/*.yaml` `ont.min_read_len` | 50,000 bp |
| k-mer size | `config/science.yaml` `kmer.k` | 21 |
| Reference genome and annotation for synteny | `config/local/anchorwave.yaml` (`config/anchorwave.example.yaml`) | B73 RefGen v5 + its GFF3 |
| Reference chromosome names | assignment tables / release naming | `chr1` … `chr10` |
| Known reference artefacts | project records | e.g. reference N-gaps and long tandem repeat arrays at which contigs end |

Add another species by changing these values, not code.

## Stage 0 — Input registration and identity

| | |
|---|---|
| entry | `python3 scripts/g2f_tool.py discover`; Snakemake `inventory`, `checksum`, `identity_stamp` (**workflow**) |
| produces | `results/<S>/inputs/inventory.json`, MD5s, identity stamp; the HiFi-only or HiFi+ONT method choice |
| gate | no duplicate or derived inputs; inputs unchanged since the results were built |

**Sample identity is a scientific question, not a file check.**
- k-mer genotype fingerprinting (`snakemake fingerprint`, [fingerprinting.md](fingerprinting.md); **workflow**)
  tests whether each read set is the genotype it is labelled as.
- Fingerprinting cannot detect two plants of the same accession pooled across sequencing runs.
- For a heterozygous or segregating accession, check per-run (per-SMRT-cell or per-flow-cell) haplotype structure
  before assembly. For read-level genotype forensics, use `scripts/g2f_read_forensics.py` (**tool**).
- HiFi and ONT libraries from different plants of one accession are possible. In that case, record which reads went
  into which assembly.

## Stage 1 — Read preparation and QC (**workflow**)

| | |
|---|---|
| HiFi | `hifi_bam_to_fastq` (pbtk), `hifi_read_stats`, `validate_hifi_reads` |
| ONT | `ont_input_stats`, `ont_filter` (≥ `min_read_len`), `validate_ont_reads`. Tar deliveries are read in place, `fastq_pass` only. |
| gate | `verify_prepared_inputs` must pass before hifiasm may start |

## Stage 2 — k-mer analysis (**workflow**)

| | |
|---|---|
| entry | `meryl_count`, `meryl_histogram`, `kmer_screen`, `kmer_assess` |
| produces | read k-mer database (reused by Merqury), histogram, screen JSON |

- The screen is **diagnostic**: it does not block assembly. A heterozygous-peak flag means the sample is not behaving
  as an inbred, and the stage 3 strategy must be chosen by a human.
- Use the homozygous (main) peak as the single-copy depth `c` in later read-depth analyses.

## Stage 3 — Assembly

| | |
|---|---|
| entry | Snakemake `hifiasm` (**workflow**), one immutable method file per parameter set in `config/methods/` |
| default | hifiasm 0.25.0 `-l0` (no purging), plus `--ul` with ONT ≥ 50 kb when present |
| produces | `p_ctg` / `a_ctg` FASTA (checked against the GFA), the hifiasm logs, and a provenance record (command, version, inputs, Git commit) |

**Choosing a strategy by heterozygosity.** This is decided by a human and recorded in the project records.

1. **Inbred or near-inbred** (no heterozygous k-mer peak; BUSCO D at the species baseline): use `-l0`.
   - Purging can remove real sequence and does not improve base accuracy.
   - Keep `-l0` unless redundancy is shown to be large *and* purging is shown not to lose genes.
2. **Heterozygous** (clear heterozygous peak; inflated size; duplicated genes far above baseline):
   - An `-l0` primary carries both alleles and redundant copies. It is neither haploid nor phased.
   - Evaluate hifiasm's dual output (hap1 + hap2) at a fixed purge level, and compare it with the `-l0` primary on:
     - pair k-mer completeness;
     - homozygous k-mers present twice or absent;
     - BUSCO D per haplotype;
     - where the `-l0` sequence ends up (sequence fate).
   - A purged primary (`-l3`) is a haplotype mosaic. Use it only as a labelled companion, not as a haplotype.
3. **hifiasm coverage-model failures.** If hifiasm misestimates the homozygous coverage, one haplotype collapses.
   - Override with `--hom-cov`, and `--pri-range` where needed.
   - Rerun from the preserved checkpoints. These runs are graph-only and fast once the `.bin` files exist.
4. **Mixed individuals** (stage 0): assemble per individual from the reads that belong to it. Do not use the pooled
   reads.

New parameter sets go in **new** method files with a new `method_id`. Each `method_id` gets its own results directory,
so earlier results are never overwritten. Methods not run through the Snakemake rule are **project script** runs. One
example is a dual-haplotype run with overrides. Record their exact command and checkpoint reuse alongside the outputs.

## Stage 4 — Standard assembly QC (**workflow**)

| check | rule | read it as |
|---|---|---|
| contiguity | `assembly_stats`, `small_contigs` | count, N50, gaps; small contigs listed for follow-up |
| gene space | `busco` (pinned dataset) | completeness. Duplication is relative to the species baseline and is never on its own evidence of haplotigs. |
| base accuracy / completeness | `merqury` | QV and k-mer completeness. For a haplotype pair, report each haplotype and the pair. |

To rerun or extend QC without reassembly, see [rerunning_qc.md](rerunning_qc.md).

## Stage 5 — Structural validation

| | entry | status |
|---|---|---|
| read support, depth classes, graph paths, pairwise sequence fate | `scripts/g2f_asmqc.py` ([assembly_validation.md](assembly_validation.md)) | **tool** |
| gene-anchored synteny to a reference | `scripts/g2f_anchorwave.py` + `g2f_anchorwave_synteny.py` ([anchorwave.md](anchorwave.md)) | **tool** (container) |
| reads across a specific breakpoint | `scripts/g2f_breakpoint_reads.py` | **tool** |
| per-junction call: SUPPORTED / CONTRADICTED / UNRESOLVED | junction evaluation over AnchorWave junctions plus HiFi and ONT spanning reads | **project script** |
| chromosome assignment, telomere scan | assignment from AnchorWave blocks; telomere-array scan | **project script** |

**Rules of evidence**
- A difference from the reference is an *observation*, not an error. Agreement with the reference is not evidence that
  a join is correct.
- Judge a junction by independent long reads that span it uniquely, by its local depth, and by the assembly graph:
  - reads contradict the join → **CONTRADICTED**;
  - unique reads cross it → **SUPPORTED**;
  - no read can cross it (e.g. a long tandem array) → **UNRESOLVED**.
- Recurrent junctions are worth checking against other lines. A junction found at the same breakpoints in several
  independent assemblies, and spanned by reads, is likely a real structural variant (or a reference artefact).
- Run sensitivity variants (AnchorWave parameters, a lower minimum number of anchors) and report calls that depend on
  them.

## Stage 6 — Structural correction (**project script**, sequence-preserving)

Only two operations are allowed. Both are implemented as a split with a verification record.
- **Break** a join that reads CONTRADICT.
- **Split** a join that no read can test, when it sits at a known artefact class (e.g. a long tandem array at a
  reference gap).

Requirements:
- No base is changed, removed or added.
- The pieces must concatenate back to the original.
- The verification record must show identical length and base composition.
- A coordinate map and AGP must relate each piece to the original contig.
- Every split or break is listed with its evidence and its rule.

Joins left UNRESOLVED are documented, not silently kept or removed. Whether to split is a **manual** decision.

## Stage 7 — Selection and human review (**manual**)

Agents prepare evidence tables. Humans decide.
- Never auto-select an assembly when QC evidence conflicts. Examples: better contiguity but more misjoin sites, or
  lower duplication but lost genes.
- Present every alternative with the same metrics, and say what each alternative omits.
- An unexplained check failure is reported as unexplained. It is never relabelled as a pass.
- Distinguish three things:
  - *technical release validation*: checksums, structure, file contracts;
  - *scientific QC*: whether the assembly is good;
  - *collaborator approval*.

## Stage 8 — Release management

| step | entry | status |
|---|---|---|
| package each selected FASTA with its QC, provenance and caveats; checksum it; validate it | the project's release builder | **project script** |
| identifier scheme `<LINE>_chrNN_ctgNNN` / `<LINE>_chrUn_ctgNNN`: chromosome only for HIGH-confidence assignments; ctg = rank by reference start (chromosome tier) or by length (chrUn) | the project's naming script | **project script** |
| reorder records (chr1…chrN numerically, then chrUn), byte-identical records | `scripts/g2f_release_fasta.py sort` | **tool** |
| independent verification of a reordered or renamed FASTA | `scripts/g2f_release_fasta.py verify` | **tool** |
| annotation-tool input gate (GeneCAD) | `scripts/g2f_genecad_input_check.py` | **tool** (re-implementation; also run the real tool's input step when permitted) |
| provenance table from the step records | `scripts/g2f_provenance_report.py` | **tool** |

**Release rules**
1. **Releases are immutable.**
   - A change, even one only to names or record order, creates a new versioned directory.
   - The previous release is re-verified (checksums and a tree snapshot) and left untouched.
   - Make the files read-only after validation.
2. **A derived release says exactly what changed.**
   - It verifies everything else independently:
     - per-record sequence SHA-256;
     - each ID present exactly once;
     - lengths, total bp, line layout;
     - byte-identical record blocks where only the order changed.
   - It carries the mapping tables, keyed by identifier.
3. **Names do not imply scaffolding.**
   - Reference-ordered numbering is an order, not an adjacency.
   - Overlapping reference spans are flagged as ambiguous.
   - Unplaced and provisional contigs carry no chromosome in their names.
4. **Potential alternative-allele contigs are flagged, not removed.** An example is a "two-copy candidate": most of its
   reference span is also covered by a stronger contig.
5. **Every release has an authoritative manifest.** It lists one FASTA per unit, with its path, sha256, record count
   and size, and the source release it was derived from.

## Extending the workflow

- **A new QC method.** Add it as a tool or rule that writes a new, versioned output next to the existing ones. Never
  replace a validated result in place.
  - Record its version and parameters in a provenance JSON.
  - Add a synthetic test.
  - Register it in `config/workflow_stages.yaml`.
- **A new assembly parameter set.** Add a new method file. Do not edit a used one; the pipeline refuses that.
- **Promote a project script.** Move it into `scripts/` with:
  - a CLI;
  - species values read from configuration;
  - a synthetic test;
  - an entry in `config/workflow_stages.yaml`.

  Keep the project's original copy and its checksum in the project records as provenance of the historical run.
- **Run nothing expensive automatically.** Assembly, read mapping, AnchorWave alignment and annotation are started by a
  human, with a dry run and a resource estimate first.
