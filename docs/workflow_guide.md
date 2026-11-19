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
| Sequence-anchor markers for reference-guided pseudomolecules | `config/local/allmaps.yaml` (`config/allmaps.example.yaml`); LG labels in `config/allmaps_science.yaml` | pan-genome anchor markers (100 bp, LG 1–10) |
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

**Foreign contamination and organellar contigs** (`scripts/g2f_decontam.py`, **tool**). Screen every assembly that
will be released: FCS-adaptor and FCS-GX for foreign sequence, and a window-based organelle screen for the species' own
plastid and mitochondrial contigs. Removing what they find is a release step; see
[Stage 8d](#stage-8d--removing-contaminant-and-organellar-records-tool).

## Stage 5 — Structural validation

| | entry | status |
|---|---|---|
| read support, depth classes, graph paths, pairwise sequence fate | `scripts/g2f_asmqc.py` ([assembly_validation.md](assembly_validation.md)) | **tool** |
| gene-anchored synteny to a reference | `scripts/g2f_anchorwave.py` + `g2f_anchorwave_synteny.py` ([anchorwave.md](anchorwave.md)) | **tool** (container) |
| whole-genome riparian overview (visual QC) | `scripts/g2f_riparian.py` from existing AnchorWave anchors; lifts anchors through split coordinate maps. `--order synteny` (GENESPACE-like, records sorted and flipped by the reference) or `--order assembly` (records grouped by the assembly's own chromosome assignment, FASTA order, released orientation; a chromosome-grouped contig view, not scaffolds); `--order both --compare LINE` draws the two side by side | **tool** |
| structural-confidence map | `scripts/g2f_structure_map.py`: synteny classes per region from a riparian run, joined to junction read calls, release edits and depth windows; read-supported reference differences vs. read-contradicted / unspanned candidates | **tool** |
| same-assembly duplicate candidates | `scripts/g2f_dup_pairs.py`: shared gene anchors, minimap2 identity / coverage (vs. each other and vs. the other haplotype), depth; rule-based reading | **tool** |
| reads across a specific breakpoint | `scripts/g2f_breakpoint_reads.py` | **tool** |
| per-junction call: SUPPORTED / CONTRADICTED / UNRESOLVED | junction evaluation over AnchorWave junctions plus HiFi and ONT spanning reads | **project script** |
| telomere identification and telomere-based integrity QC | `scripts/g2f_telomere.py` (see below) | **tool** |
| chromosome assignment | assignment from AnchorWave blocks | **project script** |

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

### Telomere integrity (`scripts/g2f_telomere.py`, **tool**)

Science: `config/telomere_science.yaml` (motif, thresholds). The motif is a species parameter: plants TTTAGGG,
vertebrates TTAGGG. Requires tidk ≥ 0.3.0; minimap2 and samtools are needed only for `tips` and `closure`. Each step
writes a new directory (built as `.partial`, renamed after validation) with a `VALIDATED.json`.

```
g2f_telomere.py check-env
g2f_telomere.py scan      --fasta UNIT.contigs.fa.gz   --label UNIT.contigs   --out scan/UNIT.contigs
g2f_telomere.py scan      --fasta UNIT.scaffolds.fa.gz --label UNIT.scaffolds --out scan/UNIT.scaffolds
g2f_telomere.py integrate --contigs-scan scan/UNIT.contigs --scaffolds-scan scan/UNIT.scaffolds --agp UNIT.agp \
                          --label UNIT --out integrate/UNIT [--join-table joins.tsv] [--split-ends splits.tsv]
g2f_telomere.py tips      --integrated integrate/UNIT --contigs C.fa.gz --scaffolds S.fa.gz --reference T2T.fa --out tips/UNIT
g2f_telomere.py closure   --integrated integrate/UNIT --alignments hifi_vs_contigs.cram --contigs C.fa.gz --out closure/UNIT
g2f_telomere.py plot      --integrated integrate/UNIT_A integrate/UNIT_B ... --out figures/
```

| step | what it does | gate |
|---|---|---|
| `scan` | `tidk ends` (primary call + sensitivity settings), `tidk search` (repeat counts per window), and a position-level array caller (every array, both strands, with coordinates) | tidk reported every record; the array caller's exact-hit totals equal `tidk search` up to hits straddling window boundaries |
| `integrate` | classes every contig and scaffold end (TEL_KB ≥ 1 kb, TEL_SHORT 200–999 bp, SUBTERMINAL, WRONG_STRAND, ABSENT); lifts every contig array ≥ 200 bp through the AGP; lists chromosome ends, telomeres that face an internal join, interstitial arrays, unplaced telomeric contigs | every lifted array is found at the identical position in the independent scaffold scan, and scaffold and contig end classes agree; otherwise it refuses (stale scan or wrong AGP) |
| `tips` | aligns each telomere-bearing contig end and each non-telomeric chromosome end (outer 1 Mb, end first) to a reference; calls REF_START / REF_END / INTERIOR_or_INVERTED in an outer (telomere-adjacent) and an inner (mostly unique) window | descriptive; prefer the inner window when the two disagree (shared subtelomeric repeats) |
| `closure` | for telomeric contig ends: reads that reach the end, reads clipped beyond it, and whether the clipped sequence is itself telomeric | the alignment reference must match the contigs (names, lengths, MD5) |
| `compare` | the same assembly before and after a scaffold revision: chromosome-end classes and kb-scale internal arrays | exit 3 if an end regressed or a new internal array appeared |

**Reading the results**
- **Counting telomeres:** count telomere-positive chromosome ends at scaffold level. Report contig-level telomeres
  that scaffolding displaced separately.
- **Expected background:** short interstitial runs (tens to ~150 bp) occur in complete genomes. Calibrate with a
  telomere-to-telomere assembly of the species before interpreting them.
- **A telomere at an internal join** means one of three things: the tip contig is reversed, sequence was appended
  beyond a real chromosome end, or a chromosome tip is fused inside a contig. Decide with `tips`, `closure`, depth
  and the join evidence. Arrays alone never decide.
- **Reference-guided scaffolds** inherit the reference's tip orientation. Where the reference itself carries an
  interstitial telomere or lacks a terminal one, check the tip against a telomere-to-telomere assembly.
- **What it does not show:** a telomere at each end does not show that a chromosome is structurally correct.

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
  `scripts/g2f_assembly_compare.py` (**tool**) builds that table from existing QC files: contiguity, NG50, size
  classes, gene-anchor statistics, BUSCO and Merqury, one row per FASTA. Compare haplotypes with haplotypes and
  primaries with primaries.
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
| remove foreign-contaminant and organellar unplaced records from a release (screens, candidate plan, apply, verify) | `scripts/g2f_decontam.py`, `scripts/g2f_release_exclude.py` ([Stage 8d](#stage-8d--removing-contaminant-and-organellar-records-tool)) | **tool** |
| short, reversible sequence IDs for annotation tools | `scripts/g2f_seqid_alias.py` | **tool** |

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

## Stage 8b — Reference-guided pseudomolecules (**tool**, optional derived product)

Entry: `scripts/g2f_allmaps.py` (`check-env`, `plan`, `run`, `status`; configuration `config/allmaps.example.yaml`,
science `config/allmaps_science.yaml`). It orders and orients released contigs into chromosome pseudomolecules from
sequence-anchor markers with a known chromosome (LG) and map position, using HISAT2 and ALLMAPS (jcvi), both from
pinned containers.

| stage | what it does | gate |
|---|---|---|
| `markers` | validates the marker FASTA and metadata (IDs match, lengths, duplicates); maps each distinct sequence once; keeps a marker-ID → unique-marker → map-position table | every marker ID accounted for; identical sequences must carry one LG and one position |
| `index:S` | `hisat2-build` of one assembly in its own stage directory | the index lists exactly the FASTA's records and lengths; optional expected sha256 |
| `map:S` | maps the unique markers; keeps a marker only if its primary alignment is unique (`NH` ≤ 1, MAPQ ≥ 60, no `ZS` second-best), unspliced, unclipped and within the edit-distance limit | counts of every rejection class recorded |
| `scaffold:S` | per-contig support from **unique** markers; contig LG conflicts; deterministic, evenly spaced ALLMAPS input; ALLMAPS `merge` + `path`; independent AGP reconstruction | every source contig used exactly once, entire, in the stated orientation; output = source + gaps, base for base |
| `scaffold_refpos:S` (optional; config `reference_positions`) | the same as `scaffold:S`, but each marker's map position is replaced by the position of its unique hit in a reference sample of the same run (e.g. the current reference assembly); markers whose reference hit is not unique or not on the chromosome of their LG are dropped | as `scaffold:S`; needs `map:REF` |

`scaffold_refpos` exists for marker maps whose positions are coordinates of an **older reference version**: blocks
inverted between reference versions would otherwise orient contigs as the older reference did. Both variants are
reference-derived; which one a product uses is a scientific decision recorded with it.

`placement.tsv` reports, per placed contig, a robust map interval (5th–95th percentile of its markers' positions) and
two review flags, neither of which changes the placement:
- `OVERLAPS_LARGER`: at least half of the contig's interval lies inside the intervals of longer contigs placed on the
  same pseudomolecule. This is typical of a redundant copy (an unpurged allele or haplotig) placed in tandem, which
  inflates the pseudomolecule.
- `ORDER_AMBIGUOUS`: the intervals of two adjacent contigs overlap by more than half of the narrower one, so their
  relative order is weakly determined.

**Rules**
1. **The product is derived and labelled as such.** The gap-free contig FASTA stays the primary product and the
   annotation input. Pseudomolecules are a separate, versioned representation with an AGP and a chain file, so
   coordinates can be lifted in both directions.
2. **Reference-guided ordering is not validation.**
   - Marker map positions are expressed in reference coordinates. Agreement with other reference-based analyses
     (gene-anchor synteny) is consistency, not independent evidence.
   - A contig placed on a chromosome is not thereby structurally correct.
   - Gaps between contigs have a fixed length and no measured size.
3. **Contigs are never split, trimmed or edited here.**
   - A contig whose markers come from several LGs is placed by its majority LG.
   - It is flagged (`MIXED_LG`, `CONFLICT_BLOCK`) and goes to junction evaluation (stage 5) for a human decision.
4. **Support is counted in unique markers and distinct map positions,** never in redundant marker records. Markers
   that share one map position carry no order information.
5. **One assembly, one index.** An index is built inside the assembly's stage directory and checked against that
   FASTA; it is never found by searching a shared directory.
6. **Haplotype assemblies are scaffolded separately,** one haploid FASTA per run. Two haplotypes are never merged
   into one pseudomolecule set.

**Differences from the NAM-style reference implementation** that this stage replaces:

| reference implementation | this stage | why |
|---|---|---|
| uniqueness by removing records with `XA`/`SA` tags | `NH`, MAPQ, `ZS` and secondary records | HISAT2 does not write `XA`/`SA`; that filter removes nothing |
| HISAT2 defaults (spliced alignment, up to 5 alignments reported) | `--no-spliced-alignment`, `-k 10`, primary record only | markers are genomic windows; repeats must be visible to be rejected |
| every marker record mapped; rows collapsed afterwards by `sort \| uniq` once the ID is cut | each distinct sequence mapped once | identical records are one piece of evidence |
| up to 100 random markers per contig (`sort -R`) | up to N markers evenly spaced along the contig | deterministic; keeps positional coverage |
| contig needs > 20 rows (secondary alignments and repeat copies included) | thresholds on unique markers and distinct map positions | non-unique hits and shared map positions are not order evidence |
| non-majority-LG markers silently removed | removed from the ALLMAPS input but reported per contig, with conflict blocks | chimeric contigs must stay visible |
| index reused if any `*.ht2` exists below `./index` | index per assembly, checked against the FASTA | avoids using another assembly's index |
| outputs written to the current directory; the per-contig map appended with `>>` | `.partial` stage directories, `VALIDATED.json`, nothing overwritten | a rerun cannot duplicate markers; resumable and auditable |
| map name taken from the output name; no check of the merged map | fixed map file name; merged-marker count checked | `allmaps merge` silently skips lines it cannot parse; marker names are split on `-` (jcvi 1.6.6 sanitises map names, other versions not checked) |
| no output check | AGP reconstruction compared with the output FASTA | proves no sequence was lost, duplicated or altered |


### Converting coordinates between contigs and scaffolds (`scripts/g2f_agp_lift.py`, **tool**)

The AGP of a reference-guided product maps every contig to its scaffold position, so features annotated on one
representation can be converted to the other:

```
python3 scripts/g2f_agp_lift.py gff --agp UNIT.agp --to object    -i genes.contigs.gff3   -o genes.scaffolds.gff3 --unmapped u.tsv
python3 scripts/g2f_agp_lift.py gff --agp UNIT.agp --to component -i genes.scaffolds.gff3 -o genes.contigs.gff3   --unmapped u.tsv
python3 scripts/g2f_agp_lift.py check --agp UNIT.agp --component-fasta contigs.fa --object-fasta scaffolds.fa --chain UNIT.chain
```

- `−` components mirror coordinates and flip the strand; IDs, attributes and CDS phase are unchanged.
- Features that overlap a gap or two components have no contig equivalent. They are written to `--unmapped` with the
  reason, and a gene hierarchy is kept or dropped as a whole.
- `chain` and `gaps` write a UCSC chain and a gap BED from the AGP.
- `check` verifies the AGP reconstruction and tests random intervals in both directions by sequence comparison.

## Stage 8c — Scaffold revision (**tool**, AGP level)

Entry: `scripts/g2f_agp_revise.py` (`orient-evidence`, `check-plan`, `apply`) and the independent verifier
`scripts/g2f_agp_revise_verify.py`. A revision changes placements only. Contigs are never edited, and every contig
stays in the FASTA.

1. **Evidence.**
   - Telomere integrity (above).
   - For tip orientation, `orient-evidence`: it aligns the leading or trailing components of chosen chromosome
     objects to a telomere-to-telomere reference. It proposes `REVERSE_ORIENTATION` for a single reversed tip
     component, or `MANUAL_REVIEW_BLOCK` when several components would change order.
2. **Plan.** A TSV with the columns `action, object, component, old_orientation, new_orientation, new_status,
   evidence_class, evidence` (an optional `unit` column). The actions are:
   - `REVERSE_ORIENTATION`;
   - `UNPLACE_TERMINAL_EXTENSION`: trailing components and the gap before each become unplaced records.

   Deciding the plan is a human decision. Label the evidence class honestly: telomere-supported versus synteny-only.
3. **`check-plan`, then `apply`.**
   - `apply` writes the new AGP, the scaffold FASTA rebuilt from the contigs, chain, gap BED, contig-to-scaffold map,
     line-level change manifest and `BUILD.json`.
   - With `--join-table`, removed joins are listed separately, and joins next to a reversed component are marked as
     not re-tested.
4. **Verify** with `g2f_agp_revise_verify.py`. It must PASS before the output is used:
   - the AGP is valid and reconstructs the FASTA;
   - the old → new AGP difference equals the plan exactly;
   - unchanged records are byte-for-byte the same sequence;
   - no sequence is lost;
   - the derived tables agree;
   - random intervals lift correctly.
5. **Re-scan** the revised scaffolds and run `g2f_telomere.py compare` against the old scan: no regressed end, no new
   internal array.
6. **Release.** A new versioned release keeps the previous AGPs and the change manifest. Coordinates of a feature
   convert old scaffold → contigs → new scaffold with `g2f_agp_lift.py`.

## Stage 8d — Removing contaminant and organellar records (**tool**)

Entry: `scripts/g2f_decontam.py` (screens and a candidate plan), `scripts/g2f_release_exclude.py` (applies the plan)
and the independent verifier `scripts/g2f_release_exclude_verify.py`. Science: `config/decontam_science.yaml`.

Whole records are removed only when they are unplaced singletons: the record is the sole, whole, + component of an AGP
object of its own name. Removing anything placed in a scaffold would be a structural edit (Stage 6 / 8c). Organellar
sequence inside nuclear contigs (NUMTs and NUPTs) is authentic nuclear sequence and stays.

```bash
g2f_decontam.py check-env --fcs-dir /path/to/fcs
# foreign sequence: FCS-adaptor + FCS-GX (one assembly at a time; GX wants its database in shared memory)
g2f_decontam.py fcs       --fasta UNIT.contigs.fa.gz --agp UNIT.agp --label UNIT --tax-id TAXID \
                          --gx-db /dev/shm/gxdb/all --fcs-dir /path/to/fcs --threads 48 --out fcs/UNIT
# the species' organelles: windows vs plastid + mitochondrial genomes; nuclear test against nuclear chromosomes
g2f_decontam.py organelle --fasta UNIT.contigs.fa.gz --agp UNIT.agp --label UNIT --organelles organelles.fa \
                          --plastid-ids PLASTID_ID --nuclear nuclear_chromosomes.fa --threads 16 --out organelle/UNIT
# candidate plan for all assemblies
g2f_decontam.py plan      --fcs fcs/* --organelle organelle/* --out plan
```

1. **Screen.**
   - `fcs` runs NCBI FCS-adaptor (`--euk`) and FCS-GX (`screen genome` with the species' NCBI taxonomy ID). Download
     and check the GX database once with FCS's own `fcs.py db get` and `db check` (about 500 GB), and record its build
     date. `fcs-summarize` reads FCS output run elsewhere. Calls are the tools' own: GX `EXCLUDE` / `TRIM` / `FIX`
     count as foreign; `REVIEW` and `INFO` are ambiguous and only reported; adaptor trims are reported per range.
   - `organelle` cuts each unplaced record into 5-kb windows and aligns each window to the organelle references
     (minimap2 `asm5`, secondary alignments kept). Coverage is the union of blocks at ≥ 90 % identity. Aligning windows
     separately avoids a chaining artefact in which a short organellar stretch inside a long contig is dropped. Records
     with ≥ 50 % coverage get a nuclear test: their non-organelle stretches are aligned to the nuclear reference.
     - **CONFIDENT_ORGANELLE**: coverage ≥ 0.90 and < 5 kb of nuclear-aligned non-organelle sequence (plastid or
       mitochondrion by the larger share, with `--plastid-ids`).
     - **ORGANELLE_DOMINANT_MIXED**: coverage ≥ 0.50 otherwise. Never excluded automatically.
     - Use references for every cytotype of the species: mitochondrial genomes differ between cytoplasms.
2. **Plan.** `plan` writes `EXCLUSION_PLAN.candidate.tsv` (`unit, record, reason_class, evidence`) with only
   whole-record calls on unplaced singletons: `FOREIGN_CONTAMINANT` (GX `EXCLUDE`) and `ORGANELLE_PLASTID` /
   `ORGANELLE_MITOCHONDRION` (`CONFIDENT_ORGANELLE`). Everything else goes to `REVIEW.tsv`: partial calls, ambiguous GX
   calls, mixed organelle records, a record with two different whole-record calls, and calls on placed contigs.
3. **Decide.** The plan is a candidate. A human reviews it with `REVIEW.tsv` and decides which records to exclude. An
   adaptor call inside a contig may be native sequence; check it (e.g. read support across the call) before masking.
4. **Apply.** `g2f_release_exclude.py check-plan`, then `apply`. Kept records are copied byte for byte, and the AGP,
   chain, contig-to-scaffold map and gap BED lose exactly the excluded records. `EXCLUDED_RECORDS.tsv` lists each
   removed ID with its length and sequence hash, so annotation features on those IDs can be found and dropped.
5. **Verify** with `g2f_release_exclude_verify.py`, which shares no code with the builder. It must PASS before release.
6. **Re-screen** the result: the organelle screen of the remaining unplaced records should find no
   `CONFIDENT_ORGANELLE` record.
7. **Release.** A new versioned release with the plan, the screens and the verifier output as evidence.

### Short sequence IDs for annotation tools (`scripts/g2f_seqid_alias.py`, **tool**)

Some annotation tools (EDTA, for example) truncate sequence IDs longer than 13 characters or rewrite punctuation.
`g2f_seqid_alias.py make` writes an alias FASTA with byte-identical sequences and a one-to-one map;
`g2f_seqid_alias_verify.py` checks it independently; `translate` converts GFF3, BED, TSV and RepeatMasker `.out`
results back to the original IDs and stops on an unknown ID.

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
