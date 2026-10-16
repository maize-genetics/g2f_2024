# Assembly structural validation (`scripts/g2f_asmqc.py`)

Evidence-based structural QC of genome assemblies, designed for comparing
assemblies of the same reads made with different settings (for example hifiasm
`-l0` vs `-l3` purging). It goes beyond size, N50, BUSCO and k-mer metrics
and asks, for each assembly:

* **Is every part supported by the reads?** HiFi reads are aligned back to
  the assembly (minimap2 `map-hifi`). Per window: depth, unique (MAPQ ≥ 10)
  depth, read continuity (unique reads extending ≥ 1 kb on both sides of a
  position), mixed sites (NucFreq-style: positions where 20-80 % of unique
  reads disagree with the assembly), likely consensus errors, and clusters of
  clipped / split reads with the location where the split reads continue.
* **Does it agree with a reference, and where it does not, do the reads?**
  Whole-genome alignment (minimap2 `asm20`) gives contig-to-chromosome
  assignment, synteny segments, junctions (inter-chromosomal, inversion,
  large jump), reference coverage and dot plots. A junction is a *suspected
  misassembly* only when read continuity also breaks there; a junction the
  reads span continuously is reported as a *likely biological structural
  variant* (or reference difference).
* **What does the assembly graph say?** From hifiasm `*.noseq.gfa` files:
  unitig topology (tips, branching ends, simple bubbles, components), each
  contig's path through the unitig graph (via the read names in the A lines),
  transitions between unlinked unitigs, and local graph neighbourhoods
  (GFA + SVG/PDF) around the top candidates.
* **How do two assemblies of the same reads differ?** Assembly A is aligned to
  B (minimap2 `asm5`, base level) and to B's other outputs (e.g. hifiasm
  hap1/hap2). Every base of A is classed as shared, redundant (an extra copy
  of a B region, near-identical or divergent), only in B's other outputs, or
  unique to A, and cross-tabulated with read depth in both assemblies. A
  read-retention check uses the hifiasm contig read lists.

No result from a single signal is reported as a confirmed error. Every
candidate carries its evidence, a confidence level and the number of
independent evidence types.

## Candidate categories

| category | meaning |
|---|---|
| `suspected_misassembly` | clipped/split reads at one position with (almost) no read spanning it, or a reference junction where read continuity breaks |
| `likely_biological_sv` | reference junction spanned continuously by reads, or a breakpoint spanned by about half the reads (heterozygous SV) |
| `collapsed_repeat_candidate` | depth ≥ 1.4× single-copy; high confidence when ≥ 1.75× with enriched mixed sites over ≥ 20 kb |
| `false_duplication_candidate` | half depth with mostly multi-mapping reads: near-identical copies splitting the reads |
| `retained_haplotypic_copy` | half depth with unique reads: one allele of a heterozygous region (expected in unpurged assemblies of heterozygous samples) |
| `low_support` | little or no read support |
| `unresolved` | conflicting signals, or a repeat longer than the reads at the site |

Depth classes are relative to the expected single-copy depth *c* (the
depth-histogram mode, cross-checked against the k-mer peak, or set
explicitly): `none` < 0.15c ≤ `low` < 0.30c ≤ `half` < 0.70c ≤ `normal` <
1.40c ≤ `dup` < 2.50c ≤ `multi`.

Graph depth classes use the hifiasm `rd` tag relative to hifiasm's own
homozygous peak (`peak_hom` in its log), not the `rd` mode: in a sample that is
heterozygous over much of the genome, the mode is the one-allele depth.

## Requirements

minimap2 (≥ 2.24), samtools (≥ 1.15), seqkit; Python ≥ 3.9 with numpy,
pysam, matplotlib, networkx and PyYAML. No other tools are needed.
`check-env` verifies all of them and every input path.

## Configuration

* `config/asmqc.example.yaml`: run configuration template (samples,
  assemblies, references, pairs, resources). Copy it to `config/local/`.
* `config/asmqc_science.yaml`: scientific parameters, version-controlled.
* `config/tools/asmqc_biohpc.yaml`: optional fixed tool paths for one site.

Adding an assembly or a sample is a configuration change only.

## Running

```bash
python3 scripts/g2f_asmqc.py -c CONFIG check-env   # tools, modules, inputs, free space
python3 scripts/g2f_asmqc.py -c CONFIG plan        # states and resource estimates; runs nothing
python3 scripts/g2f_asmqc.py -c CONFIG run         # all planned stages, one at a time
python3 scripts/g2f_asmqc.py -c CONFIG status
python3 scripts/g2f_asmqc.py -c CONFIG report      # new report directory per sample
```

`--stages` and `--samples` restrict any command. `bash scripts/run_asmqc.sh
CONFIG` does check-env → plan → run → report. Run it detached so the session
ends with the run; for example with a helper that starts a screen session
whose only process is the job:

```bash
run_detached asmqc /path/to/logs -- bash scripts/run_asmqc.sh config/local/asmqc.yaml
```

Stages, in execution order (cheap analyses first): `prepare`, `graph`,
`ref_index`, `ref_align`, `ref_analysis`, `pair_align`, `map_hifi`,
`read_support`, `regions`, `evidence`, `pair_analysis`.

## States, resumption and safety

Each stage runs in `<dir>.partial` and becomes `<dir>` only after its outputs
validate (for example, primary alignment records must equal the configured
read count). `<dir>/VALIDATED.json` records the command, tool versions,
inputs (size, mtime), parameters, code digest, host, timing and validation.

* **validated** stages are skipped.
* **stale** stages (inputs, parameters or code changed) are never overwritten;
  `--rerun-stale` keeps the old result as `<dir>.previous-<stamp>`.
* **failed** attempts keep `FAILED.json`; the next attempt moves them to
  `<dir>.attempt-<stamp>`. Nothing is deleted.
* A resource gate refuses to start a stage whose disk / RAM estimate does not
  fit (with a reserve).
* Only one runner per output root (lock file).
* External commands can be pinned to a CPU list and niced.
* Input assemblies are only read: FASTA files are symlinked into the output
  tree and indexed there.

## Outputs (per assembly unless noted)

| stage | key files |
|---|---|
| prepare | `stats.json` (size, contigs, N50/L50/N90, assembler log values), `contigs.tsv` |
| map_hifi | sorted CRAM (no base qualities) + index, `flagstat.json` |
| read_support | `windows.tsv.gz`, `clip_clusters.tsv.gz`, `depth_hist.tsv` |
| regions | `regions.tsv`, `breakpoints.tsv`, `window_classes.tsv.gz`, depth histogram / class plots |
| graph | `summary.json`, `contig_utg_path.tsv`, `path_transitions.tsv`, `utg_segments.tsv` |
| ref_analysis | `contig_assignment.tsv`, `synteny_segments.tsv`, `junctions.tsv`, `reference_coverage.tsv`, `dotplot_genome.svg/pdf` |
| evidence | `candidates.tsv` (ranked), `neighbourhoods/*.gfa/svg/pdf` |
| pair_analysis (per pair) | `summary.json` (class bp, depth cross-tab), `redundant_copies.tsv`, `a_unaligned.tsv`, `read_retention_*` |
| report (per sample) | `metrics.tsv` (long format with source), `report.md`, `report.html`, `figures/` |

## Interpreting with care

* Maize and other plant genomes differ from any reference by large
  presence/absence variation, inversions and transposon insertions.
* Reads shorter than a repeat cannot resolve it; such sites are unresolved,
  not errors.
* Depth and continuity can be biased by sequence composition.
* BUSCO duplication alone does not show haplotypic duplication in genomes with
  a history of whole-genome duplication; use the read-depth and pairwise
  classes instead.

## Tests

`python3 -m unittest tests.test_asmqc` builds a 500 kb synthetic diploid genome
with a planted misjoin, a collapsed repeat, a retained divergent allele and a
near-identical false duplicate, runs the real tools end to end (~30 s), and
checks that each is found, and that skip, stale, failure and read-count
validation behave as described.
