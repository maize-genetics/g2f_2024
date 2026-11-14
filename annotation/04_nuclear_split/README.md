# 04 · Nuclear genome split

Labels every scaffold of each assembly as nuclear, plastid or mitochondrion, then writes the nuclear
genome and its gene models. Steps 05 (B73 synteny) and 06 (TE annotation) use
these nuclear files.

| Script | What it does |
|---|---|
| `organelle_screen.sbatch` | aligns each assembly to the maize plastid and mitochondrial genomes with minimap2 and labels every scaffold (runs `classify_scaffolds.py`) |
| `classify_scaffolds.py` | the labelling rule, one assembly at a time |
| `split_nuclear.py` | writes the nuclear genome FASTA, GFF3, CDS and proteins of each assembly, and a table labelling every gene by where it sits |

## Settings

- **Organelle reference:** maize plastid (NC_001666) and mitochondrion (NC_007982) genomes in one
  FASTA (`ORGANELLE_FASTA` in `config.sh`).
- **minimap2** `-x asm5 -c`. Alignments with ≥ 95% identity and ≥ 500 bp are kept; the fraction
  of each scaffold they cover (overlaps counted once) decides its label:

  | Covered by | Label |
  |---|---|
  | plastid ≥ 80% | plastid |
  | mitochondrion ≥ 80% | mitochondrion |
  | plastid + mitochondrion ≥ 30% | review |
  | otherwise | nuclear |

  A `review` scaffold needs a decision in `scaffold_classes/review_overrides.tsv` (unit, scaffold,
  final_class, confidence, reason); rerun the screen afterwards (the alignments are reused).
- **Nuclear genome:** the 10 chromosome scaffolds, found by name (`chr1`…`chr10`, `chr01`,
  `<unit>_chr01`, and split chromosomes such as `chr4_part1`). Each assembly must have all 10, and
  each must be labelled nuclear. Unplaced scaffolds are left out.
- Set `MINIMAP2` in `config.sh` if minimap2 is not on your PATH.

## Run

From the repository folder, after step 03:

```bash
source config.sh && mkdir -p logs
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 04_nuclear_split/organelle_screen.sbatch
# when the screen has finished and no scaffolds are left to decide, write the nuclear files:
python3 04_nuclear_split/split_nuclear.py --project "$PROJECT_DIR"
```

## Inputs and outputs

| | Path (under `PROJECT_DIR`) |
|---|---|
| assemblies (input) | `inputs/<unit>.scaffolds.fa` |
| gene models (input) | `raw/<unit>/<unit>.scaffolds.helixer.psauron.gff3`, `.cds.fa`, `.proteins.fa` |
| scaffold labels | `scaffold_classes/<unit>.scaffold_classification.tsv`, summary in `scaffold_classes/scaffold_summary.tsv` |
| nuclear genome | `nuclear_genome/<unit>/<unit>.nuclear.fa` |
| nuclear gene models | `nuclear_genome/<unit>/<unit>.nuclear.helixer.gff3`, `.cds.fa`, `.proteins.fa` |
| gene labels | `raw/<unit>/<unit>.gene_classes.tsv` (every gene: scaffold, position, group) |
| gene summary | `nuclear_genome/gene_classes_summary.tsv` (genes per group, per unit) |

Gene groups in `gene_classes.tsv`: `nuclear` (on the chromosomes), `unplaced_nuclear`, `plastid`,
`mitochondrion`, `review`.

**Note:** the nuclear GFF3 files include the PSAURON scores (`psauron_score` and
`psauron_is_protein` on every gene and mRNA line), because the split starts from the GFF3 written
in step 03.
