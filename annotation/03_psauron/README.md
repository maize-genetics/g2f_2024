# 03 · PSAURON gene-model scores (whole genome)

Scores every Helixer gene model with PSAURON, which predicts from the CDS sequence alone how
likely it is to encode a real protein, and adds the score to the GFF3.

| Script | What it does |
|---|---|
| `psauron_run.sbatch` | runs PSAURON on the CDS of each genome, one Slurm array task per genome |
| `psauron_add_to_gff3.py` | adds `psauron_score` and `psauron_is_protein` to the gene and mRNA lines of the whole-genome GFF3 and writes a summary per genome |

## Settings

- **PSAURON v1.1.3**, CPU mode (`-c`), on the CDS from step 02 (stop codon included).
- PSAURON scores each CDS in its reading frame (`in_frame_score`, 0–1) and in the five other
  frames. A model is called protein-coding (`psauron_is_protein=True`) when the in-frame score is
  ≥ 0.5 and the mean out-of-frame score is < 0.5.
- `psauron_score` in the GFF3 is the in-frame score (4 decimals). CDS too short for PSAURON to
  score get `NA`.
- Set `PSAURON` in `config.sh` to however PSAURON is run on your system (installed, or through its
  container). Each genome used 8 CPUs and 32 GB memory.

## Run

From the repository folder, after step 02:

```bash
source config.sh && mkdir -p logs
N=$(wc -l < "$PROJECT_DIR/units.tsv")
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 03_psauron/psauron_run.sbatch       # test one genome
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 03_psauron/psauron_run.sbatch    # the rest
# when all have finished:
python3 03_psauron/psauron_add_to_gff3.py --project "$PROJECT_DIR"
```

## Inputs and outputs

| | Path (under `PROJECT_DIR`) |
|---|---|
| CDS (input) | `raw/<unit>/<unit>.scaffolds.helixer.cds.fa` |
| PSAURON scores | `psauron/<unit>.helixer.cds.psauron.csv` (one row per transcript; the genome's overall PSAURON score is in the first lines) |
| GFF3 with scores | `raw/<unit>/<unit>.scaffolds.helixer.psauron.gff3` |
| summary | `psauron/psauron_summary.tsv` (genes, scored, not scored, % protein-coding, mean and median score) |

Example gene and mRNA lines after the scores are added:

```
chr1  Helixer  gene  ...  ID=GENOTYPE1_chr1_000001;psauron_score=0.9731;psauron_is_protein=True
chr1  Helixer  mRNA  ...  ID=GENOTYPE1_chr1_000001.1;Parent=GENOTYPE1_chr1_000001;psauron_score=0.9731;psauron_is_protein=True
```
