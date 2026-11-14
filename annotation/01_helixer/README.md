# 01 · Helixer gene prediction (whole genome)

Predicts protein-coding gene models on each whole genome assembly with Helixer, a deep-learning
gene predictor that needs no RNA-seq or protein evidence.

| Script | What it does |
|---|---|
| `helixer_prep_inputs.sbatch` | copies (or decompresses) each assembly listed in `ASSEMBLIES_TSV` (unit, source, FASTA), records its number of sequences, total bp and sha256, and writes `units.tsv` |
| `helixer_array.sbatch` | runs Helixer on each assembly, one Slurm array task (one GPU) per genome |

## Settings

- **Helixer v0.3.3**, `--lineage land_plant`, which uses the default model
  `land_plant_v0.3_a_0080` (md5 `f0e00efcbea83c66b69258d11119a691`) and the default window size.
  `--species` is set to the unit name, so gene IDs start with it.
- Set `HELIXER` in `config.sh` to however Helixer is run on your system (installed, or through its
  container), and download the `land_plant` model first (`fetch_helixer_models.py`).
- Helixer writes large temporary files; if the default location is too small, add
  `--temporary-dir <folder>` to the Helixer command in `helixer_array.sbatch`.
- Each genome was run with one A100 GPU, 8 CPUs and 128 GB memory (24 h limit).

## Run

From the repository folder, after filling in `config.sh`:

```bash
source config.sh && mkdir -p logs
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 01_helixer/helixer_prep_inputs.sbatch
# when it reports "All inputs ready":
N=$(wc -l < "$PROJECT_DIR/units.tsv")
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 01_helixer/helixer_array.sbatch        # test one genome
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N%4 01_helixer/helixer_array.sbatch   # the rest, 4 at a time
```

## Inputs and outputs

| | Path (under `PROJECT_DIR`) |
|---|---|
| input table | `ASSEMBLIES_TSV` (see `../assemblies_example.tsv`) |
| prepared assemblies | `inputs/<unit>.scaffolds.fa`, checks in `inputs/prep_check.tsv` |
| list of units | `units.tsv` (used by every later step) |
| gene models | `raw/<unit>/<unit>.scaffolds.helixer.raw.gff3` |
| run record | `raw/<unit>/<unit>.run_info.txt` (Helixer version, GPU, model md5, parameters) |
| feature counts | `raw/<unit>/<unit>.summary_gff3.txt` |

Helixer predicts one transcript per gene, with exons, CDS and 5′/3′ UTRs. The GFF3 is written under
a temporary name and renamed only when Helixer finishes, so a crashed run never leaves a file that
looks complete; a resubmitted task skips genomes that are already done.
