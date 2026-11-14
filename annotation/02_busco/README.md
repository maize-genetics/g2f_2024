# 02 · BUSCO gene-space completeness (whole genome)

Extracts the protein and CDS sequences of every Helixer gene model and measures annotation
completeness with BUSCO against the Poales lineage.

| Script | What it does |
|---|---|
| `busco_extract_seqs.sbatch` | writes protein and CDS FASTA files from each Helixer GFF3 with gffread, removes the terminal stop from proteins, and checks that every gene has one protein and one CDS with matching IDs |
| `busco_run.sbatch` | runs BUSCO on the proteins of each genome, one Slurm array task per genome |

## Settings

- **gffread** (`-y` proteins, `-x` CDS) on the whole-genome GFF3 from step 01.
- **BUSCO v6.1.0**, protein mode (`-m proteins`), lineage **poales_odb12** (6,282 BUSCOs),
  `--offline`.
- Download the lineage once before running (for example `busco --download poales_odb12`) and set
  `BUSCO_LINEAGE` in `config.sh` to the `poales_odb12` folder. Set `GFFREAD` and `BUSCO` if the
  tools are not on your PATH.
- Each BUSCO run used 16 CPUs and 32 GB memory.

## Run

From the repository folder, after step 01 has finished:

```bash
source config.sh && mkdir -p logs
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" 02_busco/busco_extract_seqs.sbatch
# when extract_summary.tsv shows OK for every unit:
N=$(wc -l < "$PROJECT_DIR/units.tsv")
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 02_busco/busco_run.sbatch       # test one genome
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 02_busco/busco_run.sbatch    # the rest
```

## Inputs and outputs

| | Path (under `PROJECT_DIR`) |
|---|---|
| gene models (input) | `raw/<unit>/<unit>.scaffolds.helixer.raw.gff3` |
| proteins | `raw/<unit>/<unit>.scaffolds.helixer.proteins.fa` (terminal stop removed) |
| CDS | `raw/<unit>/<unit>.scaffolds.helixer.cds.fa` (stop codon included) |
| sequence checks | `busco/extract_summary.tsv` |
| BUSCO results | `busco/runs/<unit>/short_summary*.txt` |

`extract_summary.tsv` reports, per genome, the number of genes, proteins and CDS, the % of
proteins starting with methionine, proteins with internal stop codons, CDS whose length is not a
multiple of three, and the median protein length. The extraction job stops with an error if any
genome's counts or IDs do not match. A resubmitted BUSCO task skips genomes that already have a
result.
