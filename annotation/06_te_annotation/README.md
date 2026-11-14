# 06 · Transposable element annotation (nuclear genome)

Annotates transposable elements (TEs) on the nuclear genome of each assembly with EDTA and
RepeatMasker, and measures assembly continuity in repeats with the LTR Assembly Index (LAI).

| Script | What it does |
|---|---|
| `te_prep_inputs.py` | copies each nuclear genome (step 04) into the TE folder with short sequence names, and writes the name map |
| `te_submit.sh` | submits the TE jobs of each genome as one chain (each job waits for the one before) |
| `te_edta_ltr.sbatch`, `te_edta_tir.sbatch`, `te_edta_helitron.sbatch`, `te_edta_line.sbatch`, `te_edta_sine.sbatch` | structural detection with `EDTA_raw.pl`, one job per TE type, run at the same time |
| `te_edta_combine.sbatch` | combines the five raw outputs into one TE library with `EDTA.pl` |
| `te_repeatmasker.sbatch` | annotates every TE copy with RepeatMasker using that library, and writes a soft-masked genome |
| `te_restore_ids.sbatch`, `te_restore_ids.py` | put the real sequence names back and collect the final files |
| `te_lai.sbatch`, `te_lai_summary.py` | LAI of each genome, and one table of all genomes |

```
te_edta_ltr      \
te_edta_tir       |
te_edta_helitron  +--> te_edta_combine --> te_repeatmasker --> te_restore_ids --> te_lai
te_edta_line      |
te_edta_sine     /
```

## Settings

- **Input:** the nuclear genome of each assembly. Sequences are renamed `s00001`, `s00002`, …
  for the run, because EDTA 2.3.0 re-codes names longer than 8 characters while EDTA 2.2.2 (TIR
  step) does not; `te_restore_ids.py` puts the real names back.
- **Structural detection:** `EDTA_raw.pl --species Maize --type <type>`, one job per type:

  | Job | TE type | EDTA (container) | Main tool inside |
  |---|---|---|---|
  | `te_edta_ltr` | LTR retrotransposons | 2.3.0 | LTRharvest (GenomeTools 1.6.5), LTR_FINDER 1.07, LTR_retriever 3.0.4 |
  | `te_edta_tir` | TIR transposons and MITEs | **2.2.2** (TIR detection fails in 2.3.0) | TIR-Learner 3.0.3 |
  | `te_edta_helitron` | Helitrons | 2.3.0 | HelitronScanner 1.0 |
  | `te_edta_line` | LINEs | 2.3.0 | EDTA LINE module |
  | `te_edta_sine` | SINEs | 2.3.0 | AnnoSINE2 2.0.9 |

  `--species Maize` adds known maize TIR elements to TIR-Learner; it is not a TE library for the
  whole annotation.
- **Library:** `EDTA.pl --species Maize --sensitive 1 --anno 0 --evaluate 0 --overwrite 0 --cds <CDS>`
  (EDTA 2.3.0). It filters false positives, runs RepeatModeler 2.0.7 (`--sensitive 1`) to find
  families missed by the structural tools, removes gene sequences using the nuclear Helixer CDS
  from step 04, removes redundancy (CD-HIT 4.8.1) and classifies families (TEsorter 1.5.1). No
  curated maize library is used.
- **Whole-genome annotation:** RepeatMasker 4.2.3 with the genome's own library,
  `-e ncbi -q -div 40 -cutoff 225 -gff`. A bare `DNA` class in the `.out` file is changed to
  `DNA/unknown` (LAI expects a subclass). The soft-masked FASTA is made from the RepeatMasker GFF
  with `bedtools maskfasta -soft`.
- **LAI:** the `LAI` program of LTR_retriever 3.0.4 (EDTA 2.3.0 image), from the intact LTR pass
  list and the RepeatMasker `.out`.
- **Containers:** set `EDTA230_SIF`, `EDTA222_SIF` and `REPEATMASKER_SIF` in `config.sh`
  (biocontainers `edta:2.3.0--hdfd78af_0`, `edta:2.2.2--hdfd78af_0`,
  `repeatmasker:4.2.3--pl5321hdfd78af_0`); apptainer must be available on the compute nodes.
- **Resources per genome:** LTR 48 CPUs / 300 GB, TIR 32 CPUs / 500 GB, Helitron 48 CPUs /
  300 GB, LINE 24 CPUs, SINE 12 CPUs, combine and RepeatMasker as set in their scripts; up to
  5 days each. LAI 16 CPUs / 64 GB, and can take more than a day per genome.

## Run

From the repository folder, after step 04:

```bash
source config.sh && mkdir -p logs
python3 06_te_annotation/te_prep_inputs.py --project "$PROJECT_DIR"
bash 06_te_annotation/te_submit.sh                 # every genome; or name some: te_submit.sh UNIT1 UNIT2
# when all te_restore_ids jobs have finished:
N=$(wc -l < "$PROJECT_DIR/units.tsv")
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1-$N 06_te_annotation/te_lai.sbatch
# when all LAI tasks have finished:
python3 06_te_annotation/te_lai_summary.py --project "$PROJECT_DIR"
```

`te_submit.sh` takes the account, partition and e-mail from `config.sh`. To resume a genome part
way, set `START` to `combine`, `repeatmasker` or `restore`, e.g.
`START=repeatmasker bash 06_te_annotation/te_submit.sh UNIT1`.

## Inputs and outputs

| | Path (under `PROJECT_DIR/te_annotation/`) |
|---|---|
| genomes with short names | `inputs/<unit>.nuclear.fa`, name map `inputs/<unit>.nuclear.idmap.tsv` |
| raw detection | `runs/<unit>/edta_raw/<unit>.<type>.tar` |
| EDTA combine | `runs/<unit>/edta/<unit>.tar` (all EDTA files) |
| RepeatMasker | `runs/<unit>/repeatmasker/<unit>.tar` |
| job logs | `runs/<unit>/logs/` |
| **final files** | `runs/<unit>/results/` (below) |
| LAI | `runs/<unit>/lai/`, all genomes in `te_LAI.tsv` |

Final files in `runs/<unit>/results/`, with the real sequence names:

| File | Contents |
|---|---|
| `<unit>.nuclear.RepeatMasker.gff3` | all TE annotations (intact and fragmented copies) |
| `<unit>.nuclear.RepeatMasker.out` | the same annotation in RepeatMasker's format |
| `<unit>.nuclear.RepeatMasker.tbl` | RepeatMasker's summary of % of genome per TE class |
| `<unit>.nuclear.EDTA.intact.gff3` | structurally intact TEs |
| `<unit>.nuclear.EDTA.TElib.fa` | the genome's non-redundant, classified TE library |
| `<unit>.LTR.pass.list` | all intact LTR retrotransposons from LTR_retriever (input for LAI) |
| `<unit>.LTR.nmtf.pass.list` | the subset of intact LTRs without the canonical TGCA motif |
| `<unit>.nuclear.softmasked.fa` | nuclear genome with TEs in lower case |
