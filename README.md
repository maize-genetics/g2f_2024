# G2F maize genome resource pipeline

This is a Snakemake workflow that takes PacBio HiFi reads (and ONT reads,
where present) through to validated hifiasm assemblies and publication-ready
QC tables, for the nine Genomes to Fields maize lines.

It runs without any AI agent. Scientific decisions are recorded in
`docs/scientific_decisions.md`, the historical PHBW8 audit is in
`docs/audit_2026-10-10.md`, and the testing status of each component is in
`docs/validation_status.md`.

## Layout

```text
g2f_assemblies_2024/
├── pipeline/            # this Git repository
│   ├── config/          # config.yaml (wiring), science.yaml, methods/, samples.yaml,
│   │                    # resources/, tools/, local/ (untracked machine settings)
│   ├── workflow/        # Snakefile + rules/*.smk
│   ├── scripts/         # g2f_tool.py CLI + g2f/ helper package (stdlib + PyYAML)
│   ├── profiles/biohpc/ # Snakemake profile for running directly on a BioHPC server
│   ├── tests/           # unit, fixture, dry-run and stub end-to-end tests
│   └── docs/
├── data/<LINE>/{hifi,ont}/   # put reads here (never modified)
├── results/                  # all workflow outputs
├── tmp/                      # TMPDIR for jobs
└── hifiasm_v0.25.0/          # legacy PHBW8 runs: never written to
```

## One-time setup (BioHPC)

```bash
cd /workdir/$USER/g2f_assemblies_2024/pipeline
module load python/3.12.7            # provides Snakemake 9.26.1 + PyYAML
cp config/local/local.example.yaml config/local/local.yaml   # then edit
# pbtk is not a BioHPC module; pull its container once (~100 MB):
mkdir -p /workdir/$USER/containers
apptainer pull /workdir/$USER/containers/pbtk.sif \
    docker://quay.io/biocontainers/pbtk:3.5.0--h9ee0642_0
```

You only need pbtk if a sample has HiFi BAM inputs.

`config/tools/biohpc.yaml` lists the remaining tools: hifiasm 0.25.0,
seqkit 2.13.0, samtools 1.20, the Merqury 1.3 container (which also provides
meryl), and the BUSCO 6.1.0 container. For another site, copy that file and
point `tools:` in `config/config.yaml` at the copy.

## Inputs

Place reads in `data/<LINE>/hifi/` and, if available, `data/<LINE>/ont/`.
The directory name is the canonical sample ID.

**HiFi formats**
* `*.bam`: unaligned PacBio `hifi_reads` BAM with `READTYPE=CCS`.
* `*.fastq[.gz]`, `*.fq[.gz]`, `*.fasta[.gz]`, `*.fa[.gz]`: PacBio CCS read
  names are expected.

**ONT formats**
* Basecalled `*.fastq[.gz]` or `*.bam`.
* POD5/FAST5 files are rejected with a "basecalling required" message.

**Multiple files**
* Several files per technology are fine: multiple SMRT cells, flow cells, or
  split basecaller output.
* `.pbi`, `.bai`, `.xml`, `.md5`, `.txt` and similar files are ignored.
* Any other unrecognised file stops the workflow.

**Files that live elsewhere** (e.g. the legacy PHBW8 reads) are listed as
read-only absolute paths under `external_inputs` in
`config/local/local.yaml`. Nothing is copied or moved.

**Sample metadata and human decisions** (aliases, exclusions,
acknowledgements, k-mer reviews, method overrides) go in the tracked
`config/samples.yaml`.

### What stops the workflow

When any of the following is found, the workflow stops and prints the exact
`samples.yaml` entry that would resolve it:

* the same data appearing twice (identical copies, a FASTQ/FASTA sharing
  PacBio movies with a BAM, ONT files sharing read IDs, the same file under
  two samples);
* non-HiFi PacBio data (subreads, `fail_reads`, by-strand CCS);
* HiFi read names that cannot be verified as CCS names;
* ONT files that look already length-filtered;
* malformed, empty or unknown files;
* inputs that changed after results were built from them.

## Usage

All commands are run from `pipeline/` after `module load python/3.12.7`.

```bash
# 1. Review what will be used (reads only the first records of each file; seconds)
python3 scripts/g2f_tool.py discover
python3 scripts/g2f_tool.py discover --only PHBW8

# 2. Run the test suite (synthetic data; about a minute)
bash tests/run_tests.sh

# 3. Dry run: show every job and the exact shell commands, without running anything
snakemake --profile profiles/biohpc --dry-run --printshellcmds

# 4a. Read preprocessing, read QC and k-mer screen only (no assembly)
snakemake --profile profiles/biohpc --cores 32 reads_qc

# 4b. Full workflow (start it inside screen; hifiasm takes hours)
screen -S g2f
snakemake --profile profiles/biohpc --cores 64
```

* **Subset of samples:** add `--config only_samples=PHBW8` (comma-separated
  for several).
* **Concurrency limits:** `--cores` caps total threads. Set `resources:
  mem_mb` in the profile to match your memory reservation.

### Expected resources (from PHBW8, about 37× HiFi)

| Step | Threads | Memory | Time | Notes |
|---|---|---|---|---|
| hifiasm, HiFi only | 64 (configurable) | ~66 GB | ~3 h at 100–120 threads | |
| hifiasm, HiFi + ONT | 64 (configurable) | ~134 GB | ~4 h at 100 threads | |
| meryl k = 21 | 32 | up to 200 GB (configured) | not yet measured | |
| BUSCO | 32 | — | — | |
| Merqury | 16 | — | — | |

Disk use is roughly 35–40 GB of hifiasm intermediates per assembly, plus the
meryl database.

Thread and memory settings are in `config/resources/default.yaml`. They
never change results.

## Outputs (`results/<LINE>/`)

| Path | Contents |
|---|---|
| `inputs/inventory.json`, `inputs/md5/` | Input identities, selected method and why, full-file MD5s |
| `reads/hifi/stats/*.json`, `reads/hifi/hifi_reads.validated.json` | HiFi yield, length and quality distributions, movie census, coverage |
| `reads/ont/stats/`, `reads/ont/ul50000/` | Raw and retained ONT QC, filtered reads, merged `--ul` file, validation |
| `kmer/` | meryl database, histogram, screen, gate |
| `assembly/<method_id>/` | hifiasm outputs and logs, `<LINE>.<method_id>.{p_ctg,a_ctg}.fa` (write-protected), GFA/FASTA checks |
| `qc/<method_id>/` | seqkit stats, per-contig table, small-contig list, BUSCO, Merqury, validation records |
| `provenance/<step>/*.json` | Per-job execution record: command, tool versions, inputs/outputs, host, timing, Git commit, status (kept on failure) |
| `summary/<LINE>.<method_id>.summary.tsv` | Per-sample summary |
| `summary/<LINE>.<method_id>.validated.json` | Completion marker, written only after every validation passes |

The project table, built from validated records only, is
`results/summary/assembly_summary.tsv`.

### Work states

Snakemake tracks job states. Our part:

* **Failed jobs** keep their provenance record, with `status: failed`.
* **Incomplete jobs** are flagged by Snakemake on the next run.
* **Completed work** has its outputs present.
* **Validated work** has the `*.validated.json` marker.

`snakemake --summary` lists the state of every output.

## Resumability and safety

* **Resuming.** Re-running the same command continues where it stopped. The
  profile uses `--rerun-triggers mtime`, so editing code or config never
  silently re-runs finished jobs.
* **Changing a method.** Methods are immutable. A changed scientific setting
  means a new `method_id` file, which writes to a new results directory.
* **Write protection.** Assembly GFA and FASTA files are made write-protected
  after success, and Snakemake refuses to overwrite them.
* **Leftover hifiasm files.** hifiasm silently reuses `*.bin` files left by a
  failed attempt, so the rule refuses to start if any are present. Inspect
  them and move them aside first.
* **Changed inputs.** If the inputs for a sample change after its inventory
  was written, the workflow refuses to run. Move `results/<LINE>/` aside to
  reprocess that sample from scratch.

## Adding a maize line

1. `mkdir -p ../data/<LINE>/{hifi,ont}` and place the reads there. Omit
   `ont/` if there are no ONT data.
2. Optionally add a `samples.yaml` entry (aliases, metadata).
3. Run `python3 scripts/g2f_tool.py discover --only <LINE>` and resolve any
   reported issue.
4. Dry-run, then run, as above. The same versioned methods apply
   automatically.

## Common failures

| Message | Meaning and fix |
|---|---|
| `ONT_POSSIBLY_PREFILTERED` | Every sampled read is ≥ 50 kb. Supply the raw basecalled reads, or acknowledge the file in `samples.yaml` |
| `HIFI_POSSIBLE_DERIVED_COPY` / `DUPLICATE_COPY` | Two inputs contain the same reads. Exclude one, with a reason |
| `RAW_SIGNAL` | POD5/FAST5 found. Basecall them separately |
| `Assembly blocked pending human review` | The k-mer screen flagged possible heterozygosity. Inspect `kmer/*.hist` and record a `kmer_review` |
| `Refusing to run: hifiasm: expected version` | The tools config points at a different hifiasm build |
| `Refusing to start: intermediate hifiasm .bin files` | A previous attempt left files behind. Move them aside |
| `Write-protected output files` | A finished assembly would be overwritten. Use a new `method_id`, or deliberately move the old result |

## Downstream annotation (planned, not implemented)

Repeat annotation (EDTA), gene annotation (geneCAD) and CRE annotation
(PlantCAD) will be added as further rule files (`workflow/rules/annotation_*.smk`)
in this same workflow. Each will:

* take `results/<LINE>/assembly/<method_id>/<LINE>.<method_id>.p_ctg.fa` as
  input, gated on the sample's `validated.json` marker;
* write to `results/<LINE>/annotation/<method_id>/<tool>_<annotation_method_id>/`;
* have its own versioned method files and provenance records.

Sample IDs, inventories and provenance therefore carry through unchanged.
