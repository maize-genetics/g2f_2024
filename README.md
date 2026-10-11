# G2F maize genome resource pipeline

This is a Snakemake workflow that takes PacBio HiFi reads (and ONT reads,
where present) through to validated hifiasm assemblies and publication-ready
QC tables, for the nine Genomes to Fields maize lines.

It runs without any AI agent. Supporting documentation:

* Scientific decisions: `docs/scientific_decisions.md`.
* Historical PHBW8 audit: `docs/audit_2026-10-10.md`.
* Real-tool compatibility findings: `docs/realtools_validation_2026-10-11.md`.
* Testing status of each component: `docs/validation_status.md`.

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
├── containers/               # pinned containers not provided by BioHPC (pbtk)
├── tmp/                      # TMPDIR for jobs
└── hifiasm_v0.25.0/          # legacy PHBW8 runs: read-only, registered for QC
```

## One-time setup (BioHPC)

```bash
cd /workdir/$USER/g2f_assemblies_2024/pipeline
module load python/3.12.7            # provides Snakemake 9.26.1 + PyYAML
cp config/local/local.example.yaml config/local/local.yaml   # then edit
```

**pbtk 3.5.0** is needed only for HiFi BAM inputs. It is not a BioHPC
module; pull it by digest (19 MB, already done on cbsuxm01):

```bash
mkdir -p ../containers ../tmp/apptainer_cache
APPTAINER_CACHEDIR=$PWD/../tmp/apptainer_cache apptainer pull ../containers/pbtk-3.5.0.sif \
  docker://quay.io/biocontainers/pbtk@sha256:ab2d8dcf1a80e2d7595ab8799c37037766db6fb35d6744440b4643f27749e0ab
```

**BUSCO lineage `poales_odb12.2`** must be installed under
`busco_lineage_dir` (about 1 GB; needs network access):

```bash
apptainer exec --bind /home/$USER/busco_downloads /programs/busco-6.1.0/busco.sif \
  busco --download_path /home/$USER/busco_downloads --download poales_odb12.2
cat /home/$USER/busco_downloads/lineages/poales_odb12.2/dataset.cfg   # check against science.yaml pins
```

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

# 2. Run the test suite (synthetic data, stub tools; about 1.5 minutes)
bash tests/run_tests.sh
#    plus real-tool compatibility (real pbtk, hifiasm, meryl, Merqury, BUSCO on 300 kb data; ~75 s)
G2F_REALTOOLS=1 bash tests/run_tests.sh test_realtools

# 3. Dry run: show every job and the exact shell commands, without running anything
snakemake --profile profiles/biohpc --dry-run --printshellcmds

# 4a. Read preprocessing, read QC and k-mer screen only (no assembly)
snakemake --profile profiles/biohpc --cores 32 reads_qc

# 4b. QC of the registered historical assemblies only (no reads, no assembly)
snakemake --profile profiles/biohpc --cores 32 imported_qc

# 4c. Full workflow (start it inside screen; hifiasm takes hours)
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
| BUSCO (2.3 Gb assembly) | 32 | 64 GB requested | minutes | Legacy `liliopsida_odb10` runs took ~3 min at 40 CPUs; `poales_odb12.2` is larger, so somewhat longer is expected |
| GFA/FASTA check, seqkit stats (2.2 GB) | 1–4 | ≤ 8 GB | minutes | Reads both files once |
| HiFi read stats + MD5 (28 GB) | 4 + 1 | ≤ 8 GB | not yet measured (expected tens of minutes) | Streaming passes |
| Merqury | 16 | — | not yet measured | Needs the meryl database |

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
| `summary/<LINE>.<method_id>.summary.<VID>.<SID>.tsv` | Per-sample summary |
| `summary/<LINE>.<method_id>.validated.<VID>.<SID>.json` | Completion marker, written only after every validation passes |
| `imported/<id>/`, `qc/imported-<id>/` | Links to registered historical assemblies (read-only) and their QC |

The project tables are built from validated records only:

* `results/summary/assembly_summary.<VID>.<PID>.tsv`: assemblies produced by
  this workflow.
* `results/summary/imported_assemblies_summary.<VID>.<PID>.tsv`: historical
  assemblies, kept separate.

**Name components**

* `VID` is the validation identity: `science.yaml`, the method definitions
  and the validation code. `results/.identity/validation-<VID>.json` lists
  what it is made of.
* `SID` / `PID` hash the human decisions in `samples.yaml`, for one sample
  or for the whole set.
* The newest file is the current one. A file with an older identity was
  validated under older settings and is kept for the record.

### Work states

Snakemake tracks job states. Our part:

* **Failed jobs** keep their provenance record, with `status: failed`.
* **Incomplete jobs** are flagged by Snakemake on the next run.
* **Completed work** has its outputs present.
* **Validated work** has the `*.validated.json` marker.

`snakemake --summary` lists the state of every output.

## Resumability, reruns and safety

* **Resuming.** Re-running the same command continues where it stopped.
  Checksums and full-file read statistics are computed once and reused.
* **Resource changes** (threads, memory, tmpdir, tool paths) never re-run
  anything.
* **Science, method and validation-code changes.** The profile uses
  `--rerun-triggers mtime`, so expensive jobs never re-run silently.
  Instead:
  * A changed `science.yaml` or validation code gives a new validation
    identity. Only cheap validation and summary jobs re-run, and they write
    summaries under new names. Assemblies, k-mer counting, BUSCO and Merqury
    are not repeated.
  * Assembly parameters live in immutable method files. Editing a method
    file that has already been used is refused; create a new `method_id`
    instead, which writes to a new results directory.
  * Tool versions required by a method are checked at run time.
* **Write protection.** Assembly GFA and FASTA files are write-protected;
  Snakemake refuses to overwrite them.
* **Failed jobs.** Snakemake removes their outputs, and the provenance
  record keeps `status: failed`. hifiasm refuses to start while `.bin` files
  from a failed attempt exist; inspect them and move them aside.
* **Changed inputs.** If a sample's inputs or acknowledgements change after
  results exist, the workflow refuses to run. Move `results/<LINE>/` aside
  to reprocess that sample.
* **Historical files.** A registered historical file that changes after
  import makes the workflow refuse to run.

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
| `WARNING: k-mer screen flagged …` | A diagnostic only; assembly continues. Inspect `kmer/*.hist` and, if appropriate, record a `kmer_review` in `samples.yaml` |
| `WARNING: BUSCO lineage poales_odb12.2 is not installed` | Download it (see "One-time setup") |
| `lineage … pinned dataset requires …` | The installed BUSCO dataset is not the pinned release |
| `was changed after it was used` | A method file was edited in place. Restore it, or create a new `method_id` |
| `read statistics were produced by different read-statistics code` | Regenerate the statistics with `--forcerun hifi_read_stats ont_input_stats ont_filtered_stats` |
| `Refusing to run: hifiasm: expected version` | The tools config points at a different hifiasm build |
| `Refusing to start: intermediate hifiasm .bin files` | A previous attempt left files behind. Move them aside |
| `Write-protected output files` | A finished assembly would be overwritten. Use a new `method_id`, or deliberately move the old result |

## Downstream annotation (planned, not implemented)

EDTA (repeats), geneCAD (genes), PlantCAD/PlantCAD2 (CREs) and comparative
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
