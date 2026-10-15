# Workflow reference

Detailed behaviour of the assembly workflow. Start with the main
`README.md`; this page is for when you need the specifics. All commands run
from `pipeline/` with Snakemake on `PATH`.

Contents: [Inputs](#inputs) · [Expected resources](#expected-resources) ·
[Outputs](#outputs-resultsline) · [Resumability, reruns and safety](#resumability-reruns-and-safety) ·
[Running all samples](#running-all-samples) · [Adding a sample](#adding-a-sample) ·
[Common failures](#common-failures). Rerunning single QC components is
covered in `docs/rerunning_qc.md`.

## Inputs

Place reads in `data/<LINE>/hifi/` and, if available, `data/<LINE>/ont/`.
The directory name is the canonical sample ID.

**HiFi formats**
* `*.bam`: unaligned PacBio `hifi_reads` BAM with `READTYPE=CCS`.
* `*.fastq[.gz]`, `*.fq[.gz]`, `*.fasta[.gz]`, `*.fa[.gz]`: PacBio CCS read
  names are expected.

**ONT formats**
* Basecalled `*.fastq[.gz]` or `*.bam`.
* `*.tar`: an uncompressed archive of MinKNOW FASTQ, as delivered by
  sequencing facilities (`./fastq_pass/<flowcell>_pass_<run>_<n>.fastq.gz`).
  The archive is read in place and never extracted (`scripts/g2f/archive.py`):
  * FASTQ members are streamed in natural order (batch 0, 1, ..., 10, ...);
  * members under a `*fail*` directory (`fastq_fail/`) are excluded;
  * any non-FASTQ member (summaries, POD5, BAM, ...), or a mix of gzipped and
    plain members, stops the workflow (`BAD_ARCHIVE_CONTENT`);
  * truncation is detected while streaming.
  Discovery reads the archive's member headers and the first reads of the
  first member (run ID, flow cell, `sample_id`, basecall model). Tar archives
  are accepted only in `ont/`.
* POD5/FAST5 files are rejected with a "basecalling required" message.

**Multiple files**
* Several files per technology are fine: multiple SMRT cells, flow cells, or
  split basecaller output.
* `.pbi`, `.bai`, `.xml`, `.md5`, `.txt` and similar files are ignored.
* Any other unrecognised file stops the workflow.

**Files that live elsewhere** (e.g. on another filesystem) are listed as
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

### Empty genotype directories

A `data/<LINE>/` directory with no files (for example, reads not yet
transferred) is skipped with a warning, so it is never mistaken for a
processed line:

```text
[WARNING] NO_INPUTS (SAMPLE_C): sample directory contains no files; sample skipped
```

Naming such a line explicitly (`--config only_samples=SAMPLE_C`, or
`discover --only SAMPLE_C`) is an error (`REQUESTED_SAMPLE_NO_INPUTS`) and
nothing runs. The one exception is a line that has only a registered
historical assembly (`config/external_assemblies.yaml`): it gets a warning,
so `imported_qc` still works for it.

## Expected resources

Measured or estimated for a ~2.3 Gb inbred maize genome (about 37× HiFi).

| Step | Threads | Memory | Time | Notes |
|---|---|---|---|---|
| hifiasm, HiFi only | 64 (configurable) | ~66 GB | ~3 h at 100–120 threads | |
| hifiasm, HiFi + ONT | 64 (configurable) | ~134 GB | ~4 h at 100 threads | |
| meryl k = 21 | 32 | up to 200 GB (configured) | not yet measured | |
| BUSCO (2.3 Gb assembly) | 32 | 64 GB requested; peak not yet recorded | ~12 min (`poales_odb12.2`, HiFi+ONT `-l0` assembly) | |
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

## Running all samples

`scripts/run_all_assemblies.sh` (README, "All samples at full scale") runs
one Snakemake process per sample. Its defaults come from the `launcher:`
section of the top-level configuration; options override them
(`bash scripts/run_all_assemblies.sh --help`).

**Start** after a dry run, in a detached session that runs only the
launcher and closes by itself when the launcher ends, on success or failure:

```bash
bash scripts/run_all_assemblies.sh --dry-run
S=$(date +%Y%m%d_%H%M%S)
screen -dmS g2f_run bash -c "set -o pipefail; bash scripts/run_all_assemblies.sh 2>&1 | tee ../logs/launcher.$S.out; echo \$? > ../logs/launcher.$S.exit"
screen -r g2f_run      # optional: watch it (detach with Ctrl-a d)
```

The tmux equivalent is `tmux new-session -d -s g2f_run "<same command>"`.
Do not add `; exec bash` or start an interactive session and type the
command into it: that leaves an idle shell running after the job has
finished. Everything needed afterwards is in files: `launcher.$S.out`, the
exit status in `launcher.$S.exit` (absent = still running or killed), and
the per-run directory `../logs/assembly_<stamp>/`.

`--samples "A B C"` fixes the order; otherwise samples already assembled
come first, then samples with hifiasm checkpoints, then the smallest HiFi
data sets. Put the sample whose hifiasm run takes longest last: hifiasm
writes no checkpoint until all error-correction rounds have finished.

**Monitor:**

```bash
bash scripts/run_status.sh            # every sample: done/running/failed/interrupted/pending;
                                      # current step, hifiasm stage, disk, load, memory, processes
L=$(ls -d ../logs/assembly_2* | grep -v dryrun | tail -1)
tail -f $L/launcher.log               # START / END / waiting / warnings
tail -f $L/<SAMPLE>.log               # the running sample's Snakemake log ("N of M steps")
cat $L/disk.tsv                       # free disk every 10 min
```

States: `done` (exit 0), `running`, `failed` (exit code in `status.tsv`),
`interrupted` (launcher gone: crash, kill), `not_started` (the queue stopped
before it), `pending`.

**Disk.** Before each sample, the disk gate checks that its projected peak
fits above the reserve (`--reserve-gb`). If it does not, the sample waits up
to `--max-wait-h` and is then recorded `disk_wait_timeout` ("NOT STARTED:
insufficient disk space"), and the queue stops. `disk.tsv` records free
space every 10 min. A sample that is already running is not stopped.

**Stop.** Gracefully (the current sample finishes; no new one starts):
`touch "$L/state/STOP"`. Immediately: Ctrl-C in the session; this counts as
an interruption and resumes as below.

**Resume after a crash or interruption.**

1. `bash scripts/run_status.sh` shows which sample was interrupted.
2. If that sample was in hifiasm, decide about its `.bin` files. The
   workflow refuses to restart hifiasm next to unreviewed `.bin` files.
   ```bash
   S=<SAMPLE>; D=../results/$S/assembly/<method_id>
   ls -l --time-style=full-iso $D/*.bin
   grep -n 'bin files have been written\|Index has been written\|write_emask_t' $D/$S.asm.stderr.log
   ```
   * No `.bin` files: nothing to decide; hifiasm restarts from the beginning.
   * `ec.bin` and `ovlp.*.bin` exist, "bin files have been written." is in
     the log, and they were written well before the interruption: they can
     be reused. The same holds for `uidx.bin`/`ul.ovlp.bin` when "Index has
     been written" follows them. Move aside any file written close to the
     interruption, or whose completion is not logged. Then write the
     approval, which hifiasm consumes and records in
     `$S.asm.reused_bin_files.txt`:
     `echo "<name>, <date>: reuse ec/ovlp (logged complete <time>), <reason>" > $D/$S.asm.reuse_bins.approved`
   * Not sure: move all `.bin` files aside; hifiasm starts over.
3. Relaunch the same command.

The launcher then refuses if another launcher, or a Snakemake/hifiasm
process in this project, is alive; marks samples left `running` as
`interrupted` in their old log directory; removes the stale Snakemake lock;
re-runs interrupted jobs (`--rerun-incomplete`); skips completed work; and,
before restarting hifiasm, keeps the log of the unfinished attempt as
`<SAMPLE>.asm.stderr.attempt-<time>.log` (Snakemake deletes a job's
declared log on restart).

What an interruption cannot corrupt: assembly outputs are write-protected
once complete; a hifiasm result counts only after hifiasm exits 0, the
primary GFA has sequences and the alternate GFA exists; provenance records
keep status `running` until a step finishes; and Snakemake's
incomplete-job markers make it re-run anything that was cut off.

hifiasm checkpoints: hifiasm writes reusable `.bin` files only at fixed
points: after all error-correction rounds and the final overlap pass
(`ec.bin`, `ovlp.*.bin`), after the ultra-long-read index and overlaps
(`uidx.bin`, `ul.ovlp.bin`), and during UL re-integration (`re.*`). Reusing
them always goes through the written approval above. The thread count of a
resumed run may differ; it is recorded in the provenance.

## Adding a sample

1. `mkdir -p ../data/<SAMPLE>/{hifi,ont}` and place the reads there. Omit
   `ont/` if there are no ONT data.
2. Optionally add a `config/samples.yaml` entry (aliases, metadata).
3. Run `python3 scripts/g2f_tool.py discover --only <SAMPLE>` and resolve any
   reported issue.
4. Dry-run, then run, as above. The same versioned methods apply
   automatically.

## Common failures

| Message | Meaning and fix |
|---|---|
| `NO_INPUTS` (warning) | That `data/<LINE>/` directory is empty; the line is skipped. Add its reads when available |
| `REQUESTED_SAMPLE_NO_INPUTS` | A line named in `only_samples` has no read files. Add them, or remove the line from `only_samples` |
| `ONT_POSSIBLY_PREFILTERED` | Every sampled read is ≥ 50 kb. Supply the raw basecalled reads, or acknowledge the file in `samples.yaml` |
| `HIFI_POSSIBLE_DERIVED_COPY` / `DUPLICATE_COPY` | Two inputs contain the same reads. Exclude one, with a reason |
| `RAW_SIGNAL` | POD5/FAST5 found. Basecall them separately |
| `WARNING: k-mer screen flagged …` | A diagnostic only; assembly continues. Inspect `kmer/*.hist` and, if appropriate, record a `kmer_review` in `samples.yaml` |
| `WARNING: BUSCO lineage poales_odb12.2 is not installed` | Download it (see "Installation" in the README) |
| `lineage … pinned dataset requires …` | The installed BUSCO dataset is not the pinned release |
| `was changed after it was used` | A method file was edited in place. Restore it, or create a new `method_id` |
| `read statistics were produced by different read-statistics code` | Regenerate the statistics with `--forcerun hifi_read_stats ont_input_stats ont_filtered_stats` |
| `Refusing to run: hifiasm: expected version` | The tools config points at a different hifiasm build |
| `Refusing to start: intermediate hifiasm .bin files` | A previous attempt left files behind. Move them aside, or, to resume from hifiasm's checkpoints, inspect them and write who/why into `<sample>.asm.reuse_bins.approved` in the assembly directory (consumed on start; the reused files are listed in `<sample>.asm.reused_bin_files.txt`) |
| `another launcher is running for this project` | A launcher holds `logs/launcher.lock`. Use `bash scripts/run_status.sh`; wait for it, or STOP it |
| `Snakemake/hifiasm processes are still running in this project` | A run started outside the launcher (or an older launcher) works in this pipeline or results directory. Let it finish or stop it, then launch |
| `NOT STARTED: insufficient disk space` | The genotype's projected peak does not fit above the reserve. Free space, then relaunch (completed work is kept) |
| `Write-protected output files` | A finished assembly would be overwritten. Use a new `method_id`, or deliberately move the old result |
