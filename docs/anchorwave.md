# Gene-anchored structural validation with AnchorWave

`scripts/g2f_anchorwave.py` aligns contig assemblies to an annotated reference
with [AnchorWave](https://github.com/baoxingsong/AnchorWave) `proali`, and
`scripts/g2f_anchorwave_synteny.py` summarises the result. Use it to look at
collinearity, rearrangements and reference loci that an assembly represents
more than once. It does not purge or edit assemblies, and a difference from the
reference is an observation, not an error call.

## Workflow

| stage | command | output |
|---|---|---|
| `ref_cds` | `anchorwave gff2seq` (+ FASTA/GFF3 id, length and optional md5 checks) | `reference/<REF>/cds/cds.fa` |
| `ref_map` | `minimap2 -x splice -k 12 -a -p 0.4 -N 20` CDS -> reference | `ref_map/ref.sam`; fails if < 90% of CDS map back to their own locus |
| `query_map:S:A` | same CDS -> assembly | `query_map/query.sam`; checks sha256 / contig count / length when configured |
| `anchors:S:A:V` | `proali -n ... -ns` for each variant V | `anchors_V/anchors` (minutes) |
| `align:S:A:V` | `proali -n -o -f -b` for `align_variants` | `alignment.maf`, `alignment.f.maf`, `methods.bed` (hours) |

Stages run in `<dir>.partial` and become `<dir>` only after validation;
`VALIDATED.json` records inputs (size, mtime), parameters and their hash, tool
versions, the container sha256, the script sha256, preflight resource
snapshots, per-command wall time and peak RSS. Failed attempts keep a
`FAILED.json` and are moved to `<dir>.attempt-<stamp>` before a retry; a
parameter or input change makes a stage `stale`, rerun only with
`--rerun-stale` (old result moved to `<dir>.previous-<stamp>`). A final
directory without `VALIDATED.json` is never touched.

```
python3 scripts/g2f_anchorwave.py -c config/local/anchorwave.yaml check-env
python3 scripts/g2f_anchorwave.py -c config/local/anchorwave.yaml plan
nohup python3 scripts/g2f_anchorwave.py -c config/local/anchorwave.yaml run > LOGDIR/aw.log 2>&1 &
python3 scripts/g2f_anchorwave.py -c config/local/anchorwave.yaml status
```

Run the `run` command in a session that survives logout (`nohup`, `screen`
or `tmux`) and that ends when the run ends; `status` can be checked at any
time. Start with `--stages ref_cds,ref_map,query_map,anchors`, inspect the anchors,
then run `--stages align`. Configuration: `config/anchorwave.example.yaml`
(paths, resources) and `config/anchorwave_science.yaml` (parameters).

## Resources and safety

Every command runs under `taskset`/`nice`. Before each stage the runner
refuses to start if free disk, `MemAvailable` or idle CPUs on the chosen CPU
list are below the configured floors. With `memory_cap: systemd-scope` each
command runs in a transient systemd user scope with `MemoryMax = max_rss_gb`
and no swap: the kernel kills only that command if it exceeds the limit, and
the run refuses to start if the scope cannot be created. A watchdog (every
10 s) also stops the stage's process group if its RSS exceeds `max_rss_gb` or
system memory or disk fall below their floors. `anchorwave_mem_gib` is passed
to `proali -M`, a predictive envelope, not a hard limit. If the runner receives
SIGTERM or SIGHUP (e.g. its terminal session closes), it stops the running
command's process group, records the stage as interrupted (`FAILED.json`) and
exits with status 143; rerunning resumes.

SIGKILL (or the loss of the runner process for any other reason) cannot be
handled: the running command survives, still inside its memory cap but without
the watchdog. Each command's identity (named systemd scope, process group,
start time, runner pid) is written to `process.json` in the stage directory;
`status` reports a stage whose recorded command is still alive as `orphaned`
(with PIDs and the command to stop it), and `run` refuses to start, without
moving or signalling anything, until those processes are gone. Inspect the
listed PIDs and stop them yourself (`systemctl --user stop <unit>.scope`).

Each `run` writes `<outdir>/<profile>/runs/<run_id>/RUN.json` (git commit and
modified tracked files, container sha256 and tool versions, memory cap,
reference paths and checksums, assembly paths and expected checksums, config
and science file sha256, selected stages, final stage states) with copies of
the configuration and science file; every stage record carries the `run_id`.
Inputs are only read: outputs, logs and temporary files go to the stage
directories, `runs/` and `tmpdir`. Run from a clean checkout so that the
recorded commit identifies the code.

## Parameters

* `-R` is the maximum alignment depth of each **reference** position and `-Q`
  of each **query** position (AnchorWave README: maize reference, sorghum query,
  `-R 1 -Q 2`). For an unpurged or partly purged assembly of a heterozygous
  sample, both alleles of a locus can be kept, so the baseline is `-R 2 -Q 1`;
  `R1Q1` (strict one-to-one) and `R3Q1` (capacity probe) are run as
  anchor-only comparisons.
* `-e` limits copies of one gene on one query sequence; `-e 2` is run as a
  comparison for tandem/same-contig duplication.
* `-D` (maximum chain gap in anchors, default 25): an inversion containing
  fewer than about `D` genes is bridged by the surrounding collinear block and
  not reported (`tests/test_anchorwave.py`). The `R2Q1e1D5` anchor-only
  variant screens for such small rearrangements; it also splits blocks at
  presence/absence-variable gene clusters.
* Gene anchors cannot see rearrangements in gene-free sequence (e.g. repeat
  arrays or a breakpoint segment without annotated genes). Combine with
  sequence-level alignments.

## Summaries

```
python3 scripts/g2f_anchorwave_synteny.py --anchors .../anchors_R2Q1e1/anchors \
    --variant D5=.../anchors_R2Q1e1D5/anchors --query-fai ASM.fa.fai --ref-fai REF.fa.fai \
    [--maf .../alignment.maf] [--asmqc-assembly-dir ...] [--asmqc-pair-dir ...] \
    --zoom chr2:80000000-130000000 --out NEW_DIR
```

Outputs: `blocks.tsv`, `contig_assignment.tsv`, `junctions.tsv` (adjacent
blocks along a contig: inter-chromosomal, strand switch, order reversal,
length-discordant gap; HiFi depth/spanning from g2f_asmqc when given),
`duplicated_regions.tsv` (reference loci anchored by two assembly copies, with
depth-based candidate classes), `reference_windows.tsv`, MAF coverage and
identity, dot plots. `scripts/g2f_breakpoint_reads.py` counts reads spanning
given positions in BAM/CRAM files.

Copy classes are heuristics: two copies at one-allele depth with high identity
are *candidate alternative alleles*; two copies at full depth are *candidate
genuine duplications*; everything else is unresolved. Syntenic overlap alone
does not establish allelic redundancy.

Tests: `ANCHORWAVE_IMAGE=/path/anchorwave.sif python3 -m unittest tests/test_anchorwave.py`
(synthetic genome; the container tests are skipped without the image).
