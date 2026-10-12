# Full-scale run plan: all nine genotypes (prepared 2026-10-12; NOT started)

Launcher: `scripts/run_all_assemblies.sh`. Resources:
`config/resources/cbsuxm01_3x40.yaml`. Nothing here has been executed except
`--dry-run`.

## What runs

For every genotype, the existing workflow (target `per_sample`):
input inventory → checksums → HiFi BAM → FASTQ (pbtk) → HiFi/ONT read QC →
ONT ≥ 50 kb filter (streamed from the tar archives) → meryl + k-mer screen →
hifiasm 0.25.0 `-l0` with `--ul` ONT (method selected automatically; all nine
genotypes have ONT, so all use `hifiasm-0.25.0_hifi-ont-ul50k_l0_v1`) →
GFA → FASTA with consistency checks → assembly statistics, BUSCO 6.1.0
(`poales_odb12.2`), Merqury → validated per-sample summaries. When every
genotype has succeeded, the launcher builds the project summary table.

Dry run (2026-10-12): 340 jobs in total. Per genotype: 3AZA1 40, GEMS-0227 44,
FBLL 44, GEMS-0219 48, PHP02 32, GEMS-0061 40, GEMS-0113 36, PHB47 44,
PHBW8 28.

## Execution model

* One Snakemake process per genotype (`--config only_samples=<LINE>`), at
  most **3 at a time** (GNU parallel), **40 cores each** (120 total),
  **300 GB memory budget each** (Snakemake's `mem_mb` scheduling). The next
  genotype starts as soon as one finishes.
* Start order: largest HiFi input first (3AZA1, GEMS-0227, FBLL, ...), so
  the longest assembly does not run last.
* A genotype that fails does not stop the others (`--halt never`;
  `--keep-going` inside each genotype). Its log is
  `../logs/assembly_<stamp>/<LINE>.log`, the per-genotype exit codes are in
  `status.tsv`, and `joblog.tsv` has GNU parallel's timings.
* **Resume:** rerun the same command. Snakemake skips completed outputs
  (`--rerun-triggers mtime`, protected assemblies). A failed hifiasm leaves
  its `.bin` files and refuses to restart until a person moves them aside
  (existing safety rule).
* Concurrent processes in one working directory were tested
  (`test_concurrent_per_sample_processes`). The shared validation-identity
  stamp is created once before the genotypes start.

## Resource audit

### Measured now (2026-10-12)

| resource | value |
|---|---|
| CPUs | 256 logical (2 × AMD EPYC 7702, 64 cores × 2 threads) |
| RAM | 2.0 TB; 942 GB "available" plus ZFS ARC cache 1,075 GB (reclaimable down to 67 GB) |
| Disk `/local/workdir` (ZFS, no quota) | 22 TB, 95 % used, **1.13 TB free** (`zfs available` 1.10 TiB). Not the 1.6 TB previously noted |
| Other users | load ~0.3 (idle) before our pilot; the disk is shared with other users of this server |

### Estimates per genotype

The basis is PHBW8: 84.3 Gb of HiFi from 50.8 GB of BAM (0.60 bases/byte,
confirmed on a small FBLL BAM read in full: 7.78 Gb / 4.76 GB). Its
HiFi+ONT `-l0` hifiasm took 1,007,460 CPU-s and peaked at 133.7 GB RSS. Its
`.bin` intermediates were 52 GB (ec 21.4, ovlp 17.6, ONT uidx 13.2). Its
GFA+FASTA outputs are ~10 GB. meryl/Merqury sizes are rough estimates
(not yet measured at maize scale).

| genotype | HiFi BAM GB | est. HiFi Gb (cov.) | converted FASTQ GB (temp) | hifiasm .bin GB (removed after success) | meryl+Merqury GB | kept GB | hifiasm CPU-h | est. hifiasm wall @ 40 thr | est. hifiasm peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 3AZA1 | 84.6 | 140 (61x) | 93 | 75 | 54 | 69 | 466 | ~14 h | ~223 GB |
| GEMS-0227 | 58.6 | 97 (42x) | 64 | 56 | 41 | 56 | 323 | ~9.5 h | ~154 GB |
| FBLL | 58.0 | 96 (42x) | 64 | 55 | 41 | 56 | 320 | ~9.4 h | ~153 GB |
| GEMS-0219 | 48.1 | 80 (35x) | 53 | 48 | 36 | 51 | 265 | ~7.8 h | ~127 GB |
| PHP02 | 46.0 | 76 (33x) | 51 | 47 | 35 | 50 | 253 | ~7.5 h | ~121 GB |
| GEMS-0061 | 45.7 | 76 (33x) | 50 | 46 | 35 | 50 | 252 | ~7.4 h | ~120 GB |
| GEMS-0113 | 43.5 | 72 (31x) | 48 | 45 | 34 | 49 | 240 | ~7.0 h | ~114 GB |
| PHB47 | 41.2 | 68 (30x) | 45 | 43 | 33 | 48 | 227 | ~6.7 h | ~108 GB |
| PHBW8 | (FASTA in place) | 84 (37x) | 0 | 50 | 37 | 52 | 280 | ~8.2 h | ~134 GB |
| **total** | | ~790 Gb | **468** | **465** | **345** | **480** | **2,625** | | |

Wall time assumes 85 % parallel efficiency. Preprocessing and QC add
~3-5 h per genotype (checksums, pbindex + bam2fastq, read statistics,
meryl, BUSCO, Merqury). Expected total: **~35-45 h** for all nine with
three slots.

### Disk: feasible only with cleanup

* **Without cleanup:** 468 + 465 + 480 ≈ **1.41 TB**. This exceeds the
  1.13 TB free. **Not feasible.**
* **With the two cleanup switches** in `cbsuxm01_3x40.yaml` (below):
  persistent ≈ 480 GB at the end. The transient load is at most three
  genotypes' converted FASTQ + `.bin` files: ≤ 410 GB for the three largest,
  ~230 GB for the last three. Peak ≈ **710 GB**, worst case (estimate errors)
  ≈ **890 GB**, against 1.13 TB free. That leaves a ~240-400 GB margin,
  which other users of the disk can consume. The launcher refuses to start
  below 800 GB free (`--min-free-gb`).
* Further reduction if needed: remove each genotype's read meryl DB after
  Merqury (~25-40 GB each, regenerable). This is **not** implemented and
  would need approval.

### Memory

* hifiasm is the largest consumer: three concurrent largest ≈
  223 + 154 + 153 ≈ **530 GB**. meryl is configured for 120 GB (it adapts
  its batches).
* Per-genotype budget of 300 GB: hifiasm (240 GB declared) and meryl
  (120 GB) of one genotype never run at the same time. Three genotypes
  ≤ 900 GB declared, below the ~940 GB available without touching the ZFS
  cache.
* Snakemake memory values are scheduling hints, not hard limits. If 3AZA1
  needs more than estimated, it simply uses it; there is ample headroom.

### CPU and overlap

Each genotype is capped at 40 cores, so heavy QC (BUSCO 32 threads,
Merqury 16, meryl 40) of one genotype can overlap another genotype's
hifiasm without exceeding 120 in total. Within a genotype, jobs are
scheduled inside its 40 cores. No oversubscription is possible beyond
120 threads (one extra decompression helper thread per streaming job aside).

## Cleanup rules (require approval)

`config/resources/cbsuxm01_3x40.yaml`:

| switch | what is removed | when | safety |
|---|---|---|---|
| `temp_converted_hifi: true` | `results/<LINE>/reads/hifi/converted/*.fastq.gz` | by Snakemake, once read stats, meryl and hifiasm have used it | regenerable bit-for-bit from the original BAM; BAM never touched |
| `remove_hifiasm_bins_after_success: true` | `results/<LINE>/assembly/<method>/<LINE>.asm.*.bin` | immediately after hifiasm exits 0 | listed first in `<LINE>.asm.removed_bin_files.txt`; failed runs keep their bins |

Nothing pre-existing is deleted: raw data, PHBW8's preprocessed inputs, the
legacy PHBW8 assemblies (`hifiasm_v0.25.0/`), and earlier results are never
touched. Both switches are tested (`test_cleanup_switches`).

## Preflight (automatic, in the launcher)

1. Snakemake ≥ 8, GNU parallel, every configured tool and container exists,
   BUSCO lineage `poales_odb12.2` present.
2. Free disk ≥ 800 GB; memory and load reported.
3. Whole-project dry run (discovery, configuration, DAG) must succeed.
4. `launch.json` records git commit (and uncommitted changes), settings,
   genotypes and host.

## Blockers and decisions before launch

1. **GEMS-0219 ONT run `10474838_9`** (52 reads, 180 kb, longest 21.5 kb) has
   no reads ≥ 50 kb. The workflow then stops GEMS-0219 at the ONT gate by
   design ("does not silently drop ONT inputs"). Decide one of:
   * exclude it (recommended; it contributes nothing to `--ul`):
     ```yaml
     # config/samples.yaml
     GEMS-0219:
       exclude_inputs:
         10474838_9_PAQ86658_8f1582f4_fastq.tar: >
           Aborted 52-read ONT run (180 kb, longest read 21.5 kb); no reads >= 50 kb.
           Excluded <date> by <name>.
     ```
   * or acknowledge it (`acknowledge: {<file>: {no_qualifying_reads: <reason>}}`)
     so it is kept in the yield accounting.
2. **Fingerprinting results** (`docs/fingerprinting_results_2026-10-12.md`):
   review any ONT run that is not `supported` before assembly. HiFi+ONT
   assembly uses every ONT input of a genotype.
3. **Approve the cleanup switches** (above), or run with
   `--resources config/resources/<copy without cleanup>.yaml` after freeing
   disk.
4. Commit the pipeline changes so `launch.json` records a clean commit.

## Start and monitor

```bash
cd /local/workdir/naj49/g2f_assemblies_2024/pipeline
bash scripts/run_all_assemblies.sh --dry-run          # review: ../logs/assembly_<stamp>.dryrun/
mkdir -p ../logs
nohup bash scripts/run_all_assemblies.sh > ../logs/launcher.$(date +%Y%m%d_%H%M%S).out 2>&1 &
# (or run it inside `screen`)
cat ../logs/assembly_*/status.tsv                      # finished genotypes and exit codes
tail -f ../logs/assembly_*/launcher.log                # START/END events
tail ../logs/assembly_*/<LINE>.log                     # one genotype's Snakemake log
df -h /local/workdir                                   # disk
```

To stop: `kill` the launcher's process group (or Ctrl-C inside screen). Then
rerun the same command to resume. Snakemake marks interrupted jobs as
incomplete, and `--rerun-incomplete` regenerates them. An interrupted hifiasm
keeps its `.bin` files and will refuse to restart until they are moved aside.
