# Full-scale run plan: all nine genotypes (approved 2026-10-12)

Launcher: `scripts/run_all_assemblies.sh` (per genotype:
`scripts/run_one_genotype.sh`). Resources and retention:
`config/resources/cbsuxm01_3x40.yaml`.

## Approved decisions (project lead, 2026-10-12)

* **HiFi+ONT for all nine genotypes:** hifiasm 0.25.0 `-l0`
  (`hifiasm-0.25.0_hifi-ont-ul50k_l0_v1`), with every valid HiFi and ONT run.
  There are no HiFi-only comparison assemblies.
* **ONT ultra-long filter unchanged:** reads ≥ 50,000 bp (`seqkit seq -m 50000`,
  D3; the historical name is "> 50 kb").
* **GEMS-0219 `10474838_9`** (aborted, 52 reads) is excluded in
  `config/samples.yaml`. Runs `_9b` and `_9c` are used.
* **GEMS-0061 and GEMS-0113:** ONT assignments are accepted provisionally
  (fingerprinting "ambiguous"; cause not determined). This is recorded in
  `samples.yaml` metadata and `docs/fingerprinting_results_2026-10-12.md`.
* **Three genotypes × 40 threads**, at most 120 in total.
* **hifiasm `.bin` cleanup** after verified success is approved.
* **Local raw BAM/tar copies** are deleted **manually** by the project lead
  after verified preprocessing. The pipeline only reports eligibility.

## Order of operations per genotype

1. **Raw inputs validated:** discovery (formats, identity, duplicates),
   full-file MD5 checksums, BAM primary-record counts.
2. **Preprocessing:**
   * HiFi: BAM → FASTQ (pbtk);
   * ONT: tar members streamed → ≥ 50 kb filter, with no extracted copy;
   * read QC of both.
3. **Prepared-input gate** (`verify_prepared_inputs`):
   * one full pass over every prepared file (gzip integrity, MD5, records);
   * HiFi FASTQ records = read-QC count = BAM primary records;
   * ONT filtered records = filtered-read QC count, all reads ≥ 50 kb,
     nonzero yield;
   * every hifiasm input is a verified file;
   * runs reconciled with `data/g2f_9genotypes_raw_reads_manifest.tsv`: no
     missing, no unexpected, excluded runs accounted for.

   Output: `results/<LINE>/inputs/prepared_inputs.verified.json` (per run:
   original path, manifest source path, raw size and MD5, prepared path,
   size, MD5, records, bases). hifiasm depends on this file.
4. **Raw-copy eligibility:** at the end of each genotype, the launcher runs
   `g2f_tool.py raw-release` (report only).
5. **hifiasm** reads only the prepared files: converted HiFi FASTQ, filtered
   ONT FASTQ, and PHBW8's derivatives in place.
6. **Assembly verified, then `.bin` cleanup.** Inside the hifiasm job,
   after exit 0, the p_ctg GFA must contain sequence lines and the a_ctg GFA
   must exist. Only then are the `.bin` files listed (in
   `<LINE>.asm.removed_bin_files.txt` and the genotype log) and removed.
   Failed or interrupted runs keep them.
7. **QC:** GFA→FASTA consistency, assembly statistics, BUSCO 6.1.0
   (`poales_odb12.2`), meryl k-mer screen, Merqury → validated summaries.
8. **Retention** as below. When all genotypes succeed, the project summary
   table is built.

## Retention policy

| priority | files | policy |
|---|---|---|
| 1 | assemblies (GFA, FASTA), QC, logs, provenance, configuration, inventories, verification records | always kept |
| 2 | prepared hifiasm inputs: converted HiFi FASTQ (`reads/hifi/converted/`), filtered ONT FASTQ (`reads/ont/ul50000/`) | ONT always kept (small). HiFi kept unless the disk gate decides otherwise (below) |
| 3 | regenerable intermediates: hifiasm `.bin`, the pbtk index (`*.pbtk/`), Snakemake temp files | `.bin` removed after verified success; others small |

Categories: **needed by an active job** (raw files during preprocessing;
prepared files during hifiasm and meryl; `.bin` during hifiasm), **needed for
reproducibility** (priorities 1 and 2), **disposable** (priority 3).

**Prepared HiFi FASTQ retention** (`--retain-hifi auto`, the default).
Before each genotype starts, the disk gate keeps its FASTQ only if, after
reserving that genotype's full peak, the remaining space still covers every
not-yet-started genotype's final footprint and the largest concurrent
transient needs. Otherwise the FASTQ is made temporary: Snakemake deletes it
once read QC, the gate, meryl and hifiasm have used it. The decision is
recorded per genotype (`status.tsv`, `gate.jsonl`). Deleting local raw
copies frees space, so later genotypes then keep theirs.
`--retain-hifi keep|temp` forces one behaviour.

## Raw copies: when are they safe to delete?

`python3 scripts/g2f_tool.py raw-release` (run by the launcher after each
genotype and at the end; safe to rerun) writes
`results/raw_release_report.tsv` and `results/<LINE>/inputs/raw_release.json`.
A local raw BAM or tar is **ELIGIBLE** only if all of these hold:

1. the genotype's prepared inputs passed the gate (all runs accounted for);
2. its prepared file is **still present and unchanged** since verification
   (same size and mtime). So a genotype whose HiFi FASTQ was not retained
   keeps its BAMs: they are then the only local source;
3. its raw MD5 was recorded;
4. it is in the manifest, **and** `scripts/check_remote_copies.sh` confirmed
   the original exists on the source node at the manifest size. That check
   needs your interactive SSH login: key-based SSH from this server to
   cbsublfs1 is refused, and `/data4`, `/data5` are not mounted here.

hifiasm and QC depend only on the prepared files. After an eligible raw
file is deleted, the workflow reconstructs that input from the recorded
inventory and resumes with nothing to redo (tested:
`test_raw_release_lifecycle`). Deleting any **non-eligible** file changes
the inventory, and the workflow stops with an explanation. Excluded or
in-place inputs (GEMS-0219 `_9`; PHBW8's derivatives) are never listed as
eligible; keep them.

```bash
bash scripts/check_remote_copies.sh          # once (asks for your password); writes ../data/remote_copies.tsv
python3 scripts/g2f_tool.py raw-release      # refresh; prints ELIGIBLE / keep per file with reasons
```

## Disk safeguards

* **Projection** (`scripts/g2f/diskplan.py`; factors from PHBW8): peak
  during hifiasm = converted FASTQ + `.bin` + meryl/Merqury + outputs +
  filtered ONT; final footprint = the same without `.bin` (and without the
  FASTQ if not retained).

  | genotype | peak GB | final GB (FASTQ kept) | final GB (FASTQ deleted after use) |
  |---|---|---|---|
  | 3AZA1 | 240 | 165 | 72 |
  | GEMS-0227 | 180 | 124 | 59 |
  | FBLL | 178 | 123 | 59 |
  | GEMS-0219 | 154 | 106 | 53 |
  | PHP02 | 151 | 104 | 54 |
  | GEMS-0061 | 149 | 102 | 52 |
  | GEMS-0113 | 142 | 98 | 50 |
  | PHB47 | 137 | 94 | 49 |
  | PHBW8 | 102 | 52 | 52 |
  | **total** | | **969** | **500** |

  Worst-case peak with every FASTQ deleted after use ≈ 500 + the three
  largest transients ≈ 0.93 TB. Free at launch ≈ 1.21 TB, plus ~530 GB when
  the raw copies are deleted.
* **Start gate** (`g2f_tool.py disk-gate`, under a file lock): a genotype
  starts only if free space − (remaining growth of running genotypes) −
  reserve (150 GB) ≥ its peak. Otherwise it waits (rechecked every 10 min,
  at most 48 h). A resumed genotype only needs what is not already on disk.
* **Monitor:** free space and results size every 10 min in `disk.tsv`. A
  warning goes to `launcher.log` when free space drops below the reserve.
* **Drain:** `touch <logdir>/state/STOP` starts no new genotypes; running
  ones finish.

## Resources (unchanged from the audit)

256 logical CPUs; 2 TB RAM (ZFS cache reclaimable). hifiasm ≈ 2,600 CPU-h
in total, with peak RSS ≈ 108-223 GB per genotype; 300 GB scheduling budget
per genotype. Expected wall time ≈ 35-45 h.

## Start, monitor, resume

```bash
cd /local/workdir/naj49/g2f_assemblies_2024/pipeline
bash scripts/run_all_assemblies.sh --dry-run
nohup bash scripts/run_all_assemblies.sh > ../logs/launcher.$(date +%Y%m%d_%H%M%S).out 2>&1 &

L=$(ls -dt ../logs/assembly_2*[0-9] | head -1)
cat $L/status.tsv              # finished genotypes: exit code, FASTQ retention, times
tail -f $L/launcher.log        # START / waiting / END events, disk warnings
tail $L/<LINE>.log             # one genotype's Snakemake log
cat $L/disk.tsv                # free space over time
ls $L/state                    # <LINE>.running / .done / .failed
```

**Resume:** rerun the same `nohup` command. Completed outputs are not
redone (tested in a sandbox: second launch "Nothing to be done"). Jobs
interrupted mid-way are rerun (`--rerun-incomplete`). An interrupted hifiasm
leaves `.bin` files and refuses to restart until they are moved aside:

```bash
mkdir -p ../results/<LINE>/assembly/<method>/aside && mv ../results/<LINE>/assembly/<method>/*.bin $_
```

To stop everything: `pkill -f run_all_assemblies.sh; pkill -f "snakemake per_sample"`.
Snakemake stops its jobs on SIGTERM. Then resume as above.
