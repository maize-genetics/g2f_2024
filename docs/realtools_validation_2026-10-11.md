# Real-tool compatibility validation (2026-10-11)

This records which real BioHPC tools were run, on what data, and what we
found.

* **Data:** tiny synthetic inputs only, built by `tests/realtools_data.py`:
  * a random 300 kb genome;
  * 600 error-free HiFi-like reads of about 15 kb (~30×);
  * 35 ONT-like reads: 10 of about 20 kb and 25 of 60–90 kb.
* **Not run:** no PHBW8 read data and no production assembly.
* **Repeat the test:**

  ```bash
  G2F_REALTOOLS=1 bash tests/run_tests.sh test_realtools   # ~75 s, 2 threads
  ```

  This runs the complete Snakemake workflow, from discovery to the summary
  table, with the real tools named in `config/tools/biohpc.yaml`.

## Findings

| Tool (as installed) | Checked | Finding | Change made |
|---|---|---|---|
| **pbtk 3.5.0** (`bam2fastq`, `pbindex`). Container pulled by digest `sha256:ab2d8dcf…e0ab` (19 MB) to `../containers/pbtk-3.5.0.sif` | Version, `--help`, conversion of a 100-read PacBio-style BAM | `bam2fastq` **requires a `.pbi` index** next to the BAM and exits 1 without one. With an unusable index it **exits 0 and writes an empty FASTQ**. `-j` defaults to all CPUs | The rule symlinks the BAM into `results/…/<id>.pbtk/`, reuses an existing `.pbi` or builds one there with `pbindex` (`data/` untouched), and passes `-j {threads}`. Comparing the BAM record count with the converted read count, which was already in place, catches the silent empty output |
| **meryl 1.3** (inside the Merqury 1.3 container) | `count`, `histogram`, `threads=`, `memory=` | Works. `threads=` and `memory=` are honoured (e.g. "Estimated to require 39 MB memory out of 2048 MB allowed"). Uses meryl 1.3, not BioHPC's native meryl 1.4.1, so the database format matches Merqury | Container meryl kept; `OMP_NUM_THREADS` also set |
| **Merqury 1.3** (container) | Run on a real meryl database | **`merqury.sh` fails unless `MERQURY` is set**; the container does not set it. Merqury's internal meryl calls take no thread option but honour `OMP_NUM_THREADS` (2 threads observed on a 256-CPU host). It symlinks inputs into the working directory, so the assembly directory is never written to. Output formats confirmed: `<out>.qv` (name, asm-only k-mers, total, QV, error rate) and `<out>.completeness.stats` (name, `all`, solid in asm, solid total, %). The assembly name is the FASTA basename minus `.fa`. QV can be `inf` | `--env MERQURY=/opt/conda/share/merqury` added to the tools config; `OMP_NUM_THREADS` set. The parser was already correct and accepts `inf` |
| **BUSCO 6.1.0** (container) | `--version`, `--help`, source code, a full offline run with miniprot on the synthetic genome against the installed `poales_odb12` | Flags confirmed (`--miniprot`, `--skip_bbtools`, `--offline`, `-f`, `--download`). The default `datasets_version` is `odb12.2`. **With a local lineage directory, the dataset version comes from the `name=` line of its `dataset.cfg`**: the installed `poales_odb12` (OrthoDB 12.1, created 2026-05-22, 6,282 BUSCOs) runs as `odb12`, *not* `odb12.2`. **`poales_odb12.2` is not installed.** BUSCO prints a usage-statistics notice unless `--opt-out-run-stats` is given (the historical runs used it). Output layout and JSON confirmed; the parser works on the real output | Lineage pinned to `poales_odb12.2` with an identity check (see D7). `--opt-out-run-stats` added. A missing-lineage warning is printed at DAG build |
| **hifiasm 0.25.0-r726** | `--version`, `-h`, real runs with `-l0` and with `-l0 --ul` on synthetic data | `-l0` writes `*.bp.p_ctg.gfa` and `*.bp.a_ctg.gfa`. **`a_ctg` can be an empty (0-byte) file** for homozygous input. `--ul` accepts comma-separated files. On perfectly homozygous synthetic reads hifiasm still printed `peak_het: 26`, so `peak_het` lines are not evidence of heterozygosity | Empty `a_ctg` is accepted by validation. ONT files go to `--ul` as a comma-separated list, so no merged copy is written. The PHBW8 command rendered by the dry run equals the historical one except `-t` |
| **seqkit 2.13.0**, **samtools 1.20** | Used in all tests | `seqkit stats` on an empty FASTA reports `GC(%) 0.00` | GC range check skipped for empty contig sets |
| **Snakemake 9.26.1** (`module load python/3.12.7`) | Profile, DAG, resource requests, `ancient()`, rerun behaviour | Works. Under `--rerun-triggers mtime`, a changed input *set* is detected ("set of input files changed") but not acted on while outputs exist | Validation identity and per-sample decision IDs are content-addressed into marker and summary names (D8) |

## Result of the real-tool end-to-end run (synthetic data)

| Metric | Value |
|---|---|
| Method selected | `hifiasm-0.25.0_hifi-ont-ul50k_l0_v1` |
| HiFi reads | 600 (500 from FASTQ + 100 converted by pbtk) |
| ONT reads | 35 supplied, 25 retained (≥ 50 kb) |
| Contigs | 1 (299,228 bp, out of a 300,000 bp genome) |
| a_ctg | Empty |
| Merqury | QV `inf`, k-mer completeness 100% |
| BUSCO (`poales_odb12`) | 0% complete, as expected for random sequence |

**k-mer screen v1 false positive.** The first version of the k-mer screen
flagged this perfectly homozygous dataset. The cause was sampling noise in a
sparse histogram, which produced unsmoothed local maxima; the height ratio
was only 0.04. The screen now smooths the spectrum and requires secondary
maxima to reach 5% of the main peak. The real spectrum is kept as a
regression fixture (`tests/fixtures/meryl_k21_synthetic_homozygous_300kb.hist`).

## Follow-up checks (2026-10-11, later the same day)

### `poales_odb12.2` installed and verified

* **Download:** `busco --download_path ~/busco_downloads --download
  poales_odb12.2`, run from the BUSCO 6.1.0 container.
  * The first attempt failed with "Read-only file system" because the
    working directory was not bound into the container. It succeeded with
    `--bind /workdir,/local/workdir`.
  * BUSCO's downloader checks the archive MD5 (`c2cf3a2db7a23f9b8c0d3d9af54c5132`)
    against `file_versions.tsv` and stops on a mismatch, so the successful
    download means it matched.
  * The four existing lineage directories were left unchanged.
* **`dataset.cfg` contents:** `name=poales_odb12.2`, `OrthoDB_version=12.2`,
  `dataset_version=01`, `creation_date=2026-05-13`, `number_of_BUSCOs=6188`,
  `species=maize`.
  * The directory holds 6,188 HMM profiles and 6,188 score cutoffs, matching
    `number_of_BUSCOs`.
  * All of these values are pinned in `science.yaml`.
* **Compatibility:** an offline BUSCO 6.1.0 run with miniprot on the
  synthetic genome reported `datasets_version: odb12.2` and n = 6188, and
  passed every pin check.

### hifiasm consumes every file of a comma-separated `--ul` list

* **Setup:** two disjoint synthetic ONT FASTQs, A (10 reads) and B (15 reads),
  with no shared read IDs.
* **Result:** hifiasm 0.25.0-r726 reported `# passed UL reads` as:

  | `--ul` value | Passed UL reads |
  |---|---|
  | `A` | 10 |
  | `B` | 15 |
  | `A,B` | 25 |
  | `B,A` | 25 |

* **Conclusion:** both files are read; it is not just a clean exit.
* **Kept as a test:** `HifiasmMultipleUlFilesTest` in
  `tests/test_realtools.py`.
* **Recorded in summaries:** the passed-UL count from each hifiasm log, as
  `hifiasm_ul_passed_reads`. For example, the legacy PHBW8 HiFi+ONT run
  reports 12,416.

## Not validated by these tests

* Behaviour and resource use at maize scale (2.3 Gb, 84 Gb of HiFi).
* `poales_odb12.2` itself, which still has to be downloaded.
* pbtk on real Revio BAMs; the PHBW8 source BAMs have not been located.
* Merqury R plotting (not required for the parsed outputs).
