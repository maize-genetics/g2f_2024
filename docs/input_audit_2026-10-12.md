# Sequencing input audit (2026-10-12)

Read-only audit of `data/` against `data/g2f_9genotypes_raw_reads_manifest.tsv`.
Nothing was moved, renamed or deleted.

**Verdict: the files are correctly organized.** Every manifest file for the
eight non-PHBW8 genotypes is present in `data/<genotype>/{hifi,ont}/` at
exactly its manifest size. No file is misplaced, duplicated, truncated or
unlisted. PHBW8 holds only its previously prepared derivatives, as intended.
The manifest agrees with the summary table and with the sample names
embedded in the files (HiFi `@RG SM`, ONT `sample_id`). One ONT run is a
near-empty aborted run that needs a decision before assembly (below).

## Manifest format

* `g2f_9genotypes_raw_reads_manifest.tsv`: tab-delimited, LF line endings,
  header + 38 rows, 7 columns. Every row has 7 fields; no quoting,
  whitespace padding or trailing delimiters. Columns: `genotype`, `sample`
  (EB1-EB10; no EB6), `technology` (HiFi | ONT), `dataset_id`, `source_type`
  (`pacbio_hifi_bam` | `ont_fastq_tar`), `source_path` (on cbsublfs1),
  `size_bytes`.
* `g2f_9genotypes_summary.tsv`: tab-delimited, **CRLF** line endings
  (harmless; parse with stripping). It agrees with the manifest for all nine
  genotypes (file counts and HiFi byte totals).
* `data/scp_reads.sh`: the transfer script (rsync `--partial
  --append-verify`, size check after each file; PHBW8 skipped by design).
  No transfer was running at audit time.
* Nine genotypes: FBLL (EB1), PHBW8 (EB2), 3AZA1 (EB3; the directory is
  `3AZA1`, not "3AZ1"), PHB47 (EB4), PHP02 (EB5), GEMS-0113 (EB7),
  GEMS-0061 (EB8), GEMS-0219 (EB9), GEMS-0227 (EB10).

## Audit table

| Genotype | Sample | Tech | Expected (manifest) | Observed | Observed bytes | Complete | Manifest agreement | Concerns |
|---|---|---|---|---|---|---|---|---|
| FBLL | EB1 | HiFi | 3 BAM (58.0 GB) | 3 | 58.0 GB | yes | names, sizes, movies, `SM:Maize EB1-FBLL` | — |
| FBLL | EB1 | ONT | 2 tar (16.4 GB) | 2 | 16.4 GB | yes | names, sizes, run IDs, `sample_id=10474838_1/_1b` | — |
| PHBW8 | EB2 | HiFi | 2 BAM (50.8 GB) | 1 FASTA.gz | 29.2 GB | n/a | derived file contains exactly the 2 manifest movies (1,153,903 + 3,894,183 reads) | raw BAMs not on this server (by design) |
| PHBW8 | EB2 | ONT | 2 tar (21.2 GB) | 1 FASTQ.gz (≥ 50 kb) | 0.72 GB | n/a | 12,416 reads from runs `10474838_2` (10,382) and `_2b` (2,034), flow cell PAQ81223 = the 2 manifest runs | pre-filtered (acknowledged in `samples.yaml`) |
| 3AZA1 | EB3 | HiFi | 3 BAM (84.6 GB) | 3 | 84.6 GB | yes | yes (`SM:Maize EB3-3AZA1`) | largest HiFi set (~61x) |
| 3AZA1 | EB3 | ONT | 1 tar (11.8 GB) | 1 | 11.8 GB | yes | yes (`10474838_3`) | — |
| PHB47 | EB4 | HiFi | 2 BAM (41.2 GB) | 2 | 41.2 GB | yes | yes (`SM:Maize EB4-PHB47`) | — |
| PHB47 | EB4 | ONT | 3 tar (7.9 GB) | 3 | 7.9 GB | yes | yes (`10474838_4/_4b/_4c`) | `_4b` and `_4c` share flow cell PAQ86263 (reused flow cell, different runs; not duplicates) |
| PHP02 | EB5 | HiFi | 1 BAM (46.0 GB) | 1 | 46.0 GB | yes | yes (`SM:Maize EB5 - PHP02`) | single HiFi run: only split-half fingerprint controls |
| PHP02 | EB5 | ONT | 1 tar (22.8 GB) | 1 | 22.8 GB | yes | yes (`10474838_5`) | basecalled with **SUP** v4.2.0 (all other runs HAC v4.2.0) |
| GEMS-0113 | EB7 | HiFi | 1 BAM (43.5 GB) | 1 | 43.5 GB | yes | yes (`SM:Maize EB7 - GEMS-0113`) | single HiFi run |
| GEMS-0113 | EB7 | ONT | 2 tar (5.4 GB) | 2 | 5.4 GB | yes | yes (`10474838_7/_7b`) | — |
| GEMS-0061 | EB8 | HiFi | 1 BAM (45.7 GB) | 1 | 45.7 GB | yes | yes (`SM:Maize EB8 - GEMS-0061`) | single HiFi run; barcode bc2009 also used for a GEMS-0219 run (different SMRT cell: normal) |
| GEMS-0061 | EB8 | ONT | 3 tar (11.5 GB) | 3 | 11.5 GB | yes | yes (`10474838_8/_8b/_8c`) | `_8`/`_8b` share flow cell PAQ87161 (reuse) |
| GEMS-0219 | EB9 | HiFi | 3 BAM (48.1 GB) | 3 | 48.1 GB | yes | yes (`SM:Maize EB9 GEMS-0219`) | — |
| GEMS-0219 | EB9 | ONT | 3 tar (11.5 GB) | 3 | 11.5 GB | yes | yes (`10474838_9/_9b/_9c`) | **`10474838_9` is a 204,800-byte tar: 1 member, 52 reads, 180 kb, longest 21.5 kb** (aborted run; complete, not truncated) |
| GEMS-0227 | EB10 | HiFi | 2 BAM (58.6 GB) | 2 | 58.6 GB | yes | yes (`SM:Maize EB10 GEMS-0227`) | manifest dataset ID `Maize_EB10-GEMS-227` (source folder name) |
| GEMS-0227 | EB10 | ONT | 3 tar (17.5 GB) | 3 | 17.5 GB | yes | yes (`10474838_10/_10b/_10c`) | `_10b`/`_10c` share flow cell PAQ86663 (reuse) |

Totals: HiFi 16 BAM files, 425.8 GB (+ PHBW8 FASTA 29.2 GB). ONT 18 tar
files, 104.8 GB (+ PHBW8 0.72 GB). No technology is missing for any
genotype: all nine have HiFi and ONT.

## Checks performed (cheapest first)

1. **Inventory:** recursive listing of every genotype directory. Only
   `hifi/` and `ont/`, no nested directories or stray files.
2. **Names and sizes vs manifest:** 34/34 present with byte-exact sizes,
   0 unlisted (except PHBW8's two derivatives), no duplicate basenames or
   dataset IDs. Exact size equality after rsync `--append-verify` is strong
   evidence of complete transfers. No source checksums exist in the
   manifest, and full-file MD5s are computed by the workflow before
   assembly anyway (`checksum` rule), so no extra checksum pass was made.
3. **HiFi BAMs (16):** `samtools quickcheck -u` (the `-u` flag accepts
   unmapped BAMs, which legitimately have no `@SQ` targets; without it,
   quickcheck's "no targets" message is not evidence of corruption), the
   BGZF EOF marker (present in all), `@SQ` count = 0 (unaligned), `@RG
   PL:PACBIO`, `READTYPE=CCS`, `PU` = the movie in the file name, `SM` =
   the expected EB number and genotype in every file. No `.pbi` indexes
   were transferred, so the workflow builds them in `results/`.
4. **ONT tars (18):** full header listing (`tar tvf`; exits 0 = archive
   structure intact to the end). Each contains only
   `./fastq_pass/*.fastq.gz` (no `fastq_fail`, no summaries). Every member
   name carries the flow cell and protocol-run hash of its tar. Batch
   numbers are contiguous (0..N-1, no gaps). The **last member of every tar
   decompresses completely** (4,000 reads each; 52 for `10474838_9`).
5. **ONT technology and identity:** read headers have ONT UUIDs, `runid`
   = the acquisition hash in the member names, `flow_cell_id` = the tar's
   flow cell, `sample_id` = the manifest dataset ID, and
   `basecall_model_version_id` = R10.4.1 400 bps HAC v4.2.0 (PHP02: SUP).
   Mean read quality ~Q14 (PHB47 `_4c`: AvgQual 13.8, N50 10.8 kb).
6. **Workflow discovery** (`g2f_tool.py discover`, read-only) accepts all
   nine genotypes with no issues once ONT tar support was added. Before
   that, all 18 tars were rejected as `UNKNOWN_FILE`.

Not done (unnecessary or deferred to the workflow): full-file MD5 and gzip
integrity of every member (the workflow reads every byte during
checksumming and read QC and fails on any error), and full BAM record
decoding (done by `bam_record_count` and `bam2fastq`).

## Required decision

* **GEMS-0219 `10474838_9_PAQ86658_8f1582f4_fastq.tar`**: complete but
  essentially empty (52 reads, none ≥ 50 kb). The ONT validation gate will
  stop GEMS-0219 ("no reads ≥ 50000 bp; human review required").
  Recommended smallest safe correction: **do not move the file**; record
  an `exclude_inputs` entry in `config/samples.yaml` (snippet in
  `docs/full_scale_run_plan.md`). The genotype's other two ONT runs
  (11.5 GB) are unaffected.

## Observations (no action needed)

* Sample identity metadata (BAM `SM`, ONT `sample_id`, manifest) is
  internally consistent for every file, but it is still metadata. The
  sequence-based check is the fingerprinting report.
* PHBW8's raw BAMs and tars stay on the source server; its derivatives are
  traceable to exactly the manifest runs (above).
