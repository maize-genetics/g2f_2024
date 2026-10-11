# Validation status of workflow components (2026-10-11)

What each status means:

* **Fixture-tested:** exercised by `tests/` on synthetic data, or on small
  excerpts of real PHBW8 outputs.
* **Stub-tested:** the plumbing (commands, file wiring, validation gates and
  provenance) runs end to end, but the real tool is replaced by a stub in
  `tests/stubs/`.
* **Awaiting real-data validation:** not yet run on real data with the real
  tool.

A Snakemake dry run that succeeds does not show that the workflow works on
real data.

| Component | Implementation | Status |
|---|---|---|
| Input discovery, format checks, duplicate/derived detection | `scripts/g2f/discovery.py` | Fixture-tested (20 tests). Also run read-only on the real PHBW8 files: halts correctly on the pre-filtered ONT file |
| Method selection | `discovery.select_method` | Fixture-tested |
| Config loading and validation | `scripts/g2f/config.py` | Fixture-tested |
| HiFi BAM → FASTQ (pbtk `bam2fastq`) | `rules/hifi.smk` | Stub-tested (stub uses `samtools fastq`). **Real pbtk not yet run**; its container still has to be pulled |
| BAM record-count check | samtools 1.20 | Fixture-tested with a real tool on a tiny BAM |
| Read QC (seqkit `fx2tab` + aggregator) | `readstats.py` | Fixture-tested; real seqkit in the end-to-end test |
| ONT raw QC, filtering, merging, validation | `rules/ont.smk`, `validate.py` | Fixture-tested; real seqkit in the end-to-end test |
| MD5 checksums | `md5sum` | Fixture-tested (real tool) |
| meryl count / histogram | `rules/kmer.smk` | Stub-tested. **Real meryl in the Merqury 1.3 container not yet run** |
| k-mer screen and gate | `kmer.py` | Fixture-tested on synthetic spectra; **thresholds uncalibrated** |
| hifiasm invocation and version check | `rules/assembly.smk` | Stub-tested (command line, `--ul`, `-l0`, version check, stale-`.bin` guard, protected outputs) |
| GFA → FASTA and consistency check | awk + `gfa.py` | Fixture-tested; real awk in the end-to-end test |
| Assembly statistics | seqkit `stats`/`fx2tab` | Real tool in the end-to-end test; parser checked on the real PHBW8 seqkit output |
| BUSCO | `rules/qc.smk` | Stub-tested; parser checked on a real BUSCO 6.1.0 JSON. **Real BUSCO container not run** |
| Merqury | `rules/qc.smk` | Stub-tested; output format taken from Merqury documentation. **Not run**; thread control inside the container unverified |
| hifiasm log parsing | `validate.parse_hifiasm_log` | Fixture-tested on a real PHBW8 log excerpt |
| Provenance records | `provenance.py` | Stub-tested (contents checked in the end-to-end test) |
| Resumability, input-change refusal, write protection | Snakemake plus a parse-time check | Tested in the end-to-end test |
| Per-sample and project summary | `summary.py` | Stub-tested |

Before PHBW8 production, real-data validation must cover:

* pbtk version reporting;
* meryl memory and thread use;
* Merqury output names and threads;
* BUSCO's `--miniprot` / `--skip_bbtools` flags and output layout;
* hifiasm output naming at full scale.
