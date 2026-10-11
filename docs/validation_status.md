# Validation status of workflow components (2026-10-11)

What each status means:

* **Real-tool (synthetic):** the real BioHPC tool was run through the workflow
  on tiny synthetic data (`tests/test_realtools.py`;
  `docs/realtools_validation_2026-10-11.md`).
* **Fixture-tested:** logic exercised on synthetic data or on small excerpts
  of real PHBW8 outputs.
* **Stub-tested:** the plumbing runs end to end, but the tool is a stub
  (`tests/stubs/`).
* **Awaiting real-data validation:** not yet run at maize scale.

**Real-data validation so far:** QC of the historical PHBW8 HiFi+ONT `-l0`
assembly (`imported_qc`, 2026-10-11; `docs/qc_PHBW8_imported_2026-10-11.md`).
No PHBW8 sequencing reads have been processed yet. A successful dry run is
not operational validation.

| Component | Status | Still to validate on real data |
|---|---|---|
| Input discovery, duplicate/derived detection | Fixture-tested (20 tests). Run read-only on the real PHBW8 files: passes, with the approved acknowledgement | — |
| Empty line directories (`NO_INPUTS` warning) and explicit requests for lines without reads (`REQUESTED_SAMPLE_NO_INPUTS` error) | Fixture-tested. Real layout (2026-10-11): 8 empty lines skipped with warnings; `only_samples=3AZA1` exits 1 | — |
| Method selection, ONT policy, pre-filtered ONT used as supplied | Fixture-tested and stub-tested. PHBW8 dry run selects HiFi+ONT and uses the supplied ONT file directly | Full-file check that every PHBW8 ONT read is ≥ 50,000 bp |
| HiFi BAM → FASTQ (pbtk 3.5.0) | Real-tool (synthetic BAM) | Real Revio BAMs (PHBW8 BAMs not located) |
| HiFi/ONT read QC, checksums, movie census, ONT ID check | Real-tool (synthetic) | Run time and memory on 28 GB of HiFi |
| meryl count / histogram | Real-tool (synthetic) | Memory and time at 84 Gb |
| k-mer screen (diagnostic) | Fixture-tested plus one real meryl spectrum; thresholds **uncalibrated** | PHBW8 spectrum, then calibration |
| hifiasm command construction and outputs | Real-tool (synthetic, `-l0`, `--ul`, and a comma-separated multi-file `--ul` in which all files were consumed); PHBW8 command rendered | A full-scale run (not authorised) |
| GFA → FASTA and consistency check | Real-tool (synthetic). **Real data:** PHBW8 legacy 2.3 GB GFA vs FASTA, 490/490 identical | — |
| Assembly statistics | **Real data:** PHBW8 legacy, reproduces the audit values exactly | — |
| BUSCO 6.1.0 | **Real data:** PHBW8 legacy with the pinned `poales_odb12.2` (728 s at 32 threads); pins enforced | Peak memory (recorded from now on) |
| Merqury 1.3 | Real-tool (synthetic) | A maize-scale run (needs the full meryl database) |
| Provenance, version checks, write protection | Stub- and real-tool-tested | — |
| Rerun policy (science/method/code/resource changes, failed-job resume, input changes) | Stub-tested (`test_workflow.py`) | — |
| Imported historical assemblies (read-only QC) | **Real data:** PHBW8 legacy HiFi+ONT `-l0`; legacy files verified unchanged | Merqury (needs the read k-mer database) |
| Summaries | Stub- and real-tool-tested | — |
