# Historical PHBW8 assemblies (manual runs, October 2026)

> **Not produced by the Snakemake workflow.** These commands were run by
> hand on `cbsuxm01` and recorded in the project README on 2026-10-07 and
> 2026-10-08, before the workflow existed. They are kept as the historical
> record of how the first PHBW8 assemblies were made and how we arrived at
> the standardized pipeline. To
> produce new assemblies, use the workflow (see the main `README.md`), not
> these commands.

The [original record](#original-record-verbatim) below is the project's
first GitHub README (commit `f809f9f`, 2026-10-08). It is copied byte for
byte. Paths, thread counts and comments have not been edited, even where
they are now obsolete or disagree with the run logs; the differences are
listed here instead.

Related records:

* `docs/legacy/PHBW8_legacy_runs.yaml`: what the hifiasm logs, GFAs and
  BUSCO JSONs of these runs show (verified 2026-10-10).
* `docs/audit_2026-10-10.md`: audit of the historical outputs.
* `docs/qc_PHBW8_imported_2026-10-11.md`: the first workflow QC of the
  historical HiFi+ONT `-l0` assembly, registered read-only.

## Historical steps and their workflow equivalents

| Step | Historical (manual) | Workflow |
|---|---|---|
| HiFi BAM → reads | `bam2fasta -c 1` (pbtk, version not recorded); two BAMs merged into one `EB2_PHBW8.hifi.fasta.gz` | `bam2fastq` (pbtk 3.5.0 container) per BAM, no merging; hifiasm receives every file |
| ONT reads | `fastq_pass/*.fastq.gz` extracted from two tar archives (flow cell PAQ81223), then `seqkit seq -m 50001` | `seqkit seq -m 50000` per input file (method `…_hifi-ont-ul50k_l0_v1`); for PHBW8 the historical pre-filtered file is used as supplied, after checking every read is ≥ 50,000 bp |
| hifiasm, defaults | `-t 100`, no other options | Not a workflow method (kept for comparison only) |
| hifiasm, `-l0` | HiFi only, `-l0` | Method `hifiasm-0.25.0_hifi_l0_v1` (`-l0`) |
| hifiasm, HiFi + ONT, `-l0` | `--ul PHBW8_ONT_gt50kb.fastq.gz -l0` | Method `hifiasm-0.25.0_hifi-ont-ul50k_l0_v1` (`--ul` + `-l0`) |
| GFA → FASTA | `awk '/^S/{print ">"$2;print $3}'` | Same `awk` conversion (`gfa_to_fasta`), then checked record by record against the GFA (`check_gfa_fasta`) |
| Contiguity | `seqkit stats -a` (version not recorded) | `seqkit stats` 2.13.0, plus per-contig table and validation |
| BUSCO | `~/busco.sif`, `-l liliopsida_odb10 -m genome --offline --skip_bbtools -c 40` | BUSCO 6.1.0 (`/programs/busco-6.1.0/busco.sif`), `--miniprot`, lineage `poales_odb12.2` pinned and checked |

hifiasm was `/programs/hifiasm-0.25.0/hifiasm`, which reports
`0.25.0-r726` in all three logs.

## Where the record differs from, or goes beyond, what was verified

* **Paths are obsolete.**
  * The reads were in `/workdir/naj49/g2f_assemblies_2024/` and are now in
    `data/PHBW8/hifi/` and `data/PHBW8/ont/`.
  * The run directories are now `hifiasm_v0.25.0/hifi_only/default`,
    `hifi_only/l0` and `ont_hifi`, not the `default/` layout below.
* **HiFi `-l0` threads.** The record shows `-t 100`, but the hifiasm log of
  that run shows `-t 120` with a relative input path
  (`docs/legacy/PHBW8_legacy_runs.yaml`). The log is authoritative.
* **Default run input path.** The log shows `../../EB2_PHBW8.hifi.fasta.gz`
  rather than the absolute path below.
* **ONT length filter.** `seqkit seq -m 50001` keeps reads ≥ 50,001 bp,
  i.e. strictly longer than 50 kb.
  * This is the only record of the filtering command. It was not captured
    in a log, and the tar archives have not been located.
  * The workflow's checks agree with it: the minimum length in the sampled
    reads is 50,012 bp.
  * `config/samples.yaml` still lists `filtering_command` and
    `length_boundary_convention` as unknown. Recording them from this
    document is a decision for the project lead.
* **HiFi conversion.** The source BAMs have not been located. The
  `bam2fasta` version and the checksums of its output were not recorded.
* **BUSCO.** The BUSCO JSONs show version 6.1.0 in mode `euk_genome_min`
  (miniprot 0.18-r281); the container's origin is otherwise unrecorded.
  Because `liliopsida_odb10` (n = 3,236) is a different marker set from the
  workflow's `poales_odb12.2` (n = 6,188), the scores are not directly
  comparable.
* **GFA → FASTA.** For the HiFi+ONT `-l0` assembly, the workflow verified
  the output of this command: 490 of 490 records are identical to the GFA
  segments.

## Historical results (from the logs and outputs)

| Run | Contigs | Total (bp) | N50 (bp) | BUSCO `liliopsida_odb10` |
|---|---:|---:|---:|---|
| HiFi, defaults | 1,236 | 2,280,334,098 | 121,136,081 | C 99.1 [S 82.5, D 16.6], F 0.6, M 0.2 |
| HiFi, `-l0` | 1,492 | 2,315,072,677 | 110,301,032 | C 99.1 [S 81.9, D 17.3], F 0.6, M 0.2 |
| HiFi + ONT, `-l0` | 490 | 2,266,382,187 | 121,136,088 | C 99.1 [S 82.2, D 17.0], F 0.6, M 0.2 |

Source: `docs/legacy/PHBW8_legacy_runs.yaml`, primary contigs (`p_ctg`).
The workflow QC of the HiFi+ONT `-l0` assembly with `poales_odb12.2` is in
`docs/qc_PHBW8_imported_2026-10-11.md`.

---

## Original record (verbatim)

# g2f_2024
Genomes to fields 2024 genome assemblies

Machine used for compute: ```cbsuxm01```  
Files stored on: ```cbsublfs1```   

## Preprocess data
#### Combine HiFi bams into single fastq.gz file:
```
nohup bam2fasta -c 1 \
-o EB2_PHBW8.hifi \
m84109_240806_024441_s4.hifi_reads.bc2027.bam \
m84109_240611_191050_s2.hifi_reads.bc2088.bam \
> bam2fasta.log 2>&1 &
```
#### Extract and combine nanopore FASTAs
```
set -o pipefail

(
  tar -xOf 240520-10474838_2_PAQ81223_a7f454ae_fastq.tar \
    --wildcards './fastq_pass/*.fastq.gz' &&

  tar -xOf 240520-10474838_2b_PAQ81223_10d1bef7_fastq.tar \
    --wildcards './fastq_pass/*.fastq.gz'

) | gzip -dc \
  | seqkit seq -m 50001 \
      -o PHBW8_ONT_gt50kb.fastq.gz
```

## Assembly
#### Make directory for hifiasm runs:
```
mkdir -p hifiasm_v0.25.0
cd hifiasm_v0.25.0
```
#### Run hifiasm:
Defaults:
```
mkdir -p default
cd default
```
```
nohup /programs/hifiasm-0.25.0/hifiasm \
    -o EB2_PHBW8.asm \
    -t 100 \
    /local/workdir/naj49/g2f_assemblies_2024/EB2_PHBW8.hifi.fasta.gz \
    > EB2_PHBW8.asm.stdout.log \
    2> EB2_PHBW8.asm.stderr.log &

echo $!
```
-l 0 (no duplication purging for inbred line):
```
nohup /programs/hifiasm-0.25.0/hifiasm \
    -o EB2_PHBW8.asm \
    -t 100 \
    -l0 \
    /local/workdir/naj49/g2f_assemblies_2024/EB2_PHBW8.hifi.fasta.gz \
    > EB2_PHBW8.asm.stdout.log \
    2> EB2_PHBW8.asm.stderr.log &

echo $!
```
With ONT l0
```
nohup /programs/hifiasm-0.25.0/hifiasm \
    -o PHBW8_HiFi_Nano_l0 \
    --ul /local/workdir/naj49/g2f_assemblies_2024/PHBW8_ONT_gt50kb.fastq.gz \
    -t 100 \
    -l0 \
    /local/workdir/naj49/g2f_assemblies_2024/EB2_PHBW8.hifi.fasta.gz \
    > PHBW8_HiFi_Nano_l0.stdout.log \
    2> PHBW8_HiFi_Nano_l0.stderr.log &

echo $!
```
Convert hifiasm GFA file output into FASTA (from hifiasm Github)
```
awk '/^S/{print ">"$2;print $3}' \
    PHBW8_HiFi_Nano_l0.bp.p_ctg.gfa \
    > PHBW8_HiFi_Nano_l0.bp.p_ctg.fa
```
seqkit stats
```
seqkit stats -a PHBW8_HiFi_Nano_l0.bp.p_ctg.fa \
> PHBW8_HiFi_Nano_l0_seqstats.txt
```
Run BUSCO
```
nohup singularity exec \
  --bind "$PWD:$PWD:rw" \
  ~/busco.sif busco \
  -i PHBW8_HiFi_Nano_l0.bp.p_ctg.fa \
  -o PHBW8_HiFi_Nano_l0_busco_liliopsida_odb10 \
  -m genome \
  -l liliopsida_odb10 \
  --download_path ~/busco_downloads \
  --offline \
  --skip_bbtools \
  --opt-out-run-stats \
  -c 40 \
  > PHBW8_HiFi_Nano_l0.busco.stdout.log \
  2> PHBW8_HiFi_Nano_l0.busco.stderr.log &

echo $!
```



## Future Dev Rough Ideas (could be complete nonsense)
- Rough potential structure for pipeline
- Keep modular for future
- Agent skill or similar for easy updating, and maybe usage by user
- Agent generates snakemake input file and run command for user to check over and run themself
- ?User sends agent some biological information to generate a snakemake file? Either input JSON or plain text paragraph on info? JSON probably easier, which is essentially bulleted list of ploidy, genome length, heterozygosity, etc? Or does pipeline calculate this for you and then does a run? Either way, user does the runs.
- 
```
genome-pipeline/
├── AGENTS.md
├── CLAUDE.md
├── README.md
├── config/
│   ├── species/
│   │   ├── zea_mays.yaml
│   │   └── arabidopsis_thaliana.yaml
│   ├── samples/
│   │   └── PHBW8.yaml
│   └── profiles/
│       └── cornell_biohpc.yaml
├── docs/
│   ├── decisions.md
│   └── methods.md
├── scripts/
│   └── assembly/
├── runs/
│   └── PHBW8/
│       ├── commands.sh
│       ├── manifest.yaml
│       └── results_summary.tsv
├── workflow/             # Later
│   ├── Snakefile
│   └── rules/
│       ├── assembly.smk
│       ├── qc.smk
│       ├── repeats.smk
│       ├── genes.smk
│       └── cres.smk
└── tests/
```
