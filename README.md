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
    EB2_PHBW8.asm.bp.p_ctg.gfa \
    > EB2_PHBW8.asm.bp.p_ctg.fa
```
seqkit stats
```
seqkit stats -a EB2_PHBW8.asm.bp.p_ctg.fa
```
Run BUSCO
```
nohup singularity exec ~/busco.sif busco     -i EB2_PHBW8.asm.bp.p_ctg.fa     -o EB2_PHBW8_busco_liliopsida_odb10     -m genome     -l liliopsida_odb10     --download_path ~/busco_downloads     --offline     --skip_bbtools     --opt-out-run-stats     -c 40     > EB2_PHBW8.busco.stdout.log     2> EB2_PHBW8.busco.stderr.log &
```



## Future Dev Rough Ideas (could be complete nonsense, I will know later when I think about it more)
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
