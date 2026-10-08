# g2f_2024
Genomes to fields 2024 genome assemblies

Machine used for compute: ```cbsuxm01```  
Files stored on: ```cbsublfs1```   

## Assembly
#### Combine HiFi bams into single fastq.gz file:
```
nohup bam2fasta -c 1 \
-o EB2_PHBW8.hifi \
m84109_240806_024441_s4.hifi_reads.bc2027.bam \
m84109_240611_191050_s2.hifi_reads.bc2088.bam \
> bam2fasta.log 2>&1 &
```
#### Make directory for hifiasm run:
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
    ../../EB2_PHBW8.hifi.fasta.gz \
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
    ../EB2_PHBW8.hifi.fasta.gz \
    > EB2_PHBW8.asm.stdout.log \
    2> EB2_PHBW8.asm.stderr.log &

echo $!
```
Convert hifiasm GFA file output into FASTA (from hifiasm Github)
```
awk '/^S/{print ">"$2;print $3}' \
    EB2_PHBW8.asm.bp.p_ctg.gfa \
    > EB2_PHBW8.asm.bp.p_ctg.fa
```
Run BUSCO
```
nohup singularity exec ~/busco.sif busco     -i EB2_PHBW8.asm.bp.p_ctg.fa     -o EB2_PHBW8_busco_liliopsida_odb10     -m genome     -l liliopsida_odb10     --download_path ~/busco_downloads     --offline     --skip_bbtools     --opt-out-run-stats     -c 40     > EB2_PHBW8.busco.stdout.log     2> EB2_PHBW8.busco.stderr.log &
```
