# g2f_2024
Genomes to fields 2024 genome assemblies

Machine used for compute: ```cbsuxm01```  
Files stored on: ```cbsublfs1```   

## Assembly
Combine HiFi bams into single fastq.gz file:
```
nohup bam2fasta -c 1 \
-o EB2_PHBW8.hifi \
m84109_240806_024441_s4.hifi_reads.bc2027.bam \
m84109_240611_191050_s2.hifi_reads.bc2088.bam \
> bam2fasta.log 2>&1 &
```
Make directory for hifiasm run:
```
mkdir -p hifiasm_v0.25.0
cd hifiasm_v0.25.0
```
Run hifiasm:
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
