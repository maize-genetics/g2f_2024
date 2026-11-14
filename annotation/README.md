# Maize genome annotation

Scripts for gene and transposable element (TE) annotation of maize genome assemblies on a Slurm
cluster: gene prediction with Helixer, completeness with BUSCO, gene-model scores with PSAURON,
separation of the nuclear genome, syntenic orthologs to B73 with MCScan, and TE annotation with
EDTA and RepeatMasker.

Gene prediction and BUSCO are run on the whole genome; the later steps use the nuclear genome.

```
01 Helixer ─► 02 BUSCO ─► 03 PSAURON ─► 04 nuclear split ─┬─► 05 B73 synteny
   (whole genome)                                         └─► 06 TE annotation + LAI
```

| Step | Folder | What it does | Main tools |
|---|---|---|---|
| 01 | [`01_helixer/`](01_helixer/) | gene prediction on the whole genome | Helixer 0.3.3 (`land_plant_v0.3_a_0080`) |
| 02 | [`02_busco/`](02_busco/) | protein and CDS sequences; gene-space completeness | gffread, BUSCO 6.1.0 (`poales_odb12`) |
| 03 | [`03_psauron/`](03_psauron/) | scores every gene model and adds the score to the GFF3 | PSAURON 1.1.3 |
| 04 | [`04_nuclear_split/`](04_nuclear_split/) | labels scaffolds as nuclear, plastid or mitochondrion; writes the nuclear genome and its gene models | minimap2 |
| 05 | [`05_b73_synteny/`](05_b73_synteny/) | one-to-one syntenic orthologs to B73 (Zm00001eb.1) | MCScan (JCVI 1.6.8), LAST |
| 06 | [`06_te_annotation/`](06_te_annotation/) | TE annotation and LTR Assembly Index | EDTA 2.3.0 / 2.2.2, RepeatMasker 4.2.3, LTR_retriever 3.0.4 |

Each folder has its own README with the settings used, the commands to run, and the files written.

## Outputs

Everything is written under `PROJECT_DIR`:

| Folder | Contents |
|---|---|
| `inputs/` | the assemblies, one FASTA per unit, and their checks |
| `raw/<unit>/` | whole-genome Helixer GFF3 (also with PSAURON scores), proteins, CDS, gene labels |
| `busco/` | BUSCO results |
| `psauron/` | PSAURON scores and summary |
| `scaffold_classes/` | scaffold labels from the organelle screen |
| `nuclear_genome/<unit>/` | nuclear genome FASTA, GFF3 (with PSAURON scores), CDS, proteins, B73 orthologs |
| `b73_synteny/` | MCScan runs and ortholog tables |
| `te_annotation/` | TE annotation per unit (`runs/<unit>/results/`) and LAI |
