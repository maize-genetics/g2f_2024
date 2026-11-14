# 05 · B73 syntenic orthologs (nuclear genome)

Finds one-to-one syntenic orthologs between the B73 reference annotation and the nuclear
genes of each assembly with MCScan, and combines the results across assemblies.

| Script | What it does |
|---|---|
| `b73_prep.py` | writes the MCScan input files (`.bed` and `.cds`) for B73 and for each assembly |
| `b73_mcscan.sbatch` | runs MCScan for B73 against each assembly, one Slurm array task per assembly, then `b73_table.py` |
| `b73_table.py` | turns one assembly's MCScan result into ortholog tables (called by `b73_mcscan.sbatch`) |
| `b73_combined.py` | one table of all assemblies, the number of B73 genes found in each number of assemblies and lines, and a chromosome check |

## Settings

- **B73 reference:** Zm-B73-REFERENCE-NAM-5.0, annotation Zm00001eb.1, canonical transcripts only
  (one per gene), from MaizeGDB. Genes on organelle sequences are left out. Set `B73_GFF3`,
  `B73_CDS` and `B73_CANONICAL` in `config.sh`.
- **Assemblies:** the nuclear Helixer genes from step 04 (one mRNA per gene).
- **MCScan (JCVI v1.6.8)**, `jcvi.compara.catalog ortholog` on CDS: LAST alignment, C-score
  filter 0.99 (reciprocal best matches), synteny blocks, liftover of missed anchors, and a quota
  1:1 screen (`--quota=1:1`). The quota step needs the GLPK solver (`glpsol`).
- Inside a block a gene can still have more than one candidate partner. `b73_table.py` keeps one
  partner per gene (synteny anchors before lifted pairs, then higher score) and lists the others in
  `<unit>.B73_orthologs.dropped.tsv`. The best match is chosen by ranking all candidate pairs on
  these criteria and keeping a pair only if neither gene has already been paired.
- Set `JCVI_PYTHON` in `config.sh` to the python that has jcvi (LAST and GLPK must be on the PATH).
  Each assembly used 8 CPUs and 32 GB memory.

## Run

From the repository folder, after step 04, with the jcvi environment active:

```bash
source config.sh && mkdir -p logs
python3 05_b73_synteny/b73_prep.py --project "$PROJECT_DIR" \
    --b73-gff3 "$B73_GFF3" --b73-cds "$B73_CDS" --b73-canonical "$B73_CANONICAL"
N=$(wc -l < "$PROJECT_DIR/units.tsv")
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=1 05_b73_synteny/b73_mcscan.sbatch       # test one assembly
sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" --array=2-$N 05_b73_synteny/b73_mcscan.sbatch    # the rest
# when all have finished:
python3 05_b73_synteny/b73_combined.py --project "$PROJECT_DIR"
```

## Inputs and outputs

| | Path (under `PROJECT_DIR`) |
|---|---|
| nuclear genes (input) | `nuclear_genome/<unit>/<unit>.nuclear.helixer.gff3`, `.cds.fa` |
| MCScan inputs | `b73_synteny/shared/B73.bed`, `B73.cds`; `b73_synteny/runs/<unit>/<unit>.bed`, `<unit>.cds` |
| MCScan result | `b73_synteny/runs/<unit>/B73.<unit>.lifted.1x1.anchors` |
| ortholog pairs | `b73_synteny/runs/<unit>/<unit>.B73_orthologs.tsv` (B73 gene, assembly gene), also copied to `nuclear_genome/<unit>/` |
| pairs with details | `b73_synteny/runs/<unit>/<unit>.B73_orthologs.full.tsv` (transcripts, positions, block, anchor type, score) |
| per-assembly summary | `b73_synteny/runs/<unit>/<unit>.B73_orthologs.summary.tsv` |
| all assemblies | `b73_synteny/B73_orthologs_combined.tsv` (one row per B73 gene), also copied to `nuclear_genome/` |
| sharing | `b73_synteny/B73_orthologs_sharing.tsv` (B73 genes found in 0…N assemblies and lines) |
| chromosome check | `b73_synteny/chromosome_check.tsv`, `chromosome_check_discordant.tsv` |

In the sharing and combined tables, the two haplotypes of a line (`<line>_hap1`, `<line>_hap2`)
count as one line. The chromosome check reports, per assembly, how many ortholog pairs lie on the
same chromosome in B73 and in the assembly. When the assemblies were scaffolded against B73, it
measures agreement between scaffolding and synteny rather than an independent test.
