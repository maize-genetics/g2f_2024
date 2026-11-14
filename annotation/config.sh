# config.sh: settings for your system. Fill these in before running anything.
# Every script reads this file; submit jobs from the repository's top folder so they can find it:
#     cd maize-genome-annotation && sbatch -A "$SLURM_ACCOUNT" 01_helixer/helixer_array.sbatch
# (source config.sh first so $SLURM_ACCOUNT etc. are set in your shell)

# --- project -----------------------------------------------------------------
PROJECT_DIR=""          # where inputs and outputs are written (needs a few hundred GB)
ASSEMBLIES_TSV=""       # table of genome assemblies to annotate (see assemblies_example.tsv)

# --- Slurm (used on the sbatch command line, see README) ----------------------
SLURM_ACCOUNT=""        # your Slurm account
SLURM_PARTITION=""      # partition / queue
EMAIL=""                # address for job notifications (leave empty for none)

# --- software ----------------------------------------------------------------
HELIXER=""              # command that runs Helixer 0.3.3, e.g. "Helixer.py" or "apptainer exec --nv helixer_v0.3.3.sif Helixer.py"
MINIMAP2=""             # minimap2 executable (or leave empty if minimap2 is on your PATH)
GFFREAD=""              # gffread executable (or leave empty if gffread is on your PATH)
BUSCO=""                # busco 6.1.0 executable (or leave empty if busco is on your PATH)
BUSCO_LINEAGE=""        # path to the downloaded poales_odb12 lineage folder
JCVI_PYTHON=""          # python that has jcvi 1.6.8 (LAST and GLPK must be on the PATH too), or leave empty for "python"
PSAURON=""              # command that runs PSAURON 1.1.3, e.g. "psauron" or "apptainer exec psauron.sif psauron"
EDTA230_SIF=""          # EDTA 2.3.0 image
EDTA222_SIF=""          # EDTA 2.2.2 image (TIR step)
REPEATMASKER_SIF=""     # RepeatMasker 4.2.3 image

# --- reference data ------------------------------------------------------------
ORGANELLE_FASTA=""      # maize plastid + mitochondrion (NC_001666, NC_007982), for the scaffold screen
B73_GFF3=""             # B73 RefGen_v5 annotation, Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3(.gz), for MCScan
B73_CDS=""              # B73 Zm00001eb.1 CDS FASTA(.gz), for MCScan
B73_CANONICAL=""        # B73 Zm00001eb.1 canonical transcript list(.gz), one ID per line, for MCScan
