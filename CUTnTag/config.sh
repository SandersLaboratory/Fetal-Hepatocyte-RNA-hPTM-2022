# config.sh — pipeline configuration (Kraken, Anders's CUT&Tag reanalysis)
# Resolved location of the pipeline itself (must come first: used below)
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES_TSV="$SRC_DIR/samples.tsv"

# ===== Paths for this analysis =====
FASTQ_DIR=~/BrownRIH/fastq   # dir with *_R1_001.fastq.gz / *_R2_001.fastq.gz
OUTDIR=~/BrownRIH/cnt_output # all outputs go here
GENOME_DIR=~/BrownRIH/genome # rn7 fasta/gtf/index land here
THREADS=16
USE_SPIKEIN=false                    # true if CUT&Tag-IT E. coli spike-in was used
DEG_CSV="$SRC_DIR/lcm_degs.csv"   # columns: ensembl_id,gene_symbol,logFC,FDR (colony vs host)
# ======================

export FASTQ_DIR OUTDIR GENOME_DIR THREADS USE_SPIKEIN DEG_CSV SAMPLES_TSV SRC_DIR
mkdir -p "$OUTDIR"/{qc,trimmed,bam,peaks,bigwig,deseq2,integration,motifs,logs}
