#!/bin/bash
# 00_setup.sh — reference genome, annotation, SEACR, chrom sizes
set -euo pipefail
source "$(dirname "$0")/../config.sh"

mkdir -p "$GENOME_DIR" && cd "$GENOME_DIR"

# --- mRatBN7.2 (rn7), Ensembl release 110 ---
FA=ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna
GTF=ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus
[ -f Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa ] || {
  wget -q "$FA/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
  gunzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
}
[ -f Rattus_norvegicus.mRatBN7.2.110.gtf ] || {
  wget -q "$GTF/Rattus_norvegicus.mRatBN7.2.110.gtf.gz"
  gunzip Rattus_norvegicus.mRatBN7.2.110.gtf.gz
}
[ -f rn7.1.bt2 ] || bowtie2-build --threads "$THREADS" \
  Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa rn7

samtools faidx Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
cut -f1,2 Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.fai > rn7.chrom.sizes

# --- SEACR (no-control CUT&Tag peak caller) ---
# NOTE: SEACR_1.3.sh is only a wrapper — it invokes the companion SEACR_1.3.R
# from its own directory. BOTH files must be present, or every call dies at
# the final step ("cannot open file '.../SEACR_1.3.R'") after the AUC work.
[ -f SEACR_1.3.sh ] || {
  wget -q https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.sh
  chmod +x SEACR_1.3.sh
}
[ -f SEACR_1.3.R ] || wget -q https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.R

# --- E. coli spike-in index (only used if USE_SPIKEIN=true) ---
if [ "$USE_SPIKEIN" = true ]; then
  [ -f ecoli.1.bt2 ] || {
    wget -q ftp.ensemblgenomes.org/pub/release-110/bacteria/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz
    gunzip -c Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz > ecoli.fa
    bowtie2-build --threads "$THREADS" ecoli.fa ecoli
  }
fi
# --- R / Bioconductor packages via BiocManager (not conda) ---
# MGB's conda-forge mirror does not carry Bioconductor-built R packages,
# so install them from Bioconductor directly. Idempotent: skips what's present.
#
# Compiler note: some Bioconductor C code (e.g. rtracklayer's bundled UCSC
# sources) predates C23, now the default in current conda-forge compilers.
# Force C17 for R package compilation as belt-and-braces.
mkdir -p ~/.R
grep -q "std=gnu17" ~/.R/Makevars 2>/dev/null || echo 'CFLAGS += -std=gnu17' >> ~/.R/Makevars
Rscript -e '
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
need <- c("DESeq2", "apeglm", "Rsubread", "GenomicRanges")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
cat("R packages OK\n")
'

echo "setup done: $GENOME_DIR"
