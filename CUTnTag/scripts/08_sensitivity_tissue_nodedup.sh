#!/usr/bin/env bash
# 08_sensitivity_tissue_nodedup.sh — SENSITIVITY CHECK (tissue only).
#
# Question: does strict duplicate removal drive the null tissue result?
# The main pipeline's dedup BAMs were made with
#   samtools view -f 2 -F 1280 -q 10    (1280 = duplicate 1024 + secondary 256)
# which discards 62–89% of tissue fragments. This script rebuilds the tissue
# BAMs from the .markdup.bam files with the IDENTICAL filters except
# duplicates retained (-F 256), then the companion R script recounts the
# existing consensus atlases, re-runs the TA-vs-TF DESeq2 contrast, and
# recomputes the cells<->tissue concordance summaries. Main pipeline
# outputs are NOT overwritten; everything lands in $OUTDIR/sensitivity/.
#
# Standalone (not part of run_all.sh):  bash scripts/08_sensitivity_tissue_nodedup.sh
set -euo pipefail
source "$(dirname "$0")/../config.sh"
# standalone script: activate the pipeline env (run_sbatch.sh normally does
# this for run_all.sh; without it Rscript resolves to the base R, which has
# no Rsubread/DESeq2)
source /homes2/aohman/mambaforge/etc/profile.d/conda.sh
conda activate cutntag
CHROMS=$(cut -f1 "$GENOME_DIR/rn7.chrom.sizes" | grep -v -i -E '^(MT|chrM)' | tr '\n' ' ')
mkdir -p "$OUTDIR/sensitivity/bam" "$OUTDIR/logs"
for sample in $(awk -F'\t' 'NR>1 && $6=="tissue" {print $1}' "$SAMPLES_TSV"); do
  B="$OUTDIR/bam/$sample"
  O="$OUTDIR/sensitivity/bam/$sample.nodedup.bam"
  if [ ! -s "$O" ]; then
    samtools view -b -@ "$THREADS" -f 2 -F 256 -q 10 -o "$O" "$B.markdup.bam" $CHROMS
    samtools index "$O"
  fi
done
Rscript "$SRC_DIR/scripts/08_sensitivity_tissue_nodedup.R" 2>&1 | tee "$OUTDIR/logs/08_sensitivity.log"
echo "08 sensitivity done"
