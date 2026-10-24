#!/bin/bash
# 06_metagenes.sh — deeptools signal profiles around TSS for DEG gene sets, per mark
set -euo pipefail
source "$(dirname "$0")/../config.sh"

MG="$OUTDIR/integration/metagenes"; mkdir -p "$MG"
GS="$OUTDIR/integration/genesets"

for mark in 2m1 2m3 3m1 3m2; do
  # bigWigs for this mark only
  BWS=$(tail -n +2 "$SAMPLES_TSV" | awk -F'\t' -v m="$mark" '$4==m{print $1}' \
        | while read -r s; do echo -n "$OUTDIR/bigwig/${s}.bw "; done)
  LABELS=$(tail -n +2 "$SAMPLES_TSV" | awk -F'\t' -v m="$mark" '$4==m{print $1"_"$7"_"$8}' \
        | tr '\n' ' ')
  for gs in colony_up host_up leading_edge background; do
    TSS="$GS/${gs}.tss.bed"
    [ -s "$TSS" ] || { echo "skip empty $gs"; continue; }
    M="$MG/${mark}_${gs}.gz"
    computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
      -R "$TSS" -S $BWS --samplesLabel $LABELS \
      --skipZeros -p "$THREADS" -o "$M" 2>>"$OUTDIR/logs/metagenes.log"
    plotProfile -m "$M" -o "$MG/${mark}_${gs}.profile.png" \
      --plotTitle "${mark} @ TSS — ${gs}" 2>>"$OUTDIR/logs/metagenes.log"
    plotHeatmap -m "$M" -o "$MG/${mark}_${gs}.heatmap.png" \
      --colorMap RdBu_r --whatToShow 'heatmap and colorbar' 2>>"$OUTDIR/logs/metagenes.log"
  done
done
echo "metagenes done"
