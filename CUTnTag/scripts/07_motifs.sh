#!/bin/bash
# 07_motifs.sh — HOMER motif enrichment on fetal-retained promoter peaks
set -euo pipefail
source "$(dirname "$0")/../config.sh"

FA="$GENOME_DIR/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
for bed in "$OUTDIR"/integration/homer_beds/*.fetal_retained.bed; do
  [ -s "$bed" ] || { echo "skip empty $bed"; continue; }
  tag=$(basename "$bed" .fetal_retained.bed)
  out="$OUTDIR/motifs/$tag"; mkdir -p "$out"
  findMotifsGenome.pl "$bed" "$FA" "$out" -size 200 -mask -p "$THREADS" \
    2>>"$OUTDIR/logs/motifs.log"
done
echo "motifs done"
