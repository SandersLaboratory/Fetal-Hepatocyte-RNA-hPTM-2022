#!/bin/bash
# 02_align.sh — bowtie2 to rn7, filter, markdup, dedup BAM; optional E. coli spike-in
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# chromosomes to keep (drop mitochondrial contigs)
CHROMS=$(cut -f1 "$GENOME_DIR/rn7.chrom.sizes" | grep -v -i -E '^(MT|chrM)' | tr '\n' ' ')

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample r1 r2 rest; do
  echo "== $sample =="
  R1="$OUTDIR/trimmed/${sample}_R1.fq.gz"; R2="$OUTDIR/trimmed/${sample}_R2.fq.gz"
  B="$OUTDIR/bam/$sample"

  # bowtie2 -> name-sort -> fixmate -m (adds MC/MQ mate tags) -> coordinate-sort.
  # markdup hard-fails without the fixmate step ("no MC tag" error), and
  # fixmate requires name-ordered input, hence the explicit name-sort.
  bowtie2 -p "$THREADS" --very-sensitive -X 2000 --no-mixed --no-discordant \
    -x "$GENOME_DIR/rn7" -1 "$R1" -2 "$R2" 2>"$B.bowtie2.log" \
    | samtools view -bS - \
    | samtools sort -n -@ "$THREADS" -o "$B.nsort.bam" -
  samtools fixmate -m -@ "$THREADS" "$B.nsort.bam" "$B.fixmate.bam"
  samtools sort -@ "$THREADS" -o "$B.sorted.bam" "$B.fixmate.bam"
  rm -f "$B.nsort.bam" "$B.fixmate.bam"
  samtools index "$B.sorted.bam"

  # mark duplicates (kept for QC), then produce deduped analysis BAM:
  # proper pairs only, MAPQ>=10, drop secondary/dup, drop chrM
  samtools markdup -@ "$THREADS" "$B.sorted.bam" "$B.markdup.bam" 2>"$B.markdup.log"
  samtools index "$B.markdup.bam"
  samtools view -b -@ "$THREADS" -f 2 -F 1280 -q 10 -o "$B.dedup.bam" "$B.markdup.bam" $CHROMS
  samtools index "$B.dedup.bam"

  samtools flagstat "$B.dedup.bam" > "$B.flagstat.txt"
  samtools stats "$B.dedup.bam" > "$B.stats.txt"

  # optional spike-in quantification (ChIP-Rx normalization)
  if [ "$USE_SPIKEIN" = true ]; then
    bowtie2 -p "$THREADS" --very-sensitive -x "$GENOME_DIR/ecoli" \
      -1 "$R1" -2 "$R2" 2>/dev/null \
      | samtools view -bS - \
      | samtools view -c -f 2 -F 256 - > "$B.spikein.count"
  fi
done
echo "alignment done"
