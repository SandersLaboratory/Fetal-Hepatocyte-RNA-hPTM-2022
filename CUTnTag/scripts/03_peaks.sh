#!/bin/bash
# 03_peaks.sh — SEACR peak calls (no-control), bigWigs, FRiP, per-sample QC table
set -euo pipefail
source "$(dirname "$0")/../config.sh"

METRICS="$OUTDIR/qc/metrics.tsv"
echo -e "sample\tmark\ttype\tpopulation\traw_reads\tmapped_pct\tdup_pct\tunique_frags\tspikein_frac\tn_peaks\tfrip" > "$METRICS"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample r1 r2 mark label type pop biorep pool; do
  echo "== $sample =="
  B="$OUTDIR/bam/$sample"
  DEDUP="$B.dedup.bam"

  # fragments -> bedgraph for SEACR (no IgG control available)
  # CRITICAL: bedtools bamtobed -bedpe REQUIRES queryname-sorted input — it
  # pairs consecutive reads blindly. The dedup BAM is coordinate-sorted, so
  # feeding it directly would silently pair unrelated reads, corrupting every
  # fragment (and therefore every SEACR peak, FRiP value, and all downstream
  # results). Name-sort to a temp file first, remove it after.
  samtools sort -n -@ "$THREADS" -o "$B.nsorted.bam" "$DEDUP"
  bedtools bamtobed -bedpe -i "$B.nsorted.bam" \
    | awk '{s=$2<$5?$2:$5; e=$3>$6?$3:$6; print $1"\t"s"\t"e}' \
    | sort -k1,1 -k2,2n \
    | bedtools genomecov -bg -i stdin -g "$GENOME_DIR/rn7.chrom.sizes" \
    > "$OUTDIR/peaks/${sample}.fragments.bedgraph"
  rm -f "$B.nsorted.bam"

  bash "$GENOME_DIR/SEACR_1.3.sh" \
    "$OUTDIR/peaks/${sample}.fragments.bedgraph" 0.01 non stringent \
    "$OUTDIR/peaks/${sample}.seacr" > "$OUTDIR/logs/seacr_${sample}.log" 2>&1 \
    || echo "SEACR failed for $sample — see $OUTDIR/logs/seacr_${sample}.log; continuing with 0 peaks"
  # SEACR writes ${prefix}.stringent.bed — weak samples may produce nothing
  PEAKS="$OUTDIR/peaks/${sample}.seacr.stringent.bed"
  if [ -f "$PEAKS" ]; then
    NPEAKS=$(wc -l < "$PEAKS")
    # FRiP
    INPEAKS=$(samtools view -c -L "$PEAKS" "$DEDUP")
  else
    NPEAKS=0; INPEAKS=0
  fi

  # bigWig (CPM; spike-in scale factor applied later in R if enabled)
  bamCoverage -b "$DEDUP" -o "$OUTDIR/bigwig/${sample}.bw" \
    --normalizeUsing CPM --binSize 10 -p "$THREADS" 2>>"$OUTDIR/logs/bw.log"

  TOTAL=$(samtools view -c "$DEDUP")
  FRIP=$(awk "BEGIN{print ($TOTAL>0)?$INPEAKS/$TOTAL:0}")

  # metrics from logs/stats — parsed defensively (|| true + NA/0 fallbacks):
  # a metrics miss must never kill the run under set -euo pipefail.
  RAW=$(( $(zcat "$FASTQ_DIR/$r1" | wc -l) / 4 ))
  # NOTE: bowtie2 logs "<pct>% overall alignment rate" — it never says
  # "mapped", so the old 'mapped...' regex matched nothing and the failed
  # grep (exit 1) aborted the whole script via set -e.
  MAPPED=$(grep -oP '\K[0-9.]+(?=% overall alignment rate)' "$B.bowtie2.log" | head -1 || true)
  if [ -z "$MAPPED" ]; then MAPPED="NA"; fi
  # dup rate from flagstat on the markdup BAM (version-independent; avoids
  # depending on samtools markdup's stderr prose format)
  samtools flagstat "$B.markdup.bam" > "$B.markdup.flagstat.txt" 2>/dev/null
  DUP=$(awk '/ in total/{t=$1} /duplicates/{d=$1} END{if(t>0) printf "%.2f", 100*d/t; else print "NA"}' "$B.markdup.flagstat.txt")
  SPIKE="NA"
  [ "$USE_SPIKEIN" = true ] && SPIKE=$(cat "$B.spikein.count")

  echo -e "${sample}\t${label}\t${type}\t${pop}\t${RAW}\t${MAPPED}\t${DUP}\t${TOTAL}\t${SPIKE}\t${NPEAKS}\t${FRIP}" >> "$METRICS"
done
echo "peaks done — see $METRICS"
