#!/bin/bash
# 01_qc.sh — FastQC on raw reads, fastp trim (Nextera auto-detect), MultiQC
set -euo pipefail
source "$(dirname "$0")/../config.sh"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample r1 r2 rest; do
  echo "== $sample =="
  # NOTE: MGB's conda mirror ships fastqc-rs (Rust FastQC rewrite) as `fastqc`;
  # it targets FastQC report parity (multiqc-safe) but uses --outdir, not -o.
  fastqc -t "$THREADS" --outdir "$OUTDIR/qc" "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
    2>>"$OUTDIR/logs/fastqc.log"
  fastp -i "$FASTQ_DIR/$r1" -I "$FASTQ_DIR/$r2" \
    -o "$OUTDIR/trimmed/${sample}_R1.fq.gz" -O "$OUTDIR/trimmed/${sample}_R2.fq.gz" \
    --detect_adapter_for_pe -w "$THREADS" \
    -h "$OUTDIR/qc/${sample}_fastp.html" -j "$OUTDIR/qc/${sample}_fastp.json" \
    2>>"$OUTDIR/logs/fastp.log"
done
multiqc -o "$OUTDIR/qc" "$OUTDIR/qc" 2>>"$OUTDIR/logs/multiqc.log"
echo "qc done"
