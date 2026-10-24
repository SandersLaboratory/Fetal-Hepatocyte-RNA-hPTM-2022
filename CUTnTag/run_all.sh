#!/bin/bash
# run_all.sh — full pipeline driver
# Usage: bash run_all.sh [resume_after]
#   resume_after = name of the last successfully completed step; it and all
#   earlier steps are skipped. Step names:
#     00_setup 01_qc 01b_poolmap 02_align 03_peaks 04_deseq2 05_integration 06_metagenes 07_motifs
#   e.g. `bash run_all.sh 01_qc` resumes at 01b_poolmap (QC outputs are kept).
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
chmod +x scripts/*.sh

RESUME_AFTER="${1:-}"
skip=true
[ -z "$RESUME_AFTER" ] && skip=false

should_run() {
  # should_run <step> — returns 0 (run it) or 1 (skip it)
  local step="$1"
  if $skip; then
    echo ">>> skipping $step (resume after $RESUME_AFTER)"
    [ "$step" = "$RESUME_AFTER" ] && skip=false
    return 1
  fi
  return 0
}

should_run 00_setup     && bash scripts/00_setup.sh     2>&1 | tee "$OUTDIR/logs/00_setup.log"
should_run 01_qc        && bash scripts/01_qc.sh        2>&1 | tee "$OUTDIR/logs/01_qc.log"
should_run 01b_poolmap  && bash scripts/01b_poolmap.sh  2>&1 | tee "$OUTDIR/logs/01b_poolmap.log"
should_run 02_align     && bash scripts/02_align.sh     2>&1 | tee "$OUTDIR/logs/02_align.log"
should_run 03_peaks     && bash scripts/03_peaks.sh     2>&1 | tee "$OUTDIR/logs/03_peaks.log"
echo ">>> Inspect $OUTDIR/qc/metrics.tsv now: check dup rates, FRiP, spike-in."
echo ">>> Adjust MIN_FRAGS in scripts/04_consensus_deseq2.R if needed, then continue."
should_run 04_deseq2     && Rscript scripts/04_consensus_deseq2.R 2>&1 | tee "$OUTDIR/logs/04_deseq2.log"
should_run 05_integration && Rscript scripts/05_integration.R    2>&1 | tee "$OUTDIR/logs/05_integration.log"
should_run 06_metagenes  && bash scripts/06_metagenes.sh  2>&1 | tee "$OUTDIR/logs/06_metagenes.log"
should_run 07_motifs     && bash scripts/07_motifs.sh     2>&1 | tee "$OUTDIR/logs/07_motifs.log"
echo "ALL DONE — start at $OUTDIR/qc/multiqc_report.html and $OUTDIR/integration/"
