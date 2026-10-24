#!/usr/bin/env bash
# run_all_figures.sh — one-shot setup + run for the CUT&Tag figure scripts.
#
#   1. Restores the Ensembl release-110 mRatBN7.2 gene annotation (GTF)
#      into the genome folder. ONLY the annotation is needed for the
#      figures — no reference FASTA / Bowtie2 index (nothing is being
#      re-aligned). Skipped if the GTF is already there.
#   2. Makes sure the Python packages the figure scripts need are in
#      the cutntag conda environment (pyBigWig, pandas, matplotlib,
#      numpy), installing any that are missing. pandas/matplotlib/numpy
#      go in via mamba (prebuilt binaries) - pip source builds fail on
#      this cluster's system compiler, and current pip wheels target a
#      newer glibc than the cluster has.
#   3. Runs 01_regen_panels.py, then 02_fig6b_browser_tracks.py.
#      Figures land in $CNT_OUT/figures_final/ (PDF + PNG).
#
# Usage: copy this file next to 01_regen_panels.py and
# 02_fig6b_browser_tracks.py, then:   bash run_all_figures.sh
# Paths below default to the original cluster layout; override with
# environment variables if the folders live somewhere else:
#   CNT_OUT=/path/to/cnt_output GENOME_DIR=/path/to/genome bash run_all_figures.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

CNT_OUT="${CNT_OUT:-~/BrownRIH/cnt_output}"
GENOME_DIR="${GENOME_DIR:-~/BrownRIH/genome}"
GTF_NAME="Rattus_norvegicus.mRatBN7.2.110.gtf"
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/${GTF_NAME}.gz"

source /homes2/aohman/mambaforge/etc/profile.d/conda.sh
conda activate cutntag

echo "== [1/3] Gene annotation (GTF) =="
mkdir -p "$GENOME_DIR"
if [ -f "$GENOME_DIR/$GTF_NAME" ]; then
  echo "already present: $GENOME_DIR/$GTF_NAME"
else
  wget -q "$GTF_URL" -O "$GENOME_DIR/${GTF_NAME}.gz"
  gunzip "$GENOME_DIR/${GTF_NAME}.gz"
  echo "downloaded: $GENOME_DIR/$GTF_NAME"
fi
export CNT_GTF="$GENOME_DIR/$GTF_NAME"
export CNT_OUT

echo "== [2/3] Python packages =="
# pyBigWig for script 02; pandas/matplotlib/numpy for the plotting scripts.
python -c "import pyBigWig" 2>/dev/null || pip install pyBigWig
if ! python -c "import pandas, matplotlib, numpy" 2>/dev/null; then
  mamba install -y -c conda-forge pandas matplotlib numpy
fi
python -c "import pyBigWig, pandas, matplotlib, numpy; print('imports OK')"

echo "== [3/3] Figures =="
python "$SCRIPT_DIR/01_regen_panels.py"
python "$SCRIPT_DIR/02_fig6b_browser_tracks.py"
echo "DONE - figures in $CNT_OUT/figures_final/"
