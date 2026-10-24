#!/bin/bash
# run_sbatch.sh — submit the full CUT&Tag pipeline as one Slurm job (Kraken)
# Usage: sbatch run_sbatch.sh [resume_after]
#   e.g. sbatch run_sbatch.sh 01_qc   # resume after QC (see run_all.sh)
# Monitor: squeue -u $USER ; tail -f ~/BrownRIH/cnt_output/logs/slurm_<jobid>.out
#SBATCH --job-name=cutntag
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 3-00:00:00
#SBATCH --output=~/BrownRIH/cnt_output/logs/slurm_%j.out
#SBATCH --error=~/BrownRIH/cnt_output/logs/slurm_%j.err

set -euo pipefail
# NOTE: adjust this path if your conda/mamba lives somewhere else
source /home/aohman/mambaforge/etc/profile.d/conda.sh
conda activate cutntag
cd ~/BrownRIH/cutntag-reanalysis
bash run_all.sh "$@"
