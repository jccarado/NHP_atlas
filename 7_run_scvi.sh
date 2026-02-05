#!/bin/bash
#SBATCH --job-name=process_samples
#SBATCH --output=sample_%a.out
#SBATCH --error=sample_%a.err
#SBATCH --array=1-3
#SBATCH --mem=300G
#SBATCH --cpus-per-task=12
#SBATCH --partition=gpu

SAMPLES_FILE="/peng_2/peng_lab/ngs_data/private/rhesus_atlas/cell_types.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES_FILE)
CELLTYPE=$(echo "$LINE" | awk -F'\t' '{print $1}')
python run_scvi.py "$CELLTYPE"

