#!/bin/bash

#SBATCH --job-name=rds_conversion
#SBATCH --mem=250G
#SBATCH --output=slurm-logs/rds_int_%a.out
#SBATCH --error=slurm-logs/rds_int_%a.err
#SBATCH --array=1-3
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu

SAMPLES_FILE="/peng_2/peng_lab/ngs_data/private/rhesus_atlas/cell_types.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES_FILE)
CELLTYPE=$(echo "$LINE" | awk -F'\t' '{print $1}')
python anndata_to_rds_integration.py "$CELLTYPE"
