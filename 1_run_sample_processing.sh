#!/bin/bash
#SBATCH --job-name=process_samples
#SBATCH --output=sample_%a.out
#SBATCH --error=sample_%a.err
#SBATCH --array=1-58
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --partition=gpu

SAMPLES_FILE="/peng_2/peng_lab/ngs_data/private/rhesus_atlas/rhesus_atlas_metadata.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES_FILE)
BATCH_ID=$(echo "$LINE" | awk -F'\t' '{print $1}')
echo $BATCH_ID
DIRECTORY=$(echo "$LINE" | awk -F'\t' '{print $2}')
# Call your python script with these arguments
python /peng_2/peng_lab/scripts/rhesus_atlas/Atlas/run_preprocessing.py "$BATCH_ID" "$DIRECTORY"