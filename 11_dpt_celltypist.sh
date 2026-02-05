#!/bin/bash

#SBATCH --mem=300GB
#SBATCH --cpus-per-task=10
#SBATCH --output=dpt.out
#SBATCH --error=dpt.err
#SBATCH --partition=gpu

python 11_dpt_celltypist.py
