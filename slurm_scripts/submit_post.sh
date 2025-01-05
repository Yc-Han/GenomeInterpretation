#!/bin/bash

###############################################################################
# Slurm job options (name, compute options, output, etc.)
###############################################################################
#SBATCH --job-name=post_cl_4c            # Job name
#SBATCH --partition=cpu             # Partition to use
#SBATCH --cpus-per-task=10            # Number of CPU cores per task
#SBATCH --ntasks=1                   # Total number of tasks
#SBATCH --mem-per-cpu=20G                   # Memory per node (adjust as needed)
#SBATCH --time=2-00:00:00            # Time limit: DD-HH:MM:SS
#SBATCH --output=post_cl_%j.out      # Standard output
#SBATCH --error=post_cl_%j.out       # Standard error

###############################################################################
# Load any required modules
###############################################################################
eval "$(/home/yhan/miniconda3/bin/conda shell.bash hook)"
conda activate tf
###############################################################################
# Execute your R script
###############################################################################
Rscript permutation/script/post_cl.R