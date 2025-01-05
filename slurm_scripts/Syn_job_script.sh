#!/bin/bash
#SBATCH --job-name=post_single           # Job name
#SBATCH --partition=cpu                 # Partition (adjust if needed)
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=1                   # Number of CPU cores per task
#SBATCH --time=48:00:00                     # Time limit
#SBATCH --mem=10G                           # Memory per node
#SBATCH --output=synthetic/post_single_%A_%a.out  # Standard output
#SBATCH --error=synthetic/post_single_%A_%a.out   # Standard error

# Load conda environment
eval "$(/home/yhan/miniconda3/bin/conda shell.bash hook)"
conda activate tf

# Run the R script with the file name as an argument
Rscript permutation/script/post_single.R "$FILE_NAME"
