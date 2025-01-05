#!/bin/bash
#SBATCH --job-name=hmm_training
#SBATCH --output=hmm_training.out
#SBATCH --error=hmm_training.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8  # Adjust based on available cores
#SBATCH --time=12:00:00    # Adjust based on expected runtime
#SBATCH --mem=32G          # Adjust based on memory requirements
#SBATCH --partition=standard  # Adjust based on your cluster's partitions

# Load the R module (adjust based on your environment)
module load R

# Run the R script
Rscript permutation/script/hmm.R
