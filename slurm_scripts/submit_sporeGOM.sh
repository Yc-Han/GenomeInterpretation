#!/bin/bash
#SBATCH --job-name=sporeGOM_job           # Job name
#SBATCH --partition=lrz-v100x2            # Partition with GPU (adjust if needed)
#SBATCH --nodes=1                         # Number of nodes (1 in this case)
#SBATCH --ntasks=1                        # Number of tasks (single R script execution)
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --gres=gpu:1                      # Request 1 GPU
#SBATCH --time=04:00:00                   # Time limit (adjust as needed)
#SBATCH --mem=32G                         # Memory per node (32GB)
#SBATCH --output=sporeGOM_%j.out          # Standard output (with job ID)
#SBATCH --error=sporeGOM_%j.err           # Standard error (with job ID)

# Load required modules
module load R/4.4.0  # Adjust R version if necessary
module load keras    # Load Keras if available as a module
module load cuda  # Load CUDA for GPU (adjust version as needed)

# Set up environment (if necessary, activate virtual env or similar)
# source activate your_environment  # Uncomment if using a conda environment

# Run the R script
Rscript sporeGOM.r
# End of script