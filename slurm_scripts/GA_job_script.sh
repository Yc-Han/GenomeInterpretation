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

# Get the file name from the list using the SLURM_ARRAY_TASK_ID
FILE_LIST="synthetic/sourceseqs.csv"
FILE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST | tr -d '"')

# Check if FILE_NAME is empty
if [ -z "$FILE_NAME" ]; then
  echo "No file name for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Run the R script with the file name as an argument
Rscript permutation/script/post_single.R "$FILE_NAME"
