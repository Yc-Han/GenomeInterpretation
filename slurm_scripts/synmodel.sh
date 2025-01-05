#!/bin/bash
#SBATCH --job-name=syn_model            # Job name
#SBATCH --partition=gpu               # Partition to use
#SBATCH --gres=gpu:h100:1             # Number of CPU cores per task
#SBATCH --ntasks=1                    # Total number of tasks
#SBATCH --mem=100G                    # Memory per node (adjust as needed)
#SBATCH --time=2-00:00:00             # Time limit: DD-HH:MM:SS
#SBATCH --output=synmodel_%j.out       # Standard output
#SBATCH --error=synmodel_%j.out        # Standard error

eval "$(/home/yhan/miniconda3/bin/conda shell.bash hook)"

conda activate tf

Rscript synthetic/syn_model.R