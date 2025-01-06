#!/bin/bash
SEED=42
SEED_COUNTER=0
while read FILE_NAME ; do
  while [ $(squeue -u $USER -h -t pending,running -r | wc -l) -ge 10 ] ; do
    sleep 10
  done
  CURRENT_SEED=$((SEED + SEED_COUNTER))
  ((SEED_COUNTER++))
  sbatch --export=FILE_NAME="${FILE_NAME}",SEED="${CURRENT_SEED}" job_script.sh
done < synthetic/sourceseqs.csv
