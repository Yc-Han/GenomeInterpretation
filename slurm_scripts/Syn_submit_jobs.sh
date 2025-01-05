#!/bin/bash

while read FILE_NAME ; do
  while [ $(squeue -u $USER -h -t pending,running -r | wc -l) -ge 10 ] ; do
    sleep 10
  done
  sbatch --export=FILE_NAME="${FILE_NAME}" job_script.sh
done < synthetic/sourceseqs.csv
