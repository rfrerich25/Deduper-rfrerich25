#!/bin/bash

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --job-name=deduper                # optional: job name


/usr/bin/time -v ./Frerichs_deduper.py -f C1_SE_uniqAlign.sorted.sam \
 -o C1_SE_uniqAlign.sorted.deduped.sam -u STL96.txt