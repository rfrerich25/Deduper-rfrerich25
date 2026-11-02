#!/bin/bash

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --job-name=samtools                # optional: job name

/usr/bin/time -v pixi run samtools sort \
 -o C1_SE_uniqAlign.sorted.sam C1_SE_uniqAlign.sam