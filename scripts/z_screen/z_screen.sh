#!/bin/bash
#SBATCH --job-name=zScreen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=qbiol
#SBATCH --mem=20GB                # memory (MB)
#SBATCH --time=4-4:00           # time (D-HH:MM)
#SBATCH -o zScreen.%N.%j.out     # STDOUT
#SBATCH -e zScreen.%N.%j.err     # STDERR
#
#SBATCH --array=1-32

echo "Start time: "; date

module load applications/R/4.0.3

Rscript --vanilla /home/uqjengel/yuhan/z_screen/z_screen.R $SLURM_ARRAY_TASK_ID

echo "End time: "; date