#!/bin/bash
#SBATCH --job-name=popScreen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=qbiol
#SBATCH --mem=6GB                # memory (MB)
#SBATCH --time=0-24:00           # time (D-HH:MM)
#SBATCH -o popScreen.%N.%j.out     # STDOUT
#SBATCH -e popScreen.%N.%j.err     # STDERR
#
#SBATCH --array=1-32

echo "Start time: "; date

module load applications/R/4.0.3

Rscript --vanilla /home/uqjengel/yuhan/pop_screen/pop_screen.R $SLURM_ARRAY_TASK_ID

echo "End time: "; date