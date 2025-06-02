#!/bin/bash
#SBATCH --job-name=te3                 # Job name
#SBATCH --mail-type=END,FAIL           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nicholas.davies@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                     # Run on a single core
#SBATCH --mem=2gb                      # Job memory request
#SBATCH --time=168:00:00               # Time limit hrs:min:sec
#SBATCH --output=te3.%A_%a.out
#SBATCH --array=1-50%40

./te ./Runs/3-Suppression/config.cfg $SLURM_ARRAY_TASK_ID
