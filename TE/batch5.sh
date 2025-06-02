#!/bin/bash
#SBATCH --job-name=te5                 # Job name
#SBATCH --mail-type=END,FAIL           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nicholas.davies@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                     # Run on a single core
#SBATCH --mem=2gb                      # Job memory request
#SBATCH --time=168:00:00               # Time limit hrs:min:sec
#SBATCH --output=te5.%A_%a.out
#SBATCH --array=109-216%40

./te ./Runs/5-Parasite/config.cfg $SLURM_ARRAY_TASK_ID
