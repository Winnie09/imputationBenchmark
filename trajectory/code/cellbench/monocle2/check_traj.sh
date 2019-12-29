#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=2:00:00
#SBATCH --mem=110G
ml R
Rscript  02_check_traj.R

