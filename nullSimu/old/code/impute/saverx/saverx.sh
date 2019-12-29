#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --mem=50G
#SBATCH --partition=shared
#SBATCH -A hji7

ml R
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/saverx.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/saverx/$1.rds
