#!/bin/bash -l
#SBATCH --time=2:0:0
#SBATCH --mem=50G
#SBATCH --partition=shared
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/autoimpute.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/autoimpute/$1
