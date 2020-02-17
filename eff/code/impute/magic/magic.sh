#!/bin/bash -l
#SBATCH --time=72:0:0
#SBATCH --partition=lrgmem
#SBATCH --cpus=24

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/magic.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/result/impute/magic/$1.csv
