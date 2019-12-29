#!/bin/bash -l
#SBATCH --time=24:0:0
#SBATCH --partition=lrgmem
#SBATCH --cpus=24

ml R
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/scscope.sh  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/data/processed/$1  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/result/impute/scscope/$1.csv
