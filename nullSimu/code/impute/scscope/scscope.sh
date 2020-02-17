#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --partition=lrgmem
#SBATCH --ntasks-per-node=2
#SBATCH -A hji7
mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scscope_latent/
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/scscope.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scscope/$1.csv /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scscope_latent/$1.csv
