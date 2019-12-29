#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --mem=110G
#SBATCH --partition=lrgmem
#SBATCH -A hji7

mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scVI_latent
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/scVI.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scVI/$1.csv /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/scVI_latent/$1.csv
