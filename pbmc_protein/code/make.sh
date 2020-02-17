#!/bin/bash -l
#SBATCH --partition=shared
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH -A hji7
ml R
Rscript make.R
