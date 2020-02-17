#!/bin/bash -l
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=3
ml R
Rscript  makeCellBench.R 
