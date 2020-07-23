library(parallel)
library("scImpute")
packageVersion("scImpute")

scimpute( commandArgs(trailingOnly = T)[1],
         infile = "rds", 
         outfile = "rds", 
         out_dir = commandArgs(trailingOnly = T)[2],
         drop_thre = 0.5,
         Kcluster = 1,
         ncores = 24)

