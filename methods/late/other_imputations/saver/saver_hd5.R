# Title     : run saver with 8 CPUs, input is 10x_genomics
# Objective : compare performance
# Created by: rui
# Created on: 4/3/18

library('SAVER')
library(cellrangerRkit)
library(doParallel)
getwd()
sessionInfo()

fname = '10kG_62.5kC.h5'
genome = 'mm10'
GeneBCMatrix = get_matrix_from_h5(fname ,genome)
df = exprs(GeneBCMatrix)

registerDoParallel(cores = 8)  # v0.4
saver4 <- saver(df)

write.csv(saver4$estimate, file='10kG_62.5kC.h5.saver.csv')