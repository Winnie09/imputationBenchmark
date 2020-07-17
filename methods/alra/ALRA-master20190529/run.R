source('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/methods/alra/ALRA-master/alra.R')

A <- readRDS(paste0(commandArgs(trailingOnly = T)[1]))

A_norm <- normalize_data(A)

k_choice <- choose_k(A_norm)

A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]

saveRDS(A_norm_completed,paste0(commandArgs(trailingOnly = T)[2]))
