library(data.table)
allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/drimpute'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/drimpute/',f,'.rds'))
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/drimpute/",f,'.rds'))
})

