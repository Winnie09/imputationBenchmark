allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/alra'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- t(readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/alra/',f,'.rds')))
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/alra/",f,'.rds'))
})

