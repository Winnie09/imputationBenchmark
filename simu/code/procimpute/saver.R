allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saver'))
getf <- sub('.rds','',list.files("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saver/"))
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saver/',f,'.rds'))$estimate
  sexpr <- log2(sexpr + 1)
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saver/",f,'.rds'))
})
