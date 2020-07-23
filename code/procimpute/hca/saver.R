allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/saver'))
res <- lapply(allf,function(f) {
print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/saver/',f,'.rds'))$estimate
  sexpr <- log2(sexpr + 1)
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/saver/",f,'.rds'))
})

