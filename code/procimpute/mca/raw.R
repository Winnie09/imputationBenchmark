allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge')
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge/',f,'/genebycell.rds'))
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/raw/",f,'.rds'))
})

