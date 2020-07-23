library(R.matlab)
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/autoimpute'))
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge/',sub('.csv','',f),'/genebycell.rds'))
  sexpr <- readMat(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/autoimpute/',f,'.mat'))$arr[1,,]
  row.names(sexpr) <- row.names(d)
saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/autoimpute/',f,'.rds'))
})



