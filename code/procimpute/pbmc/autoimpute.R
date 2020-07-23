library(R.matlab)
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/pbmc/autoimpute'))
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/pbmc/sorted/genebycell.rds')
  sexpr <- readMat(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/pbmc/autoimpute/',f,'.mat'))$arr[1,,]
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d)
 saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/autoimpute/',f,'.rds'))
})

