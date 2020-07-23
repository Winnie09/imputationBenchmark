library(R.matlab)
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/autoimpute'))
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/',f,'/genebycell.rds'))
  sexpr <- readMat(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/autoimpute/',f,'.mat'))$arr[1,,]
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d) 
 saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/autoimpute/',f,'.rds'))
})
  
