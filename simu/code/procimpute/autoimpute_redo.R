library(R.matlab)
allf = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/code/procimpute/autoimpute_redof.txt')
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
  sexpr <- readMat(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/autoimpute/',f,'.mat'))$arr[1,,]
  row.names(sexpr) <- row.names(d)
  saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/autoimpute/',f,'.rds'))
  return(0)
})

