library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/magic'))
getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/magic/'))
runf <- setdiff(allf, getf)
res <- lapply(runf, function(f){
  print(f)
  sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/magic/',f,'.csv'),data.table = F)))
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
  sexpr[sexpr<0] <- 0
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d)
  saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/magic/',f,'.rds'))
})

