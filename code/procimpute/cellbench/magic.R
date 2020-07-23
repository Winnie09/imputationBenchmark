library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/magic'))
res <- lapply(allf, function(f){
  sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/magic/',f,'.csv'),data.table = F)))
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/',f, '/genebycell.rds'))
  sexpr[sexpr<0] <- 0
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d)
  saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/magic/',f,'.rds'))
})
