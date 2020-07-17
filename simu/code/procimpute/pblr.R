allf = sub('.csv', '',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/pblr'))
getf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/pblr'))
runf <- setdiff(allf,getf)
res <- lapply(runf, function(f){
  print(f)
  sexpr <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/pblr/',f,'.csv'),as.is=T, header=F)
  sexpr = as.matrix(sexpr)
  sexpr = log2(sexpr + 1)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
  if (nrow(sexpr) != nrow(d)){
    print(f)
  } else {
    row.names(sexpr) <- row.names(d)
    colnames(sexpr) <- colnames(d)
    saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/pblr/',f,'.rds'))    
  }
})


