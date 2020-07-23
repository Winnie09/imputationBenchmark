allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/pblr/'))
res <- lapply(allf, function(f){
  print(f)
  sexpr <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/pblr/',f,'.csv'),as.is=T,header=F)
  sexpr = as.matrix(sexpr)
  sexpr = log2(sexpr + 1)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge/',f,'/genebycell.rds'))
  row.names(sexpr) <- row.names(d)  ## lack of one row
  colnames(sexpr) <- colnames(d)
  saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/pblr/',f,'.rds'))   
})

