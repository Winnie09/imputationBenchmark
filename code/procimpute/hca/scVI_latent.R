library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/scVI_latent'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/scVI_latent/',f,'.csv'),data.table = F))
  sexpr = t(sexpr)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/',f,'/genebycell.rds')) 
  colnames(sexpr) = colnames(d)
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/scVI_latent/",f,'.rds'))
})

