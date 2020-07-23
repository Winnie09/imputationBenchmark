library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/saucie_latent'))
res <- lapply(allf, function(f){
  print(f)
  sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/saucie_latent/',f,'.csv'),data.table = F)))
  raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/',f,'/genebycell.rds'))
  colnames(sexpr) = colnames(raw)
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/saucie_latent/",f,'.rds'))  
})

