library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/scVI_latent'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/scVI_latent/',f,'.csv'),data.table = F)) 
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/scVI_latent/",f,'.rds'))
})

