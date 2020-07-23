setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result')
allf <- list.files('./impute/pbmc/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/pbmc/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/pbmc/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/pbmc/sorted/genebycell.rds')
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d)
  sexpr <- log2(sexpr + 1)
  dir.create("./procimpute/pbmc/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/pbmc/saverx/",f))
})




