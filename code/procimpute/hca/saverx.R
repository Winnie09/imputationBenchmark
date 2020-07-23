setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result')
allf <- list.files('./impute/hca/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/hca/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/hca/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/genebycell.rds')
  colnames(sexpr) = colnames(d)
  sexpr <- log2(sexpr + 1)
  dir.create("./procimpute/hca/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/hca/saverx/",f))
})

