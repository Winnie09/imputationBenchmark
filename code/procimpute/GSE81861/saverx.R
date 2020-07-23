setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result')
allf <- list.files('./impute/GSE81861/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/GSE81861/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/GSE81861/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  sexpr <- log2(sexpr + 1)
  dir.create("./procimpute/GSE81861/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/GSE81861/saverx/",f,'.rds'))
})

