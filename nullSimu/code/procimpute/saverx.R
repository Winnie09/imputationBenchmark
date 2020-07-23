setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/')
allf <- list.files('./impute/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  sexpr <- log2(sexpr + 1)
  dir.create("./procimpute/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/saverx/",f))
})

