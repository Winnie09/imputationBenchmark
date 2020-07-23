setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result')
allf <- list.files('./impute/10xcellline/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/10xcellline/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/10xcellline/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  sexpr <- log2(sexpr + 1)
  d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
  colnames(sexpr)  = colnames(d)
  dir.create("./procimpute/10xcellline/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/10xcellline/saverx/",f))
})

