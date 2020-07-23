library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scscope'))
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge/',f,'/genebycell.rds'))
  sexpr <-  as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scscope/',f,'.csv'),data.table = F))
  sexpr <- log2(sexpr + 1)
  row.names(sexpr) <- row.names(d)
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/scscope/",f,'.rds'))
})

