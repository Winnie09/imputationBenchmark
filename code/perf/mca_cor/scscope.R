library(data.table)
bexpr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/ct.rds')
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scscope'))
res <- lapply(allf,function(f) {
  print(f)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge/',f,'/genebycell.rds'))
  sexpr <-  as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scscope/',f,'.csv'),data.table = F))
  row.names(sexpr) <- row.names(d)
  be <- bexpr[,f]
  intgene <- intersect(row.names(sexpr),names(be))
  apply(sexpr[intgene,],2,cor,be[intgene],method='spearman')
})
names(res) <- allf
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/cor/scscope.rds')

