library(data.table)
bexpr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/ct.rds')
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/dca'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/dca/',f,'/mean.tsv'),data.table=F)
  row.names(sexpr) <- sexpr[,1]
  sexpr <- as.matrix(sexpr[,-1])
  be <- bexpr[,f]
  intgene <- intersect(row.names(sexpr),names(be))
  apply(sexpr[intgene,],2,cor,be[intgene],method='spearman')
})
names(res) <- allf
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/cor/dca.rds')

