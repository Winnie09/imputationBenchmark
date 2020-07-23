d <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/transfer.csv',as.is=T,header=F)
d <- d[d[,2]!='',]
library(data.table)
tab <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/metadata.tsv',data.table=F)
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/ENCFF.rds')
mexpr <- sapply(d[,1],function(i) {
  print(i)
  tmp <- strsplit(d[d[,1]==i,2],';')[[1]]
  print(intersect(colnames(expr),tab[tab[,4] %in% tmp,1]))
  rowMeans(expr[,intersect(colnames(expr),tab[tab[,4] %in% tmp,1]),drop=F])
})
saveRDS(mexpr,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/ct.rds')
