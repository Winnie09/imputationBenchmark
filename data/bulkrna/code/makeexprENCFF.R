library(data.table)
allf <- list.files('/scratch/users/whou10@jhu.edu/Wenpin/rna_imputation/data/bulkrna/tsv')
expr <- sapply(allf,function(f) {
  tmp <- fread(paste0('/scratch/users/whou10@jhu.edu/Wenpin/rna_imputation/data/bulkrna/tsv/',f),data.table = F)
  tmp <- tmp[grep('^ENSMUSG',tmp[,1]),]
  v <- tmp$TPM
  names(v) <- tmp[,1]
  v
})
expr <- expr[rowSums(expr) > 0,]
expr <- log2(expr + 1)
load('/home-4/zji4@jhu.edu/scratch/resource/gn/res/grcm38.rda')
row.names(expr) <- sub('\\..*','',row.names(expr))
geneid <- sub('\\..*','',geneid)
expr <- expr[row.names(expr) %in% geneid,]
tmp <- genename[match(row.names(expr),geneid)]
rm <- rowMeans(expr)
kid <- sapply(unique(tmp),function(i) {
  id <- which(tmp==i)
  if(length(id) == 1) {
    id
  } else {
    id[which.max(rm[id])]
  }
})
expr <- expr[kid,]
row.names(expr) <- tmp[kid]
colnames(expr) <- sub('.tsv','',colnames(expr))
saveRDS(expr,file='/scratch/users/whou10@jhu.edu/Wenpin/rna_imputation/data/bulkrna/expr/ENCFF.rds')

