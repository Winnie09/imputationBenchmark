
library(Matrix)
process10x_rmDupGenes <- function(genebycellmat){
  tb = genebycellmat
  tb <- tb[rowSums(tb) > 0,]
  gn = rownames(tb)  
  rs <- rowSums(tb)
  kid <- sapply(unique(gn),function(sid) {
    tmp <- which(gn==sid)
    if (length(tmp)==1) {
      tmp
    } else {
      tmp[which.max(rs[tmp])]
    }
  })
  tb <- tb[kid,]
  row.names(tb) <- gn[kid]
  tb <- tb[!grepl('^MT-',row.names(tb)),]
  tb = round(tb)
}

tb1 = as.matrix(readMM('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/matrix.mtx'))
g1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/genes.tsv')
barc1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/barcodes.tsv')
barc1 = paste0('293T_', barc1[,1])
rownames(tb1) = g1[,2]
colnames(tb1) = barc1
libsize1 = colSums(tb1)/1e3
mat = process10x_rmDupGenes(tb1)
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')
write.table(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.csv',sep=",",quote=F)
