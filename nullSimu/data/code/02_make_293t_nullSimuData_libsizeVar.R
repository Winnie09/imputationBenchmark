set.seed(12345)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')
## libsize with large variance
lib <- runif(1000,0.1,1.9)
saveRDS(lib,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeLargeVar/lib.rds')
simumat <- t(apply(d,1,function(i) {
      m <- mean(i)*lib
      rpois(1000,m)
}))
colnames(simumat) = paste0('cell',1:ncol(simumat))
simumat = simumat[rowMeans(simumat>0)>0.1,]
saveRDS(simumat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeLargeVar/genebycell.rds')
write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeLargeVar/genebycell.csv'), sep=",",quote=F)
saveRDS(t(simumat),'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeLargeVar/cellbygene.rds')
write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeLargeVar/cellbygene.csv'), sep=",",quote=F)

## libsize with small variance
set.seed(12345)
lib <- runif(1000,0.9,1.1)
saveRDS(lib,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/lib.rds')
simumat <- t(apply(d,1,function(i) {
      m <- mean(i)*lib
      rpois(1000,m)
}))
colnames(simumat) = paste0('cell',1:ncol(simumat))
simumat = simumat[rowMeans(simumat>0)>0.1,]
saveRDS(simumat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/genebycell.rds')
write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/genebycell.csv'), sep=",",quote=F)
saveRDS(t(simumat),'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/cellbygene.rds')
write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/cellbygene.csv'), sep=",",quote=F)
