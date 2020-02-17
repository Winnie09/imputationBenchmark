set.seed(12345)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')## libsize##
gvar = apply(d,1,var)
id = names(gvar[gvar>median(gvar)])
set.seed(12345)
lib <- runif(1000,0.9,1.1)
saveRDS(lib,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_libsizeSmallVar/lib.rds')
simumat <- t(apply(d,1,function(i) {
  m <- mean(i)*lib
  rpois(1000,m)
}))
colnames(simumat) = paste0('cell',1:ncol(simumat))
simumat = simumat[id,]
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_smallLibVar_geneVarMedian/', showWarnings = F)
saveRDS(simumat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_smallLibVar_geneVarMedian/genebycell.rds')
write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_smallLibVar_geneVarMedian/genebycell.csv'), sep=",",quote=F)
saveRDS(t(simumat),'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_smallLibVar_geneVarMedian/cellbygene.rds')
write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k_smallLibVar_geneVarMedian/cellbygene.csv'), sep=",",quote=F)

