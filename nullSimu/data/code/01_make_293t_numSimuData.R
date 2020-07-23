d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')
d <- d[apply(d,1,var) > rowMeans(d),]
simumat <- t(apply(d,1,function(i) {
      m <- mean(i)
      v <- var(i)
      rnbinom(1000,size=m^2/(v-m),mu=m)  
}))
colnames(simumat) = paste0('cell',1:ncol(simumat))
simumat = simumat[rowMeans(simumat>0)>0.1,]
saveRDS(simumat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k/genebycell.rds')
write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k/genebycell.csv'), sep=",",quote=F)
saveRDS(t(simumat),'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k/cellbygene.rds')
write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/cell1k/cellbygene.csv'), sep=",",quote=F)
