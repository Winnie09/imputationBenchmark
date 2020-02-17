dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/data/processed/'
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
d_norm = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_genebycell.rds')
for (num.cell in c(1e3, 5e3)){
      set.seed(12345)
      id = sample(ncol(d), num.cell)
      m = d[,id]
      m_norm = d_norm[,id]
      dir.create(paste0(dir,num.cell),showWarnings = F)
      
      saveRDS(m, paste0(dir, num.cell,'/genebycell.rds'))      
      write.table(m, file=paste0(dir,num.cell,'/genebycell.csv'), sep=',', quote=F)
      tm = t(m)
      saveRDS(tm, paste0(dir,num.cell,'/cellbygene.rds'))
      write.table(tm, paste0(dir, num.cell,'/cellbygene.csv'), sep=',', quote=F)
      
      saveRDS(m_norm, paste0(dir, num.cell,'/norm_genebycell.rds'))      
      write.table(m_norm, paste0(dir,num.cell,'/norm_genebycell.csv'), sep=',', quote=F)
      tm_norm = t(m_norm)
      saveRDS(tm_norm, paste0(dir,num.cell,'/norm_cellbygene.rds'))
      write.table(tm_norm, paste0(dir, num.cell, '/norm_cellbygene.csv'), sep=',', quote=F)
}

