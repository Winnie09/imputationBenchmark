library(ggplot2)
library(gridExtra)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/')
allmtd = setdiff(allmtd, 'zinbwave')
f = "cell1k_libsizeSmallVar.rds"
lib = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/', sub('.rds','',f), '/lib.rds'))      
for (i in 1:length(allmtd)){
  mtd = allmtd[i]
  print(f)
  lst = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd,'/',f))
  pc = lst$pc
  v = lst$var
  str(pc)
  if (i == 1) {
    pd = data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib)
  } else {
    pd = rbind(pd, data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib))
  }
}
pd$method = factor(pd$method, levels = c('raw', setdiff(unique(pd$method),'raw')))
colnames(pd)[4] = 'LibrarySizeFactor'
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca/fig1_1_pd_include_latent.rds')
