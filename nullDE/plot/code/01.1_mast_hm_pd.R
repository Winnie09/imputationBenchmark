dataset = as.character(commandArgs(trailingOnly = T)[[1]])
## GSE81861, 10x, hca
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
## GSE81861, hca, 10x
allmtd = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/mast/'))

df <- sapply(allmtd, function(mtd){
  rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/mast/', mtd,'/res/')
  af = list.files(rdir)
  sapply(af, function(f){
    res <- readRDS(paste0(rdir,f))
    res$fdr = p.adjust(res$pvalue, method='fdr')
    sum(res$fdr<0.05)  
  })
})
if (is.list(df)){
  df = df[sapply(df,length)>0]
  pd = as.matrix(do.call(cbind, df))  
} else {
  pd = df
}

if (grepl('_0_0.rds.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds.rds','',rownames(pd))  
}  else if (grepl('_0_0.rds',rownames(pd)[1])){
  rownames(pd) = sub('_0_0.rds','',rownames(pd))  
} else {
  rownames(pd) = sub('.rds','',rownames(pd))
}

library(reshape2)
pd = melt(pd)
colnames(pd) = c('data','method','Num')
mtdorder = names(sort(tapply(pd[,'Num'], list(pd[,'method']), mean,na.rm=T), decreasing = T))
stat = tapply(pd[,'Num'], list(pd[,'method']), mean,na.rm=T)
saveRDS(mtdorder,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/nullDE_',dataset,'_mast.rds'))
saveRDS(stat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/nullDE_',dataset,'_mast.rds'))
pd$method = factor(as.character(pd$method), levels=mtdorder)
saveRDS(pd,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/',dataset,'_mast_hm.rds'))

