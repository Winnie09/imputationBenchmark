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
colnames(pd) = c('Data','Method','Num')
mtdorder = names(sort(tapply(pd[,'Num'], list(pd[,'Method']), mean,na.rm=T), decreasing = T))
stat = tapply(pd[,'Num'], list(pd[,'Method']), mean,na.rm=T)
saveRDS(mtdorder,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/nullDE_',dataset,'_mast.rds'))
saveRDS(stat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/nullDE_',dataset,'_mast.rds'))
pd$Method = factor(as.character(pd$Method), levels=mtdorder)
library(ggplot2)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/',dataset,'_mast_hm.pdf'),height=3.8,width=4)
ggplot() + geom_tile(data=pd,aes(x=Data,y=Method,fill=Num)) + theme_hm(pd$Method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.05,0.1,0.15,0.2,0.4,0.6,0.8,1))  + 
  xlab('') + ylab('')
dev.off()
