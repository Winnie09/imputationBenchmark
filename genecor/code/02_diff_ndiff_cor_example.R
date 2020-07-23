# /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/code/diff
library(parallel)
library(ggplot2)
library(reshape2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
allf <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m))
f = allf[1]
f = '100_1000_0.1_1.rds'
print(f)
dg = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',sub('.rds','',f),'/diffgn.rds'))
raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',sub('.rds','',f),'/genebycell.rds'))
#########
cormat <- sapply(allm, function(m){
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))) {
    d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))
    dg = dg[dg%in%row.names(d)]
    m1d = d[dg,]
    m1nd = d[setdiff(row.names(d),dg),]
    sample(as.vector(corfunc(t(m1d),t(m1nd))),1e5)
  }
})
cormat = cbind(cormat, raw = as.vector(corfunc(t(raw[dg,]), t(raw[setdiff(row.names(raw),dg),]))))

pd = melt(cormat)
pd = pd[,-1]
colnames(pd) = c('method','correlation')
pd[,'method'] = factor(pd[,'method'],levels=c('raw',setdiff(unique(pd[,'method']),'raw')))
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/genecor/plot/plot/',sub('.rds','',f),'_diff_ndiff_violin.pdf'),width=13,height=4)
ggplot(data=pd,aes(x=method,y=correlation,col=method)) + geom_violin() + theme_classic() + theme(legend.position = 'none')
dev.off()
##########
cormat <- sapply(allm, function(m){
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))) {
    d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))
    dg = dg[dg%in%row.names(d)]
    m1d = d[dg,]
    sample(as.vector(corfunc(t(m1d),t(m1d))),1e5)
  }
})
cormat = cbind(cormat, raw = as.vector(corfunc(t(raw[dg,]), t(raw[dg,]))))
pd = melt(cormat)
pd = pd[,-1]
colnames(pd) = c('method','correlation')
pd[,'method'] = factor(pd[,'method'],levels=c('raw',setdiff(unique(pd[,'method']),'raw')))
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/genecor/plot/plot/',sub('.rds','',f),'_diff_diff_violin.pdf'),width=13,height=4)
ggplot(data=pd,aes(x=method,y=correlation,col=method)) + geom_violin() + theme_classic() + theme(legend.position = 'none')
dev.off()

