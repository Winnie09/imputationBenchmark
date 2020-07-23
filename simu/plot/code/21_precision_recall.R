DEMethod = commandArgs(trailingOnly = T)[1]
library(reshape2)
library(ggplot2)
allf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/')))
# allf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procdiff/'))
AUC = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/precision_recall_auc.rds'))
tapply(AUC[,'AUC'], list(AUC[,'method']), median)
apx <- sapply(allf, function(f){ ###############3
  print(f)
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',f,'.rds'))
  # res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procdiff/',f,'.rds'))
  tmp <- sapply(res, function(i) { ###########
    y = 1-i[,"Real_FDR"]
    x = i[,'Sensitivity']
    tmp = approx(x,y,xout=seq(0,1,length.out = 1000),rule=2)$y
  })
  rowMeans(tmp)
})

auc <- sapply(colnames(apx), function(i){
  a = apx[, i]
  b = seq(0,1,length.out = 1000)
  sum((a[-length(a)]+a[-1]) * diff(b)/2)
})
colnames(apx) = paste0(colnames(apx),'(',round(auc[match(colnames(apx), names(auc))],2),')')
#colnames(apx) = colnames(apx)[order(auc)]
pd = melt(apx)
pd$Var1 = rep(seq(0,1,length.out = 1000), length(unique(pd$Var2)))
colnames(pd) = c('x','method','y')
pd$method = factor(pd$method, levels = colnames(apx)[order(auc)])
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/precision_recall.pdf'),height=5,width=6.5)
ggplot(data=pd,aes(x=x, y=y, color=method)) + geom_line() +
  theme_classic() + theme(legend.position = 'top', legend.title=element_blank()) + 
  ylab('precision')+xlab('recall')
dev.off()
