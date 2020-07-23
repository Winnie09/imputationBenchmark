# allf = sub('.rds','',list.files('/Volumes/My Passport/wenpin/rna_imputation/simu/result/procdiff/'))
DEMethod = commandArgs(trailingOnly = T)[1]
allf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/')))
s <- sapply(allf, function(f){
  print(f)
  # res = readRDS(paste0('/Volumes/My Passport/wenpin/rna_imputation/simu/result/procdiff/',f,'.rds'))
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',f,'.rds'))
  sapply(res, function(i) {
    tmp = i[,"Specificity"]
    tmp = tmp[!is.na(tmp)]
    median(tmp)
  })
})
s <- do.call(cbind, s)
s = s[,colnames(s)[order(apply(s,2,median),decreasing = T)]]
library(ggplot2)
library(reshape2)
pd = melt(s)
p <- factor(sapply(pd[,1], function(i) strsplit(as.character(i),'_')[[1]][3]))
q <- factor(sapply(pd[,1], function(i) strsplit(as.character(i),'_')[[1]][4]))
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/specificity_p.pdf'),height=3,width=8)
ggplot(data=pd, aes(x=Var2, y=value,fill=p)) + geom_boxplot(outlier.size = NA, outlier.shape = NA)+
  theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1), legend.position = 'top') + 
  ylab('specificity')+xlab('')
dev.off()
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/specificity_q.pdf'),height=3,width=11)
ggplot(data=pd, aes(x=Var2, y=value,fill=q)) + geom_boxplot(outlier.size = NA, outlier.shape = NA)+
  theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1), legend.position = 'top') + 
  ylab('specificity')+xlab('')
dev.off()


