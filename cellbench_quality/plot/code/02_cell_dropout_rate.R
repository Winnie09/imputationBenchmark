allf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/'))
per = sapply(allf, function(f){
  d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/',f,'.rds'))
  perc_0 = sapply(1:ncol(d), function(i){
    sum(d[,i]==0)/nrow(d)
  })
})

data10x = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/raw/hg19.rds')
d = data10x
perc0_10x = sapply(1:ncol(d), function(i){
  sum(d[,i]==0)/nrow(d)
})

dataGSE = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds')
d = dataGSE
perc0_GSE = sapply(1:ncol(d), function(i){
  sum(d[,i]==0)/nrow(d)
})

per[['10x']] = perc0_10x
per[['fludigm']] = perc0_GSE

library(reshape2)
a = melt(per)
colnames(a) = c('perc','dataset')
a$dataset = factor(a$dataset,levels = c('10x','fludigm', setdiff(unique(a$dataset), c('10x','fludigm')))) 
library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/cellbench_quality/plot/plot/cell_dropout_rate.pdf',width=10,height=4)
ggplot(data=a,aes(x=dataset,y=perc)) + geom_violin() + geom_jitter(data=a, aes(x=dataset,y=perc,color=dataset),size=0.1,alpha=0.1,width=0.15)+
  theme_light() + xlab('') + ylab('cell dropout rate')+
  theme(axis.text.x = element_text(angle=45,hjust=1), legend.position = 'none')
dev.off()
