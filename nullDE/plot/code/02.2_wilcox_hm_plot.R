dataset = as.character(commandArgs(trailingOnly = T)[1])
pd <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/',dataset,'_wilcox_hm.rds'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
library(RColorBrewer)
library(ggplot2)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/',dataset,'_wilcox_hm.pdf'),height=3.8,width=4)
ggplot() + geom_tile(data=pd,aes(x=data,y=method,fill=Num)) + theme_hm(pd$method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.05,0.1,0.15,0.2,0.4,0.6,0.8,1))  + 
  xlab('') + ylab('')
dev.off()

