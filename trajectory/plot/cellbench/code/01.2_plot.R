trajmtd = commandArgs(trailingOnly = T)[1]
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
source('./resource/function.R')
colordf = readRDS('./resource/method_latent_color.rds')
ddir = paste0('./trajectory/result/cellbench/',trajmtd,'/cor_ov/')
pd <- readRDS(paste0('./trajectory/plot/cellbench/plot/',trajmtd,'/cor_ov_hm_pd.rds'))
pd1 = pd$pd1
pd2 = pd$pd2

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
p1 <- ggplot(pd1,aes(x=data,y=method,fill=cor)) + geom_tile() + theme_hm(pd1$method) + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.8,0.9,0.95,1),na.value = 'lightgrey')   + xlab('') + ylab('')+ labs(fill='correlation')+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.6,title.position = 'left',title.hjust=0.3,title.vjust=1)) + 
  theme(legend.position = 'bottom')

p2 <- ggplot(pd2,aes(x=data,y=method,fill=ov)) + geom_tile() + theme_hm(pd2$method) + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1),na.value = 'lightgrey') + xlab('') + ylab('') + labs(fill='overlap')+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.6,title.position = 'left',title.hjust=0.3,title.vjust=1)) + 
  theme(legend.position = 'bottom')

pdf(paste0('./trajectory/plot/cellbench/plot/',trajmtd,'/cor_ov_hm.pdf'), width=6,height=7)
grid.arrange(p1,p2,layout_matrix=matrix(1:2,1))
dev.off()


