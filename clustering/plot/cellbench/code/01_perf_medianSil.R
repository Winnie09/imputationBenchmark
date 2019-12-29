rm(list=ls())
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
setwd('/Users/wenpinhou/Dropbox/rna_imputation')
coldf = readRDS('./resource/method_latent_color.rds')
source('./resource/function.R')
# clumethod = commandArgs(trailingOnly = T)[1]
clumethod='louvein'
allf = list.files(paste0('./clustering/perf/cellbench/',clumethod,'/medianSil/'))
allf = setdiff(allf,'zinbwave.rds')
res <- lapply(allf,function(f){
      df = readRDS(paste0('./clustering/perf/cellbench/',clumethod,'/medianSil/',f))
      colnames(df)[6] = 'medianSil'
      df
})
d = do.call(rbind,res)
v = readRDS(paste0('./result/perf/rank/clustering_',clumethod,'_cellbench.rds'))
d$method = factor(coldf[match(as.character(d$method),coldf$shortName),'fullName'],levels=coldf[match(v,coldf$shortName),'fullName'])
library(ggplot2)
library(gridExtra)
library(viridis)

p1 <- ggplot(d,aes(x=data,y=method,fill=Hacc)) + geom_tile() + theme_hm2(d$method) + 
      theme(axis.text.x = element_text(angle=30,hjust=1)) + xlab('') + ylab('') +
   scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  +
   guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8,title.position = 'top',title.hjust=0.3,title.vjust=NULL)) + theme(legend.position = 'right')
   
   
p2 <- ggplot(d,aes(x=data,y=method,fill=Hpur)) + geom_tile(color='grey') + theme_hm2(d$method) + 
   theme(axis.text.x = element_text(angle=30,hjust=1)) + xlab('') + ylab('') +
   scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  +
   guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8,title.position = 'top',title.hjust=0.3,title.vjust=NULL)) + theme(legend.position = 'right')

p3 <- ggplot(d,aes(x=data,y=method,fill=ARI)) + geom_tile(color='grey') + theme_hm2(d$method) + 
   theme(axis.text.x = element_text(angle=30,hjust=1)) + xlab('') + ylab('') +
   scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  +
   guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8,title.position = 'top',title.hjust=0.3,title.vjust=NULL)) + theme(legend.position = 'right')

p4 <- ggplot(d,aes(x=data,y=method,fill=medianSil)) + geom_tile(color='grey') + theme_hm2(d$method) + 
   theme(axis.text.x = element_text(angle=30,hjust=1)) + xlab('') + ylab('') +
   scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')))  +
   guides(fill = guide_colourbar(barwidth = 0.6, barheight = 8,title.position = 'top',title.hjust=0.3,title.vjust=NULL)) + theme(legend.position = 'right')

pdf(paste0('./clustering/plot/cellbench/plot/',clumethod,'/perf_medianSil.pdf'), width=10,height=9)
grid.arrange(p1,p2,p3,p4,layout_matrix=matrix(1:4,2))
dev.off()

