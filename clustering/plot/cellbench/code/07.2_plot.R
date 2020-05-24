clumethod = as.character(commandArgs(trailingOnly = T)[1])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
pd <- readRDS(paste0('./clustering/plot/cellbench/plot/',clumethod,'/clustering_summary_pd.rds'))
pd1 = pd$pd1
pd2 = pd$pd2
pd3 = pd$pd3
pd4 = pd$pd4

colordf = readRDS('./resource/method_latent_color.rds')
source('./resource/function.R')
p1 <- ggplot() + geom_tile(data=pd1,aes(x=stat,y=method,fill=value)) +
  theme_hm(pd1$method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.8,0.9,0.95,1))  + xlab('') + ylab('') +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.6,title.position = 'top',title.hjust=0.3,title.vjust=NULL)) + theme(legend.position = 'bottom')+
  ggtitle(ifelse(clumethod=='kmeans','k-means clustering','Louvain clustering'))+
  theme(plot.title = element_text(hjust=1.2))

#########
mtd1 = 'no_imp'
f1 = 'sc_celseq2_5cl_p1.rds'
p2 <- ggplot(data=pd2, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .5) +
  theme_classic() + theme(legend.position = 'bottom', legend.title=element_blank()) +
  ggtitle(paste(mtd1,sub('.rds','',f1)))+ labs(fill='')+xlab('PC1') + ylab('PC2')+
  theme(axis.title.x = element_text(size=8),axis.title.y=element_text(size=8),plot.title = element_text(size=8), legend.text = element_text(size=8))+
  guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1)))
##########
mtd2 = 'MAGIC'
f2 = 'sc_celseq2_5cl_p1.rds'
p3 <- ggplot(data=pd3, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .5) +
  theme_classic() + theme(legend.position = 'bottom', legend.title=element_blank()) +
  ggtitle(paste(mtd2,sub('.rds','',f2))) +
  theme(axis.title.x = element_text(size=8),axis.title.y=element_text(size=8),plot.title = element_text(size=8), legend.text = element_text(size=8))+xlab('PC1') + ylab('PC2')+
  guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1)))

pd0 <- rbind(cbind(method = 'no_imp',pd2),cbind(method = 'MAGIC',pd3))
p0 <- ggplot(data=pd0, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .5) +
  theme_classic() + facet_wrap(~method, scales='free') + theme(legend.position = 'bottom', legend.title=element_blank()) +
  theme(axis.title.x = element_text(size=12),axis.title.y=element_text(size=12),plot.title = element_text(size=8), legend.text = element_text(size=7), axis.text=element_blank(),axis.ticks = element_blank())+xlab('Principal Component 1') + ylab('Principal Component 2')+
  guides(color = guide_legend(nrow = 1, byrow = TRUE,override.aes = list(size=2,alpha=1))) + 
  theme(strip.background = element_rect(color='white',fill='white'),strip.text.x = element_text(size=12))+
  theme(legend.spacing.x =unit(-0.1,'cm'))


library(ggrepel)
v <- colordf[match(pd4$mtd,colordf$fullName),'color']
names(v) <- pd4$mtd
p4 <- ggplot(pd4,aes(x=louvein,y=kmeans,label=mtd,color=mtd)) + geom_point() + geom_text_repel(size=3) + theme_bw() + theme(legend.position = 'none') +
  scale_color_manual(values=v)+ ggtitle('overall score') +
  xlab('Louvain')
  
pdf(paste0('./clustering/plot/cellbench/plot/',clumethod,'/clustering_summary.pdf'),width=6,height=5.8)
grid.arrange(p0,p1,p4,layout_matrix= cbind(matrix(1,nrow=7,ncol=3),rbind(matrix(0,ncol=4,nrow=3),matrix(4,ncol=4,nrow=4))))
dev.off()



