clumethod = as.character(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggrepel)
source('./resource/function.R')
# source('./resource/01_addFail.R')
coldf = readRDS('./resource/method_latent_color.rds')
mtdorder = readRDS(paste0('./result/perf/rank/clustering_',clumethod,'_pbmc.rds'))
mtdorder = coldf[match(mtdorder,coldf[,'shortName']),'fullName']
# mtdorder = c(setdiff(coldf$fullName, mtdorder),mtdorder)
d = readRDS(paste0('./clustering/result/stat/',clumethod,'_pbmc_statistics_scaled_summary.rds'))
# d = d[match(mtdorder,rownames(d)),]
pd  = melt(d)
colnames(pd) = c('method','data','eval')
pd = short2full(pd)
pd$method = factor(as.character(pd$method),levels=mtdorder)
dataset='pbmc'
p1 <- plotfunc(pd,napd = getNapd(coldf,pd,'clustering','kmeans',dataset=dataset),ifelse(clumethod=='kmeans','k-means clustering','Louvain clustering'),'value')+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5,title.position = 'top',title.hjust=0.3,title.vjust=NULL,label.theme = element_text(angle = 45,size=8))) + theme(legend.position = 'bottom')+ guides(color=guide_legend(title.position = 'top'))+
  theme(plot.title = element_text(hjust=1.2))
dev.off()
#########
mtd1 = 'raw'
d = readRDS(paste0('./clustering/plot/pbmc/pd/pc/',mtd1,'.rds'))
umap = readRDS(paste0('./clustering/plot/pbmc/pd/umap/',mtd1,'.rds'))
pd2 = data.frame(umap1 = umap[,1],umap2=umap[,2])

mtd2 = 'magic'
d = readRDS(paste0('./clustering/plot/pbmc/pd/pc/',mtd2,'.rds'))
umap = readRDS(paste0('./clustering/plot/pbmc/pd/umap/',mtd2,'.rds'))
pd3 = data.frame(umap1 = umap[,1],umap2=umap[,2])
ct = factor(sub(':.*','', rownames(d)),levels=c("cd34","cd4_t_helper","b_cells","cytotoxic_t","cd56_nk","naive_t","memory_t","regulatory_t","cd14_monocytes","naive_cytotoxic"))

pd23 <- rbind(cbind(method = 'no_imp',pd2, ct=ct),cbind(method = 'MAGIC',pd3,ct=ct))
getPalette = colorRampPalette(brewer.pal(10, "Paired"))
colorv = getPalette(length(unique(ct)))
p0 <- ggplot(data=pd23, aes(x = umap1, y = umap2, color = ct)) + geom_point(size = .2, alpha=0.1) +
  theme_classic() + theme(legend.position = 'bottom', legend.title=element_blank()) +
  theme(axis.title.x = element_text(size=12),axis.title.y=element_text(size=12),plot.title = element_text(size=8), legend.text = element_text(size=8), axis.text=element_blank(),axis.ticks = element_blank())+xlab('UMAP 1') + ylab('UMAP 2')+
  facet_wrap(~method, scales='free') + 
  theme(strip.background = element_rect(color='white',fill='white'),strip.text.x = element_text(size=12))+
  guides(color = guide_legend(nrow = 4, byrow = TRUE,override.aes = list(size=2,alpha=1))) +
  scale_color_manual(values=colorv) +
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))


res1 = readRDS('./result/perf/assess/clustering_kmeans_pbmc.rds')
res2 = readRDS('./result/perf/assess/clustering_louvein_pbmc.rds')
names(res1) = coldf[match(names(res1),coldf$shortName),'fullName']
names(res2) = coldf[match(names(res2),coldf$shortName),'fullName']
int = intersect(names(res1), names(res2))
res1 = res1[int]
res2 = res2[int]

pd4 = data.frame(kmeans=res1, louvein=res2,mtd=int)
p4 <- ggplot(pd4,aes(x=louvein,y=kmeans,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none') +xlim(c(0.35,0.85)) + ylim(c(0.3,0.75)) + 
  scale_color_manual(values=coldf[match(int,coldf$fullName),'color']) +
  ylab('k-means') + xlab('Louvain') + ggtitle('overall score')

library(gridExtra)
pdf(paste0('./clustering/plot/pbmc/plot/',clumethod,'/clustering_summary.pdf'),width=6.5,height=6.5)
grid.arrange(p0,p1,p4,layout_matrix= cbind(matrix(1,nrow=7,ncol=3),rbind(matrix(0,ncol=4,nrow=3),matrix(4,ncol=4,nrow=4))))
dev.off()



