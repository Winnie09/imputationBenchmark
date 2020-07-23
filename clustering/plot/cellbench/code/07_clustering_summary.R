#07_clustering_summary.R 
clumethod = as.character(commandArgs(trailingOnly = T)[1])
library(ggplot2)
library(gridExtra)
library(reshape2)
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_cellbench_statistics_scaled_summary.rds'))
v = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/clustering_',clumethod,'_cellbench.rds'))
pd0 = melt(d)
pd0$Var1 = factor(as.character(pd0$Var1), levels=v)

library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
p1 <- ggplot() + geom_tile(data=pd0,aes(x=Var2,y=Var1,fill=value)) + 
  theme_hm(pd0$Var1) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.8,0.9,0.95,1))  + 
  xlab('') + ylab('')
#########
mtd1 = 'raw'
af = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/',clumethod,'/',mtd1))
f1 = 'sc_celseq2_5cl_p1.rds'
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd1,'/',f1))
v = apply(d,1,sd)
d = d[v >= median(v), ]
d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
ct = factor(sub('.*:','', rownames(d)))
pd1 = data.frame(pc1 = d[,1], pc2 = d[,2], pc3 = d[,3], ct = ct)
p2 <- ggplot(data=pd1, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .5) + 
  theme_classic() + theme(legend.position = 'bottom', legend.title=element_blank()) + 
  ggtitle(paste(mtd1,sub('.rds','',f1)))+ labs(fill='')+xlab('PC1') + ylab('PC2')+
  theme(axis.title.x = element_text(size=8),axis.title.y=element_text(size=8),plot.title = element_text(size=8), legend.text = element_text(size=8))+
    guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1)))
##########
mtd2 = 'magic'
f2 = 'sc_celseq2_5cl_p1.rds'
af = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/',clumethod,'/',mtd2))
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd2,'/',f2))
v = apply(d,1,sd)
d = d[v >= median(v), ]
d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
ct = factor(sub('.*:','', rownames(d)))
pd2 = data.frame(pc1 = d[,1], pc2 = d[,2], pc3 = d[,3], ct = ct)
p3 <- ggplot(data=pd2, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .5) + 
  theme_classic() + theme(legend.position = 'bottom', legend.title=element_blank()) + 
  ggtitle(paste(mtd2,sub('.rds','',f2))) +
  theme(axis.title.x = element_text(size=8),axis.title.y=element_text(size=8),plot.title = element_text(size=8), legend.text = element_text(size=8))+xlab('PC1') + ylab('PC2')+
  guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1)))
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/plot/',clumethod,'/clustering_summary.pdf'),width=6,height=6)
grid.arrange(p1,p2,p3,layout_matrix=matrix(c(1,1:3),2))
dev.off()
