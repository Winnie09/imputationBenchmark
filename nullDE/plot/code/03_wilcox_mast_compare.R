setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
coldf = readRDS('./resource/method_color.rds')
dataset = as.character(commandArgs(trailingOnly = T)[1])
# dataset = '10x'
# dataset = 'hca'
# dataset = 'GSE81861'
# dataset = 'sc_10x_5cl'
library(reshape2)
library(ggplot2)
library(ggrepel)
o1 = readRDS(paste0('./result/perf/assess/nullDE_',dataset,'_mast.rds'))
o2 = readRDS(paste0('./result/perf/assess/nullDE_',dataset,'_wilcox.rds'))
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o1[int],Wilcoxon=o2[int],mtd=int)
pd$mtd = factor(coldf[match(pd$mtd, coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd$mtd), coldf[,'shortName']),'fullName'])
xmin <- ifelse(dataset=='10x',-10,ifelse(dataset=='GSE81861',-5,-10))
xmin <- ifelse(dataset=='sc_10x_5cl', -15, xmin)
x_add <- ifelse(dataset == 'hca', 4, ifelse(dataset=='sc_10x_5cl',8,0))

ymin <- ifelse(dataset=='10x',-12,ifelse(dataset=='GSE81861',-0.05,-0.1))
ymin <- ifelse(dataset=='sc_10x_5cl', -0.5, ymin)
v <- coldf[match(pd$mtd,coldf$fullName),'color']
names(v) <- pd$mtd

pdf(paste0('./nullDE/plot/plot/',dataset,'_wilcox_mast_compare.pdf'),width=4,height=4)
ggplot(pd,aes(x=MAST,y=Wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none')+
  scale_color_manual(values=v)+
  xlim(c(xmin,max(pd$MAST)+x_add)) + ylim(c(ymin,max(pd$Wilcoxon)+0.01))
dev.off()

rm(list=ls())

