setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
coldf = readRDS('./resource/method_color.rds')
library(reshape2)
library(ggplot2)
library(ggrepel)
ove = readRDS('./realDE/result/cellbench/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
ove = readRDS('./realDE/result/cellbench/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- rowMeans(ove)
int <- intersect(names(o1),names(o2))
pd <- data.frame(MAST=o1[int],Wilcoxon=o2[int],mtd=int)
pd$mtd = factor(coldf[match(pd$mtd, coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd$mtd), coldf[,'shortName']),'fullName'])
pdf('./realDE/plot/cellbench/plot/wilcox_mast_compare.pdf',width=4,height=4)
v <- coldf[match(pd$mtd,coldf$fullName),'color']
names(v) <- pd$mtd
ggplot(pd,aes(x=MAST,y=Wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_bw() + theme(legend.position = 'none') + scale_color_manual(values=v)
dev.off()
