library(reshape2)
library(ggplot2)
library(ggrepel)
ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o1 <- rowMeans(ove)
ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
o2 <- rowMeans(ove)
int <- intersect(names(o1),names(o2))
pd <- data.frame(mast=o1[int],wilcoxon=o2[int],mtd=int)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/plot/wilcox_mast_compare.pdf',width=5,height=5)
ggplot(pd,aes(x=mast,y=wilcoxon,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + theme_minimal() + theme(legend.position = 'none')
dev.off()

