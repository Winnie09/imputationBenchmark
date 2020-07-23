library(reshape2)
library(ggplot2)
ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder1 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
saveRDS(mtdorder1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_mast.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder1)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/plot/mast/overlap_hm.pdf',height=6,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(angle=90)) + labs(fill='overlap')  + ggtitle('MAST')
dev.off()


ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder2 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
saveRDS(mtdorder2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_wilcox.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/plot/wilcox/overlap_hm.pdf',height=6,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(angle=90)) + labs(fill='overlap') + ggtitle('Wilcox') 
dev.off()
