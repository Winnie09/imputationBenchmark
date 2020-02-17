library(reshape2)
library(ggplot2)
library(RColorBrewer)
ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder1 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_mast_cellbench.rds')
saveRDS(stat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/bulk_sc_diffgene_overlaps_mast_cellbench.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder1)
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast/overlap_hm_pd.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast/overlap_hm.pdf',height=3.8,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1)) + labs(fill='overlap')  + ggtitle('MAST')+
  theme(axis.text.y = element_text(size=10,color=ifelse(levels(pd$method)=='raw','red','black')),plot.title=element_text(size=10),legend.title=element_text(size=10))
dev.off()


ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder2 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
saveRDS(stat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')

pd$method = factor(as.character(pd$method), levels = mtdorder2)
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/wilcox/overlap_hm_pd.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/wilcox/overlap_hm.pdf',height=3.8,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  # scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1)) + labs(fill='overlap')  + ggtitle('Wilcoxon')+
  theme(axis.text.y = element_text(size=10,color=ifelse(levels(pd$method)=='raw','red','black')),plot.title=element_text(size=10),legend.title=element_text(size=10))
dev.off()


ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/foldchange/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder3 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder3,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
saveRDS(stat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder3)
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/foldchange/overlap_hm_pd.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/foldchange/overlap_hm.pdf',height=3.8,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  # scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1)) + labs(fill='overlap')  + ggtitle('Foldchange')+
  theme(axis.text.y = element_text(size=10,color=ifelse(levels(pd$method)=='raw','red','black')),plot.title=element_text(size=10),legend.title=element_text(size=10))
dev.off()
