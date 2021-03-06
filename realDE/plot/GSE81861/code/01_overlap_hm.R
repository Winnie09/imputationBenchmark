library(reshape2)
library(ggplot2)
library(RColorBrewer)
ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder1 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean, na.rm=T)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean, na.rm=T)
saveRDS(mtdorder1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_mast_GSE81861.rds')
saveRDS(stat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/bulk_sc_diffgene_overlaps_mast_GSE81861.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder1)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/mast/',recursive = T,showWarnings = F)
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/mast/overlap_hm_pd.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/mast/overlap_hm.pdf',height=3.8,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1)) + labs(fill='overlap') + ggtitle('MAST') +
  theme(axis.text.y = element_text(size=10,color=ifelse(levels(pd$method)=='raw','red','black')),plot.title=element_text(size=10),legend.title=element_text(size=10))
dev.off()


ove = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder2 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean, na.rm=T)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean, na.rm=T)
saveRDS(mtdorder2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/bulk_sc_diffgene_overlaps_wilcox_GSE81861.rds')
saveRDS(stat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/bulk_sc_diffgene_overlaps_wilcox_GSE81861.rds')

pd$method = factor(as.character(pd$method), levels = mtdorder2)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/wilcox/',recursive = T)
saveRDS(pd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/wilcox/overlap_hm_pd.rds')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/wilcox/overlap_hm.pdf',height=3.8,width=4)
ggplot(pd,aes(x=data, y=method, fill=prop)) + geom_tile() + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1)) + labs(fill='overlap') + ggtitle('Wilcoxon') +
  theme(axis.text.y = element_text(size=10,color=ifelse(levels(pd$method)=='raw','red','black')),plot.title=element_text(size=10),legend.title=element_text(size=10))
dev.off()


