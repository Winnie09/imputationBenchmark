library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
ove = readRDS('./realDE/result/cellbench/overlap/mast/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder1 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder1,'./result/perf/rank/bulk_sc_diffgene_overlaps_mast_cellbench.rds')
saveRDS(stat,'./result/perf/assess/bulk_sc_diffgene_overlaps_mast_cellbench.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder1)
saveRDS(pd,'./realDE/plot/cellbench/plot/mast/overlap_hm_pd.rds')


ove = readRDS('./realDE/result/cellbench/overlap/wilcox/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder2 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder2,'./result/perf/rank/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
saveRDS(stat,'./result/perf/assess/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')

pd$method = factor(as.character(pd$method), levels = mtdorder2)
saveRDS(pd,'./realDE/plot/cellbench/plot/wilcox/overlap_hm_pd.rds')


ove = readRDS('./realDE/result/cellbench/overlap/foldchange/bulk_sc_diffgene_overlaps.rds')
ove = ove[rowMeans(is.na(ove))<1, ]
pd = melt(ove)
colnames(pd) = c('method','data','prop')
mtdorder3 = names(sort(tapply(pd[,'prop'], list(pd[,'method']), mean)))
stat = tapply(pd[,'prop'], list(pd[,'method']), mean)
saveRDS(mtdorder3,'./result/perf/rank/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
saveRDS(stat,'./result/perf/assess/bulk_sc_diffgene_overlaps_wilcox_cellbench.rds')
pd$method = factor(as.character(pd$method), levels = mtdorder3)
saveRDS(pd,'./realDE/plot/cellbench/plot/foldchange/overlap_hm_pd.rds')

