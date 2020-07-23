clu = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/kmeans/magic/sorted.rds'))
ct = gsub(':.*','',names(clu))
df = data.frame(cluster = paste0('cluster',clu), celltype = ct, stringsAsFactors = F)

library(ggplot2)
tab <- table(df)
tab <- tab/rowSums(tab)
library(reshape2)
pd <- melt(tab)
pd$cluster <- factor(as.character(pd$cluster), levels = paste0('cluster', seq(1,10)))
# ggplot(pd,aes(x=cluster,y=celltype,fill=value)) + geom_tile() + theme_classic() + scale_fill_gradient2(low='grey',high='red', mid = 0.8)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/kmeans/Reviewer_magic_cluster_celltype_proportion.pdf', width = 8.5, height = 3.5)
ggplot(data = pd) +
  geom_bar(aes(x = cluster, y = value, fill = celltype), stat = 'identity', position = 'dodge') +
  theme_classic() +
  ylab('Celltype Proportion') +
  scale_fill_brewer(palette = 'Set3')
dev.off()
