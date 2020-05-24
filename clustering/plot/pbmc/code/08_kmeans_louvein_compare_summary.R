res1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_kmeans_pbmc.rds')
res2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_louvein_pbmc.rds')

library(ggplot2)
library(reshape2)
library(ggrepel)

int = intersect(names(res1), names(res2))
res1 = res1[int]
res2 = res2[int]
pd = data.frame(kmeans=res1, louvein=res2,mtd=int)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/compare/kmeans_louvein_summary.pdf', width=4,height=4)
ggplot(pd,aes(x=louvein,y=kmeans,label=mtd,color=mtd)) + geom_point() + geom_text_repel() +
  theme_bw() + theme(legend.position = 'none') 
dev.off()
