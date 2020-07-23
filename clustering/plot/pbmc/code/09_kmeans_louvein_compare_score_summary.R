v1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_kmeans_pbmc.rds')
v2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_louvein_pbmc.rds')
int = intersect(names(v1[!is.na(v2)]), names(v2[!is.na(v2)]))
pd = data.frame(kmeans = v1[int], louvein = v2[int], mtd = int)
library(ggplot2)
library(reshape2)
library(ggrepel)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/compare/kmeans_louvein_score_summary.pdf',width=5,height=5)
ggplot(pd,aes(x=louvein,y=kmeans,label=mtd,color=mtd)) + geom_point() + geom_text_repel() +
  theme_classic() + theme(legend.position = 'none') 
dev.off()

