clumethod = as.character(commandArgs(trailingOnly = T)[1])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
library(ggplot2)
library(gridExtra)
library(reshape2)
d = readRDS(paste0('./clustering/result/stat/',clumethod,'_cellbench_statistics_scaled_summary.rds'))
v = readRDS(paste0('./result/perf/rank/clustering_',clumethod,'_cellbench.rds'))
pd1 = melt(d)

colnames(pd1) <- c('method','stat','value')
colordf = readRDS('./resource/method_latent_color.rds')
v = colordf[match(v, colordf[,'shortName']),'fullName']
pd1$method = colordf[match(as.character(pd1$method), colordf[,'shortName']),'fullName']
pd1$method = factor(as.character(pd1$method),levels=v)

mtd1 = 'raw'
af = list.files(paste0('./clustering/result/cellbench/',clumethod,'/',mtd1))
f1 = 'sc_celseq2_5cl_p1.rds'
d = readRDS(paste0('./result/procimpute/cellbench/',mtd1,'/',f1))
v = apply(d,1,sd)
d = d[v >= median(v), ]
d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
ct = factor(sub('.*:','', rownames(d)))
pd2 = data.frame(pc1 = d[,1], pc2 = d[,2], pc3 = d[,3], ct = ct)


mtd2 = 'magic'
f2 = 'sc_celseq2_5cl_p1.rds'
af = list.files(paste0('./clustering/result/cellbench/',clumethod,'/',mtd2))
d = readRDS(paste0('./result/procimpute/cellbench/',mtd2,'/',f2))
v = apply(d,1,sd)
d = d[v >= median(v), ]
d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
ct = factor(sub('.*:','', rownames(d)))
pd3 = data.frame(pc1 = d[,1], pc2 = d[,2], pc3 = d[,3], ct = ct)

v1 = readRDS('./result/perf/assess/clustering_kmeans_cellbench.rds')
v2 = readRDS('./result/perf/assess/clustering_louvein_cellbench.rds')
int = intersect(names(v1), names(v2))
str(int)
pd4 = data.frame(kmeans = v1[int], louvein = v2[int], mtd = int,stringsAsFactors = F)
pd4$mtd = colordf[match(pd4$mtd,colordf$shortName),'fullName']
saveRDS(list(pd1=pd1,pd2=pd2,pd3=pd3,pd4=pd4), paste0('./clustering/plot/cellbench/plot/',clumethod,'/clustering_summary_pd.rds'))

