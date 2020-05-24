library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/kmeans/medianSil/'))
res <- lapply(allf,function(f){
  df = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/kmeans/medianSil/',f))
})
d = do.call(rbind,res)
colnames(d)[6] = 'medianSil'
res1 <- sapply(3:6, function(i){
  tapply(d[,i],list(d[,1]),mean)  
})
colnames(res1) = colnames(d)[3:6]

allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/louvein/medianSil/'))
res <- lapply(allf,function(f){
  df = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/louvein/medianSil/',f))
  colnames(df)[which(colnames(df)=='meanSil')] <- 'medianSil'
  df
})
d = do.call(rbind,res)
colnames(d)[6] = 'medianSil'
res2 <- sapply(3:6, function(i){
  tapply(d[,i],list(d[,1]),mean)  
})
colnames(res2) = colnames(d)[3:6]

int <- intersect(rownames(res1),rownames(res2))
plist = list()
for (i in 1:ncol(res1)){
  pd <- data.frame(kmeans=res1[int,i],louvein=res2[int,i],mtd=int)
  plist[[i]] = ggplot(pd,aes(x=kmeans,y=louvein,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
    theme_minimal() + theme(legend.position = 'none') + ggtitle(colnames(res1)[i])
  
}
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/plot/compare/kmeans_louvein.pdf',width=12,height=12)
grid.arrange(grobs=plist, nrow=2)
dev.off()
