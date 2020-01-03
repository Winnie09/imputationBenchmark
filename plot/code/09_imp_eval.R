allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/')
allmtd = setdiff(allmtd, c('viper','screcover'))

mtd = 'saver'
res = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/',mtd,'.rds'))

v = unlist(res)

sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
intergene = intersect(rownames(sexpr),rownames(bexpr))
colnames(bexpr)[which(colnames(bexpr) == 'H1-hESC')] = 'H1'
colnames(bexpr)[which(colnames(bexpr) == 'IMR-90')] = 'IMR90'
cl = 'H1'
library(ggplot2)
p1 <- ggplot(data=data.frame(sc = sexpr[intergene, which(sub('_.*','', colnames(sexpr)) == cl)[1]], bulk = bexpr[intergene,cl]), aes(x = sc, y = bulk)) +
  geom_point(color = 'black', size = 0.2, alpha = 0.5) + 
  theme_classic() +
  xlab('an imputed single cell') +
  ylab('bulk') +
  ggtitle(paste0('SAVER ', cl))

df <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/',mtd,'.rds'))
  res[[cl]]
})
library(reshape2)
pd = melt(df)
colnames(pd) = c('sc','method','cor')
pd$method = factor(as.character(pd$method), levels = names(sort(tapply(pd$cor, pd$method, median))))
p2 <- ggplot(data=pd, aes(x=factor(method), y = cor)) + geom_violin() + coord_flip() + theme_classic() + xlab('') + ylab('correlation between sc and bulk')


hmdf1 <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/',mtd,'.rds'))
  sapply(res, median)
})
hmdf2 <- sapply(allmtd, function(mtd){
  print(mtd)
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/10xcellline_cor/',mtd,'.rds'))){
    res = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/10xcellline_cor/',mtd,'.rds'))
    sapply(res, median)      
  } 
})
hmdf = cbind(t(hmdf1),t(hmdf2))
hmdf = melt(hmdf)
colnames(hmdf) <- c('method','ct','cor')
p3 <- ggplot(data=hmdf, aes(x=ct, y=method)) + geom_tile(aes(fill=cor)) + scale_fill_gradient(low = "black", high = "yellow") + xlab('') + ylab('') + theme(legend.position = 'bottom')
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/plot/plot/imp_eval.pdf',width=10,height=8)
grid.arrange(p1,p2,p3,layout_matrix=matrix(c(1,2,2,3,3,3),3))
dev.off()
