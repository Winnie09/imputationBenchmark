allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/cor')
res <- lapply(allf,function(f) {
  tmp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/cor/',f))
  sapply(tmp,median)
})
names(res) <- sub('.rds','',allf)
res <- res[sapply(res,length) > 20]
cn <- table(unlist(sapply(res,names)))
cn <- names(cn)[cn==length(res)]

res <- sapply(res,function(i) i[cn])
library(ggplot2)
library(reshape2)
mp <- apply(res,2,median,na.rm=T)
pd <- melt(res)
colnames(pd) <- c('cell','Method',"Correlation")
pd[,'Method'] <- factor(as.character(pd[,'Method']),levels = c('raw',setdiff(names(sort(mp,decreasing = T)),'raw')))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/plot/plot/correlation_mca.pdf',width=5,height=5)
ggplot(pd,aes(x=Method,y=Correlation,fill=Method)) + geom_boxplot() + geom_point() + theme_classic() + theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()
