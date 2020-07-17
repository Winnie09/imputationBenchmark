af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench',pattern = 'celseq|mix')
res <- NULL
for (f in af) {
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/',f,'/genebycell.rds'))
  res <- rbind(res,data.frame(CDR=colMeans(d > 0),Sample=f,stringsAsFactors = F))
}
library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/cellbench_quality/plot/plot/cell_detection_rate_violin.pdf', width=10,height=4)
ggplot() + geom_violin(data=res,aes(x=Sample,y=CDR)) + geom_jitter(data=res,aes(x=Sample,y=CDR,color=Sample), alpha=.5, size=.1,width=.15)+
  theme_light() +
  theme(axis.text.x = element_text(angle=45,hjust=1), legend.position = 'none') + xlab('') + ylab('cell detection rate')
dev.off()
