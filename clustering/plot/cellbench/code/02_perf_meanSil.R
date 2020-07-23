allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/meanSil/')

res <- lapply(allf,function(f){
  df = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/meanSil/',f))
})
d = do.call(rbind,res)

library(ggplot2)
library(gridExtra)
d <- d[d$method != 'mcimpute',]
p1 <- ggplot(d,aes(x=data,y=method,fill=Hacc)) + geom_tile() + theme_classic() + scale_fill_gradient(low='red',high='white') + theme(axis.text.x = element_text(angle=-45,hjust=-0.1)) + xlab('') + ylab('')
p2 <- ggplot(d,aes(x=data,y=method,fill=Hpur)) + geom_tile() + theme_classic() + scale_fill_gradient(low='red',high='white') + theme(axis.text.x = element_text(angle=-45,hjust=-0.1)) + xlab('') + ylab('')
p3 <- ggplot(d,aes(x=data,y=method,fill=ARI)) + geom_tile() + theme_classic() + scale_fill_gradient(high='red',low='white') + theme(axis.text.x = element_text(angle=-45,hjust=-0.1)) + xlab('') + ylab('')
p4 <- ggplot(d,aes(x=data,y=method,fill=meanSil)) + geom_tile() + theme_classic() + scale_fill_gradient(high='red',low='white') + theme(axis.text.x = element_text(angle=-45,hjust=-0.1)) + xlab('') + ylab('')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/plot/perf_meanSil.pdf', width=10,height=8)
grid.arrange(p1,p2,p3,p4,layout_matrix=matrix(1:4,2))
dev.off()
