library(ggplot2)
library(gridExtra)
res1 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plotdata/diff.rds')
q = sub('.*_','',res1$setting)
res1 = res1[which(q==1 | q == 2 | q == 10), ]
res1[,1] <- factor(res1[,1],levels=names(sort(tapply(res1$diff,res1$method,median))))
p1 <- ggplot(res1,aes(x=method,y=diff,fill=method)) + geom_boxplot(outlier.shape = NA, outlier.colour = NA) + geom_point(size=.1) + 
      theme_classic() + coord_flip() + xlab('') + geom_hline(yintercept = 0,col='red') + 
      theme(legend.position = 'none') 

res2 <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plotdata/auc.rds')
q = sub('.*_','',res2$setting)
res2 = res2[which(q==1 | q == 2 | q == 10), ]
res2[,1] <- factor(res2[,1],levels=names(sort(tapply(res2$AUC,res2$method,median))))
p2 <- ggplot(res2,aes(x=method,y=AUC,fill=method)) + geom_boxplot(outlier.shape = NA, outlier.color = NA) + geom_point(size=.1) + 
      theme_classic() + coord_flip() + xlab('') + theme(legend.position = 'none')

pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/diff_roc_boxplot.pdf'),width=6,height=3.5)
grid.arrange(p1,p2,layout_matrix=matrix(1:2,nrow=1))
dev.off()
