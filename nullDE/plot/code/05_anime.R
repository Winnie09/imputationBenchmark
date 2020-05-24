cnt = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/10_100_0_0/genebycell.rds'))
cnt = cnt[rowMeans(cnt)>1,1:50]
cnt = cnt[1:50,]

library(reshape2)
library(ggplot2)
pd = melt(cnt)
ggplot(pd,aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  theme_void() + xlab('') + ylab('')

set.seed(12345)
cnt = matrix(rnorm(900),nrow=30)
pd = melt(cnt)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/anime.pdf',width=5,height=5)
ggplot(pd,aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  theme_void() + xlab('') + ylab('')+
  theme(legend.position = 'none')
dev.off()
