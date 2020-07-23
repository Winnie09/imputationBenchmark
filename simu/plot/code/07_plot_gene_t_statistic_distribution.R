set.seed(12345)
allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/gene_t_statistic_distribution_plotdata')
mdata <- sdata <- tdata <- mdata2 <- sdata2 <- tdata2 <- NULL
for (f in allf){
  data = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/gene_t_statistic_distribution_plotdata/',f))
  id = sample(min(nrow(data$mdata), 100))
  mdata = rbind(mdata, data$mdata[id,])
  sdata = rbind(sdata, data$sdata[id,])
  tdata = rbind(tdata, data$tdata[id,])
  id = sample(min(nrow(data$mdata2), 500))
  mdata2 = rbind(mdata2, data$mdata2[id,])
  sdata2 = rbind(sdata2, data$sdata2[id,])
  tdata2 = rbind(tdata2, data$tdata2[id,])
}
mdata = cbind(mdata, type = 'diffgene')
sdata = cbind(sdata, type = 'diffgene')
tdata = cbind(tdata, type = 'diffgene')
mdata2 = cbind(mdata2, type = 'nondiffgene')
sdata2 = cbind(sdata2, type = 'nondiffgene')
tdata2 = cbind(tdata2, type = 'nondiffgene')

mdata = rbind(mdata, mdata2)
sdata = rbind(sdata, sdata2)
tdata = rbind(tdata,tdata2)
library(ggplot2)
library(gridExtra)
p1 <- ggplot(data=mdata) + geom_boxplot(aes(x=method, y = m, fill=type), outlier.size = 0.1,show.legend = T) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size=10,angle = 10), legend.position = 'bottom') +
  labs(title = '', x = '', y = 'absolute gene expr mean difference') +
  scale_y_continuous(trans='log10')

p2 <- ggplot(data=sdata) + geom_boxplot(aes(x=method, y = s, fill=type), outlier.size = 0.1,show.legend = T) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size=10,angle = 10), legend.position = 'bottom') +
  labs(title = '', x = '', y = ' standard deviation estimate') +
  scale_y_continuous(trans='log10')

p3 <- ggplot(data=tdata) + geom_boxplot(aes(x=method, y = t, fill=type), outlier.size = 0.1,show.legend = T) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size=10,angle = 10), legend.position = 'bottom') +
  labs(title = '', x = '', y = 'absolute t statistic') +
  scale_y_continuous(trans='log10')

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/gene_t_statistic_distribution.pdf', width=25,height=4)
grid.arrange(p1,p2,p3,layout_matrix=matrix(1:3,byrow = T,nrow=1))
dev.off()


