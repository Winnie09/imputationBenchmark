setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
coldf = readRDS('./resource/method_color.rds')
res = readRDS('./pbmc_protein/res/perf.rds')
colnames(res) = coldf[match(colnames(res), coldf[,'shortName']),'fullName']
library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
pd = melt(res)
colnames(pd) = c('data','method','value')
pd  = pd[, c('method','data','value')]
mtdorder = names(sort(tapply(pd[,3], list(pd$method),mean,na.rm=T)))
pd$method = factor(as.character(pd$method), levels=c(setdiff(coldf$fullName,unique(pd$method)),mtdorder))

napd = cbind(expand.grid(setdiff(coldf$fullName,mtdorder),unique(pd$data)),NA)
colnames(napd) = c('method','data','value')
pd <- rbind(pd,napd)
#levels(napd$method) = levels(pd$method)
napd$NA.reason = 'ImputationFail'
# napd$method <- as.character(napd$method)
pdf('./pbmc_protein/res/plot/plot/roc_hm.pdf',width=4,height=4)
ggplot() + geom_tile(data=pd, aes(x=data,y=method,fill=value))+
  geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
  theme_minimal() + xlab('') + ylab('')+
  theme(axis.text.x = element_text(angle=90)) + labs(fill='ROC') + 
  theme_hm(pd$method) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.9,0.95,1)) +
  scale_color_manual(values='black')+
  theme(legend.key.width=unit(0.35,'cm'),legend.key.height = unit(0.35,'cm'))
dev.off()



