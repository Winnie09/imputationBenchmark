library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
evalfail <- c('SAUCIE')
source('./resource/function.R')
coldf = readRDS('./resource/method_color.rds')
## mast
pd = readRDS('./nullDE/plot/plot/sc_10x_5cl_mast_hm.rds') ############
pdlev <- coldf[match(levels(pd$method),coldf[,'shortName']),'fullName']
pd[,'method'] = coldf[match(pd[,'method'],coldf[,'shortName']),'fullName']

namet <- setdiff(coldf[,2],pd[,2])
namet <- c(setdiff(namet,evalfail), evalfail)
napd <- cbind(expand.grid(namet,unique(pd[,1])),NA)
colnames(napd) <- c('method','data','Num')

fullpd <- rbind(pd,napd)
fullpd$method <- factor(as.character(fullpd$method),levels=c(namet,pdlev))
fullpd = fullpd[,colnames(napd)]
napd$type <- 'ImputationFail'
napd$type[napd$method %in% evalfail] <- 'DifferentialFail'

pdf('./nullDE/plot/plot/sc_10x_5cl_mast_hm.pdf',height=5,width=5)  ### #####################
ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=Num)) + geom_tile(data=napd,aes(x=data,y=method,color=type),fill='white',size=0.2) + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('') + theme_hm(fullpd$method) + scale_color_manual(values=c('grey','black'))+
  guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top')) + 
  ggtitle('MAST')
dev.off()


## wilcox
pd = readRDS('./nullDE/plot/plot/sc_10x_5cl_wilcox_hm.rds') ############
pd[,'method'] = factor(coldf[match(pd[,'method'],coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd$method),coldf[,'shortName']),'fullName'])

pdf('./nullDE/plot/plot/sc_10x_5cl_wilcox_hm.pdf',height=4.2,width=4.2)  ############
ggplot(pd,aes(x=data, y=method, fill=Num)) + geom_tile() + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
  theme_minimal() + xlab('') + ylab('')+
  theme_hm(pd$method)+
  ggtitle('Wilcoxon')
dev.off()

