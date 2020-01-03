setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd = readRDS('./plot/pd/15_imp_diff_eval_4subfig.pd.rds')
coldf = readRDS('./resource/method_color.rds')
pd1 = pd[['pd1']]
pd2 = pd[['pd2']]
pd3 = pd[['pd3']]
pd4 = pd[['pd4']]
pd3[,'Method'] = factor(coldf[match(pd3[,'Method'],coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd3$Method),coldf[,'shortName']),'fullName'])
pd4[,'Method'] = factor(coldf[match(pd4[,'Method'],coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd4$Method),coldf[,'shortName']),'fullName'])

library(reshape2)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
mtd = 'SAVER'
cl = '293T'
rawcor = cor(pd1[,'sc'], pd1[,'bulk'], method='spearman')
p1 = ggplot(data=pd1, aes(x = sc, y = bulk)) +
  geom_point(color = 'black', size = 0.2, alpha = 0.5) +
  theme_classic() + xlab('pseudobulk') + ylab('bulk (FPKM)') + 
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
  ggtitle(paste0('SCC=',round(rawcor,2)))
rawcor2 = cor(pd2[,'sc'],pd2[,'bulk'],method='spearman')
p2 <- ggplot(data=pd2, aes(x = sc, y = bulk)) +
  geom_point(color = 'black', size = 0.2, alpha = 0.5) + 
  theme_classic() +
  xlab('SAVER imputed single cells') + ylab('bulk (FPKM)') +
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
  ggtitle(paste0('SCC=',round(rawcor2,2)))
v <- coldf[match(pd3$Method,coldf$fullName),'color']
names(v) <- pd3$Method

p3 <- ggplot() + 
  geom_jitter(data=pd3,aes(x=Method,y=Correlation,color=Method),size=0.05,width=0.2,alpha=0.5)+
  geom_violin(data=pd3, aes(Method, Correlation),scale='width',fill=NA) + coord_flip() + theme_classic() + 
  geom_hline(yintercept = rawcor, color='red',linetype="dashed",size=1.5)+
  xlab('') + ylab('correlation between single cells and bulk') +
  theme(axis.text.y = element_text(size=12,color = ifelse(levels(pd3[,'Method'])=='no_imp','red','black')), legend.position = 'none',axis.text.x=element_text(size=12,color='black'), axis.title.x=element_text(size=12))+
  scale_color_manual(values=v)

p4 <- ggplot(data=pd4, aes(x=Dataset, y=Method)) + geom_tile(aes(fill=corMedian)) + 
  theme_hm(method_vec=pd4[,'Method'])+
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)) + xlab('') + ylab('') + 
  theme(legend.position = 'bottom')+
  labs(fill='Spearman correlation') +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.6,title.position = 'top',title.hjust=0.5,vjust=10))
pdf('./plot/plot/imp_diff_eval_4subfig.pdf',width=8,height=6)
grid.arrange(p1,p2,p3,p4,layout_matrix=matrix(c(1,3,3,3,2,3,3,3,rep(4,8)),4))
dev.off()



