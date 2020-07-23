library(ggplot2)
library(gridExtra)
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
pd = readRDS('./nullSimu/plot/plot/pca/fig1_1_pd_include_latent.rds')

df = readRDS('./resource/method_latent_color.rds')
pd$method = df[match(as.character(pd$method), df[,'shortName']), 'fullName']
pdf(paste0('./nullSimu/plot/plot/pca/fig1_1_include_latent.pdf'),width=9,height=4)
ggplot(data = pd, aes(x=pc1,y=pc2,color=LibrarySizeFactor)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free',nrow=3) + 
  scale_color_gradient2(low="yellow",mid='green',high="blue",midpoint = mean(pd$LibrarySizeFactor)) + 
  theme_classic() +
  theme(legend.position = 'bottom', axis.ticks=element_blank(),axis.text=element_blank())+
  xlab('Principal Component 1') +
  ylab('Principal Component 2')+
  guides(color = guide_colourbar(barwidth = 8, barheight = 0.2))
dev.off() 

