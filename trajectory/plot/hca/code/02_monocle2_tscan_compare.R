setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
ddir1 = './trajectory/result/hca/monocle2/cor_ov/'
allf1 = list.files(ddir1)
colordf = readRDS('./resource/method_latent_color.rds')
res1 <- sapply(allf1,function(f){
  data = readRDS(paste0(ddir1,f))
  data = matrix(unlist(data),ncol=2,byrow=T)
  data[,1]/(data[,1]+data[,2])
})
names(res1) = sub('.rds','',names(res1))
names(res1) = colordf[match(names(res1),colordf$shortName),'fullName']

ddir2 = './trajectory/result/hca/tscan/cor_ov/'
allf2 = list.files(ddir2)
res2 <- sapply(allf2,function(f){
  data = readRDS(paste0(ddir2,f))
  data = matrix(unlist(data),ncol=2,byrow=T)
  data[,1]/(data[,1]+data[,2])
})
names(res2) = sub('.rds','',names(res2))
names(res2) = colordf[match(names(res2),colordf$shortName),'fullName']
library(reshape2)
library(ggplot2)
library(ggrepel)
library(gridExtra)
int <- intersect(names(res1),names(res2))
pd <- data.frame(monocle2=res1[int],tscan=res2[int],mtd=int)
pdf('./trajectory/plot/hca/plot/monocle2_tscan_compare.pdf',width=4,height=4)
ggplot(pd,aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none') + ggtitle('Percentage') +
  scale_color_manual(values=colordf[match(int,colordf$fullName),'color']) +
  ylab('TSCAN') + xlab('Monocle2') + xlim(c(0.55,1)) + ylim(c(0.6,1))
dev.off()
