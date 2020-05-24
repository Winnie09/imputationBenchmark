trajmtd = 'monocle2'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
source('./resource/function.R')
colordf = readRDS('./resource/method_latent_color.rds')
ddir = paste0('./trajectory/result/cellbench/',trajmtd,'/cor_ov/')
pd <- readRDS(paste0('./trajectory/plot/cellbench/plot/',trajmtd,'/cor_ov_hm_pd.rds'))
pd1 = pd$pd1
pd2 = pd$pd2

library(ggplot2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
napd = cbind(expand.grid(unique(pd1[is.na(pd1$cor),'method']), unique(pd1$data)), NA)
colnames(napd) <- c('method','data','eval')
napd$NA.reason <- 'TrajectoryFail'


p1 <- ggplot() + geom_tile(data=pd1,aes(x=data,y=method,fill=cor,color='grey')) + 
  # geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2)+
  theme_hm(pd1$method) + 
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.8,0.9,0.95,1),na.value = 'white')   + xlab('') + ylab('')+ labs(fill='correlation')+
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.6,title.position = 'top',title.hjust=0.3,title.vjust=1)) + 
  theme(legend.position = 'bottom') +
  ggtitle(ifelse(trajmtd=='tscan','TSCAN','Monocle2'))+
  scale_color_manual(values='grey')+
  theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))

p2 <- ggplot() + geom_tile(data=pd2,aes(x=data,y=method,fill=ov,color='grey')) + theme_hm(pd2$method) + 
  # geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2)+
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1),na.value = 'white') + xlab('') + ylab('') + labs(fill='overlap')+
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.6,title.position = 'top',title.hjust=0.3,title.vjust=1)) + 
  theme(legend.position = 'bottom')+
  ggtitle(ifelse(trajmtd=='tscan','TSCAN','Monocle2'))+
  scale_color_manual(values='grey')+
  theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))


cor1 = readRDS('./result/perf/assess/trajectory_cor_monocle2.rds')
cor2 = readRDS('./result/perf/assess/trajectory_cor_tscan.rds')
ov1 = readRDS('./result/perf/assess/trajectory_ov_monocle2.rds')
ov2 = readRDS('./result/perf/assess/trajectory_ov_tscan.rds')


short2full <- function(names_vec){
  colordf[match(names_vec,colordf$shortName),'fullName']
}
names(cor1) = short2full(names(cor1))
names(cor2) = short2full(names(cor2))
names(ov1) = short2full(names(ov1))
names(ov2) = short2full(names(ov2))


int <- intersect(names(cor1),names(cor2))
pd3 <- data.frame(monocle2=cor1[int],tscan=cor2[int],mtd=int)
v1 <- colordf[match(pd3$mtd,colordf$fullName),'color']
names(v1) <- pd3$mtd

p3 <- ggplot(pd3,aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_bw() + 
  theme(legend.position = 'none') + 
  ggtitle('correlation') +xlim(c(0.1,1.09)) +ylim(c(0.3,1))+
  scale_color_manual(values=v1)

int <- intersect(names(ov1), names(ov2))
pd4 <- data.frame(monocle2 = ov1[int], tscan = cor2[int], mtd=int)
v2 <- colordf[match(pd4$mtd,colordf$fullName),'color']
names(v2) <- pd4$mtd

p4 <- ggplot(pd4, aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none') + ggtitle('overlap') + xlim(c(0.6,1.09)) + ylim(c(0.3,1))+
  scale_color_manual(values=v2)

pdf(paste0('./trajectory/plot/cellbench/plot/',trajmtd,'/summary.pdf'), width=10,height=7)
grid.arrange(p1,p2,p3,p4,layout_matrix=matrix(c(1,1,2,2,3,4),2))
dev.off()


