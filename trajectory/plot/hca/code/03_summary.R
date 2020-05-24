trajmethod = commandArgs(trailingOnly = T)[1]
# trajmethod = 'monocle2'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
suppressMessages(library('monocle'))
library(TSCAN)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggrepel)

ct = readRDS('./data/processed/hca/ct/MantonBM6.rds')
ct[is.na(ct)] = 'unknown'
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colorvec = getPalette(length(unique(ct)))

mtd = 'knnsmooth'
print(mtd)
f = 'MantonBM6.rds'
if (file.exists(paste0('./trajectory/result/hca/',trajmethod,'/cds/',mtd,'/',f))){
  cds = readRDS(paste0('./trajectory/result/hca/',trajmethod,'/cds/',mtd,'/',f))
  cds$Celltype = ct[colnames(cds@reducedDimS)]
  if (trajmethod=='monocle2'){
    p1 <- plot_cell_trajectory(cds, cell_size = 0.5, color_by = 'Celltype') + ggtitle(mtd) + theme(legend.position = 'bottom')  + scale_color_manual(values=colorvec ) +
      guides(color = guide_legend(nrow = 3, byrow = TRUE,override.aes = list(size=2,alpha=1))) +
      theme(legend.spacing.x = unit(0, 'cm'),legend.spacing.y = unit(-0.5, 'cm'))+
      labs(color='') + ggtitle('kNN-smoothing') + xlab('Monocle2 Component 1') + ylab('Monocle2 Component 2')
  } else if (trajmethod=='tscan'){
    mcl <- cds[['mcl']]
    row.names(mcl[[1]]) <- paste0(sub(':','_',row.names(mcl[[1]])),1:nrow(mcl[[1]]))
    names(mcl[[3]]) <- paste0(sub(':','_',names(mcl[[3]])),1:nrow(mcl[[1]]))
    p1 <- plotmclust(mcl) + theme(legend.position = 'none') + ggtitle(mtd)
  }  
}

####### compare
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
saveRDS(res1,'./result/perf/assess/trajectory_monocle2_hca.rds')
saveRDS(rev(names(sort(res1))), './result/perf/rank/trajectory_monocle2_hca.rds')

ddir2 = './trajectory/result/hca/tscan/cor_ov/'
allf2 = list.files(ddir2)
res2 <- sapply(allf2,function(f){
  data = readRDS(paste0(ddir2,f))
  data = matrix(unlist(data),ncol=2,byrow=T)
  data[,1]/(data[,1]+data[,2])
})
names(res2) = sub('.rds','',names(res2))
names(res2) = colordf[match(names(res2),colordf$shortName),'fullName']
saveRDS(res2,'./result/perf/assess/trajectory_tscan_hca.rds')
saveRDS(rev(names(sort(res2))), './result/perf/rank/trajectory_tscan_hca.rds')

int <- intersect(names(res1),names(res2))
pd <- data.frame(monocle2=res1[int],tscan=res2[int],mtd=int)
v <- colordf[match(pd$mtd,colordf$fullName),'color']
names(v) <- pd$mtd

p2 <- ggplot(pd,aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(legend.position = 'none') + ggtitle('percentage') +
  scale_color_manual(values=v) +
  ylab('TSCAN') + xlab('Monocle2') + xlim(c(0.55,1)) + ylim(c(0.6,1))

pdf(paste0('./trajectory/plot/hca/plot/',trajmethod, '/summary.pdf'),width = 4, height = 8)
grid.arrange(p1,p2,layout_matrix = matrix(c(rep(1,12),rep(2,12)), ncol=3,byrow=T))
dev.off()  
