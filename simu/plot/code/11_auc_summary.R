library(reshape2)
library(ggplot2)
library(gridExtra)
DEMethod = commandArgs(trailingOnly = T)[1]
AUC = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/auc.rds'))
mtdorder = names(sort(tapply(AUC[,'AUC'], list(AUC[,'method']), mean, na.rm=T)))
AUC$method = factor(as.character(AUC$method), levels=mtdorder)
cellnum <- sapply(AUC[,'setting'],function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse = '_'))
propgene <- sapply(AUC[,'setting'],function(i) strsplit(i,'_')[[1]][3])
signal <- sapply(AUC[,'setting'],function(i) strsplit(i,'_')[[1]][4])

pd1 <- melt(tapply(AUC[,'AUC'],list(cellnum,AUC[,'method']),median))
pd2 <- melt(tapply(AUC[,'AUC'],list(propgene,AUC[,'method']),median))
pd3 <- melt(tapply(AUC[,'AUC'],list(signal,AUC[,'method']),median))
pd1[,1] <- as.character(pd1[,1])
pd1[,2] <- as.character(pd1[,2])
pd2[,1] <- as.character(pd2[,1])
pd2[,2] <- as.character(pd2[,2])
pd3[,1] <- as.character(pd3[,1])
pd3[,2] <- as.character(pd3[,2])

pd1[,1] <- paste0('cn:',pd1[,1])
pd2[,1] <- paste0('pg:',pd2[,1])
pd3[,1] <- paste0('si:',pd3[,1])

pd <- rbind(pd1,pd2,pd3)
colnames(pd) <- c('Type','Method','AUC')
pd$Method = factor(pd$Method, levels = mtdorder)
p1 <- ggplot(pd,aes(x=Method,y=Type,fill=AUC)) + geom_tile() + theme_classic() + coord_flip() + 
      theme(axis.text.x = element_text(angle=45,hjust=1))+xlab('')+ylab('')


plotdata2 = AUC[which(cellnum=='100_1000'),]
p2 <- ggplot(plotdata2,aes(x=method,y=AUC,fill=method)) + 
      # geom_boxplot(outlier.shape = NA, outlier.size = NA) + 
      geom_violin(scale='width')+
      geom_point(size=.1) + 
      theme_classic() + coord_flip() + theme(legend.position = 'none') + xlab('') + ylab('')+
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold'), axis.text.y = element_text(angle=45), axis.text.x = element_text(angle=45,hjust=1)) + 
      ggtitle('cn:100_1000') 

i = 'saver'
tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds'))
cellnum1 <- sapply(names(tmp),function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse = '_'))
tmp = tmp[which(cellnum1 == '100_1000')]
res <- lapply(tmp,function(j) {
      cbind(seq(0,1,0.01),approx(j[,2],j[,1],seq(0,1,0.01),rule=2)$y)
})
res <- do.call(rbind,res)
res <- tapply(res[,2],res[,1],quantile,c(0.25,0.5,0.75))
allres <- data.frame(i,as.numeric(names(res)),do.call(rbind,res),stringsAsFactors = F)
colnames(allres) <- c('method','x','y1','y2','y3')
p3 <- ggplot(data=allres) + geom_line(data=allres,aes(x=x,y=y2),col='#009999',size=2.5) + 
      geom_line(data=allres,aes(x=x,y=y1),col='grey') + geom_line(data=allres,aes(x=x,y=y3),col='grey') + 
      theme_classic() +geom_vline(xintercept = 0.25, colour='red') + xlab('real fdr') + ylab('sensitivity')+
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + ggtitle('SAVER 100_1000') 


df = as.data.frame(tmp[['100_1000_0.3_0.5']])
p4<- ggplot(data = df) + geom_line(data=df,aes(x=Real_FDR, y=Sensitivity),col='#009999',size=2) + theme_classic() +
      geom_vline(xintercept = 0.25, colour='red') + 
      xlab('real fdr') + ylab('sensitivity') + theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + 
      ggtitle('SAVER 100_1000_0.3_0.5') 
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/auc_summary.pdf'),width=10,height=7)
grid.arrange(p4,p3,p2,p1, layout_matrix=matrix(c(1,2,rep(3,2), rep(4,12)),nrow=4))
dev.off()
