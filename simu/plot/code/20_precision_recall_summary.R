DEMethod = commandArgs(trailingOnly = T)[1]
library(reshape2)
library(ggplot2)
library(gridExtra)
AUC = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/precision_recall_auc.rds'))
AUC = AUC[AUC$method != 'zinbwave',]
mtdorder = names(sort(tapply(AUC[,'AUC'], list(AUC[,'method']), mean, na.rm=T)))
AUC$method = factor(as.character(AUC$method), levels = mtdorder)
cellnum <- sapply(AUC[,'setting'],function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse = '_'))
propgene <- sapply(AUC[,'setting'],function(i) strsplit(i,'_')[[1]][3])
signal <- sapply(AUC[,'setting'],function(i) strsplit(i,'_')[[1]][4])

mtdorder = names(sort(tapply(AUC[,'AUC'], list(AUC[,'method']), median)))
pd1 <- melt(tapply(AUC[,'AUC'],list(cellnum,AUC[,'method']),median))
pd2 <- melt(tapply(AUC[,'AUC'],list(propgene,AUC[,'method']),median))
pd3 <- melt(tapply(AUC[,'AUC'],list(signal,AUC[,'method']),median))
pd1[,1] <- as.character(pd1[,1])
pd1[,2] <- as.character(pd1[,2])
pd2[,1] <- as.character(pd2[,1])
pd2[,2] <- as.character(pd2[,2])
pd3[,1] <- as.character(pd3[,1])
pd3[,2] <- as.character(pd3[,2])

pd1[pd1[,1]=='100_10',1] = '10_100'
pd1[pd1[,1]=='1000_10',1] = '10_1000'
pd1[pd1[,1]=='1000_100',1] = '100_1000'

pd1[,1] <- paste0('cn:',pd1[,1])
pd2[,1] <- paste0('p=',pd2[,1])
pd3[,1] <- paste0('q=',pd3[,1])

pd <- rbind(pd1,pd2,pd3)
colnames(pd) <- c('Type','Method','AUC')
pd$Method = factor(pd$Method, levels = mtdorder)
p1 <- ggplot(pd,aes(x=Method,y=Type,fill=AUC)) + geom_tile() +theme_minimal() + coord_flip() + 
      theme(axis.text.x = element_text(angle=45,hjust=1))+xlab('')+ylab('')


plotdata2 = AUC[which(cellnum=='100_1000'),]
p2 <- ggplot(plotdata2,aes(x=method,y=AUC,fill=method)) + 
      geom_violin(scale='width',trim=F)+
      geom_boxplot(outlier.shape = NA, outlier.size = NA,width=0.1) + 
      geom_point(size=.1) + theme_minimal() + coord_flip() + theme(legend.position = 'none') + 
      xlab('') + ylab('')+
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold'), axis.text.y = element_text(angle=45)) + 
      ggtitle('cn:100_1000') 

i = 'saver'
tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds'))
df = as.data.frame(tmp[['100_1000_0.3_0.05']])
df = cbind(df,Precision=1-df[,'Real_FDR'])
p3<- ggplot(data = df) + geom_line(data=df,aes(x=Sensitivity, y=Precision),col='#009999',size=2) + 
      theme_minimal() +
      xlab('recall') + ylab('precision') + 
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + 
      ggtitle('SAVER 100_1000_0.3_0.05') 

pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/precision_recall_summary.pdf'),width=7,height=6)
grid.arrange(p3,p2,p1, layout_matrix=matrix(c(1,1,2,2,2,rep(3,10)),nrow=5))
dev.off()
