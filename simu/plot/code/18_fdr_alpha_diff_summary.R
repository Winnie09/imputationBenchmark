DEMethod = commandArgs(trailingOnly = T)[1]
library(reshape2)
library(ggplot2)
library(gridExtra)
diff = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/', DEMethod,'/plotdata/fdr_alpha_diff.rds'))
diff = diff[!diff[,'method'] %in% 'zinbwave',]
diff$method = as.character(diff$method)
mtdorder = names(sort(tapply(diff[,'diff'], list(diff[,'method']), mean, na.rm=T)))
cellnum <- sapply(diff[,'setting'],function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse = '_'))
propgene <- sapply(diff[,'setting'],function(i) strsplit(i,'_')[[1]][3])
signal <- sapply(diff[,'setting'],function(i) strsplit(i,'_')[[1]][4])



pd1 <- melt(tapply(diff[,'diff'],list(cellnum,diff[,'method']),median))
pd2 <- melt(tapply(diff[,'diff'],list(propgene,diff[,'method']),median))
pd3 <- melt(tapply(diff[,'diff'],list(signal,diff[,'method']),median))
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
colnames(pd) <- c('Type','Method','Diff')
pd$Method = factor(pd$Method, levels=mtdorder)
p1 <- ggplot(pd,aes(x=Method,y=Type,fill=Diff)) + geom_tile() + theme_classic() + coord_flip() + 
      theme(axis.text.x = element_text(angle=45,hjust=1))+
      xlab('')+ylab('')+
      scale_fill_gradient2(high='red4',low='black')

plotdata2 = diff[which(signal=='1'),]
plotdata2$method = factor(plotdata2$method, levels=mtdorder)
p2 <- ggplot(plotdata2,aes(x=method,y=diff,fill=method)) + 
      geom_violin(trim=FALSE, scale='width') +
      geom_boxplot(outlier.shape = NA, outlier.size = NA, width=0.1) +
      geom_point(size=.1)  + coord_flip() + theme(legend.position = 'none') + xlab('') + ylab('')+ theme_minimal() +
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold'), axis.text.y = element_text(angle=45)) + 
      ggtitle('q=1') +
      geom_hline(color='red',yintercept = 0)#+ theme_classic()

i = 'saver'
tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds'))

df = as.data.frame(tmp[['100_1000_0.3_0.5']])
p3<- ggplot(data = df) + geom_line(data=df,aes(x=Reported_FDR, y=Real_FDR),col='#009999',size=2) + 
      theme_classic() +
      geom_abline(slope=1, intercept=0, colour='red') + 
      xlab(expression(alpha)) + ylab('fdr') + 
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + 
      ggtitle(paste0('SAVER 100_1000_0.3_0.5 (', round(diff[diff$method=='saver'&diff$setting=='100_1000_0.05_10','diff'],2),')')) 

pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/fdr_alpha_diff_summary.pdf'),width=7,height=6)
grid.arrange(p3,p2,p1, layout_matrix=matrix(c(1,1,2,2,2,rep(3,10)),nrow=5))
dev.off()
