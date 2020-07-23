library(reshape2)
library(ggplot2)
library(gridExtra)
DEMethod = commandArgs(trailingOnly = T)[1]
diff = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod, '/plotdata/diff.rds'))
mtdorder = names(sort(tapply(diff[,'diff'], list(diff[,'method']), mean, na.rm=T)))
diff$method = factor(as.character(diff[,'method']), levels = mtdorder)
cellnum <- sapply(diff[,'setting'],function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse = '_'))
propgene <- sapply(diff[,'setting'],function(i) strsplit(i,'_')[[1]][3])
signal <- sapply(diff[,'setting'],function(i) strsplit(i,'_')[[1]][4])

pd1 <- melt(tapply(diff[,'diff'],list(cellnum,diff[,'method']),median, na.rm=T))
pd2 <- melt(tapply(diff[,'diff'],list(propgene,diff[,'method']),median, na.rm=T))
pd3 <- melt(tapply(diff[,'diff'],list(signal,diff[,'method']),median, na.rm=T))
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
colnames(pd) <- c('Type','Method','Diff')
pd$Method = factor(pd$Method, levels=mtdorder)
p1 <- ggplot(pd,aes(x=Method,y=Type,fill=Diff)) + geom_tile() + theme_classic() + coord_flip() + theme(axis.text.x = element_text(angle=45,hjust=1))+xlab('')+ylab('')


plotdata2 = diff[which(signal=='0.5'),]
p2 <- ggplot(plotdata2,aes(x=method,y=diff,fill=method)) + 
      # geom_boxplot(outlier.shape = NA, outlier.size = NA) + geom_point(size=.1) + 
      geom_violin(scale='width') +
      theme_classic() + coord_flip() + theme(legend.position = 'none') + xlab('') + ylab('')+
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + ggtitle('si:0.5') +
      geom_hline(color='red',yintercept = 0)

i = 'saver'
tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_', DEMethod,'/procdiff/',i,'.rds'))

df = as.data.frame(tmp[['100_1000_0.3_0.5']])
p3<- ggplot(data = df) + geom_line(data=df,aes(x=Real_FDR, y=Reported_FDR),col='#009999',size=2) + 
      theme_classic() +
      geom_abline(slope=1, intercept=0, colour='red') + 
      xlab('real fdr') + ylab('claimed fdr') + 
      theme(legend.position = 'none', plot.title = element_text(size = 8, face='bold')) + 
      ggtitle('SAVER 100_1000_0.3_0.5') 

dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/'), showWarnings = F)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/diff_summary.pdf'),width=8,height=7)
grid.arrange(p3,p2,p1, layout_matrix=matrix(c(1,2,2,2,rep(3,8)),nrow=4))
dev.off()
