trajmtd = commandArgs(trailingOnly = T)[1]
ddir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/cellbench/',trajmtd,'/cor_ov/')
pdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/cellbench/plot/',trajmtd,'/cor_ov_hm.pdf')
allf = list.files(ddir)
allf = setdiff(allf,'zinbwave.rds')
res <- lapply(allf,function(f){
  data = readRDS(paste0(ddir,f))
  data = matrix(unlist(data),ncol=2,byrow=T)
  data.frame(method = sub('.rds','',f), data = c("cellmix1.rds","cellmix2.rds","cellmix3.rds","cellmix4.rds","RNAmix_celseq2.rds","RNAmix_sortseq.rds"), 
             cor = data[,1], ov = data[,2])
  
})
df = do.call(rbind,res)
df[,'data'] = sub('.rds','',df[,'data'] )


mtdorder1 = names(sort(tapply(df[,'cor'], list(df[,'method']), mean, na.rm=T), decreasing = F))
stat = tapply(df[,'cor'], list(df[,'method']), mean, na.rm=T)
mtdorder1 = c(setdiff(unique(df$method), mtdorder1), mtdorder1)
saveRDS(mtdorder1,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/trajectory_cor_',trajmtd,'.rds'))
saveRDS(stat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/trajectory_cor_',trajmtd,'.rds'))
mtdorder2 = names(sort(tapply(df[,'ov'], list(df[,'method']), mean, na.rm=T), decreasing = F))
stat = tapply(df[,'ov'], list(df[,'method']), mean, na.rm=T)
mtdorder2 = c(setdiff(unique(df$method), mtdorder2), mtdorder2)
saveRDS(mtdorder2,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/trajectory_ov_',trajmtd,'.rds'))
saveRDS(stat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/trajectory_ov_',trajmtd,'.rds'))


#########
df = df[!df[,'method']%in%names(which(tapply(is.na(df[,'cor']),list(df[,'method']),sum)==6)),]

library(ggplot2)
library(gridExtra)
df1 = df
df1[,'method'] = factor(as.character(df1[,'method']), levels = intersect(mtdorder1,unique(df1[,'method'])))
p1 <- ggplot(df1,aes(x=data,y=method,fill=cor)) + geom_tile() + theme_classic() + scale_fill_gradient(low='black',high='yellow',limits = c(min(df[,3:4],na.rm=T),1),na.value = 'grey') + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('') + 
  theme(axis.text.y=element_text(color = ifelse(levels(df1[,'method'])=='raw','red','black')))
df2 = df
df2[,'method'] = factor(as.character(df2[,'method']), levels = intersect(mtdorder2,unique(df2[,'method'])))
p2 <- ggplot(df2,aes(x=data,y=method,fill=ov)) + geom_tile() + theme_classic() + scale_fill_gradient(low='black',high='yellow',limits=c(min(df[,3:4],na.rm=T),1),na.value = 'grey') +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + ylab('') +
  theme(axis.text.y=element_text(color = ifelse(levels(df2[,'method'])=='raw','red','black')))

pdf(pdir, width=10,height=3.5)
grid.arrange(p1,p2,layout_matrix=matrix(1:2,1))
dev.off()
