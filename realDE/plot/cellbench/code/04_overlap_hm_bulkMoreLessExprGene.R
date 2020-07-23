for (diffmtd in c('wilcox','mast')){
  rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/',diffmtd,'/')
  ove_high =  readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_bulkMoreExprGene.rds'))
  ove_low = readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_bulkLessExprGene.rds'))
  ove_high = ove_high[complete.cases(ove_high),]
  ove_low = ove_low[complete.cases(ove_low),]
  mtdorder1 = names(sort(rowMeans(ove_high,na.rm=T)))
  mtdorder2 = names(sort(rowMeans(ove_low,na.rm=T)))
  library(reshape2)
  library(gridExtra)
  library(ggplot2)
  pd_high = melt(ove_high)
  colnames(pd_high)=c('method','data','overlap')
  pd_high$method = factor(as.character(pd_high$method),levels=mtdorder1)
  
  p1 <- ggplot(data=pd_high,aes(x=data,y=method,fill=overlap)) + geom_tile()+ 
    scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
    theme_minimal() + xlab('') + ylab('')+
    theme(axis.text.x = element_text(angle=90)) + labs(fill='overlap') + ggtitle(paste0(diffmtd,',bulk high foldchange gene')) +
    theme(axis.text.y = element_text(color=ifelse(levels(pd_high$method)=='raw','red','black')))
  
  pd_low = melt(ove_low)
  colnames(pd_low)=c('method','data','overlap')
  pd_low$method = factor(as.character(pd_low$method),levels=mtdorder2)
  p2 <- ggplot(data=pd_low,aes(x=data,y=method,fill=overlap)) + geom_tile()+ 
    scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
    theme_minimal() + xlab('') + ylab('')+
    theme(axis.text.x = element_text(angle=90)) + labs(fill='overlap') + ggtitle(paste0(diffmtd,',bulk low foldchange gene')) +
    theme(axis.text.y = element_text(color=ifelse(levels(pd_low$method)=='raw','red','black')))
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/',diffmtd,'/overlap_hm_bulkMoreLessExprGene.pdf'),height=6,width=10)
  print(grid.arrange(p1,p2,nrow=1))
  dev.off()
}


