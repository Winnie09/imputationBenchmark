# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# ove_high = readRDS(paste0('./realDE/result/cellbench/overlap/mast/bulk_sc_diffgene_overlaps_moreExprGene.rds'))
# 
# v = apply(ove_high, 1, mean, na.rm=T)
# mtdorder1 = names(sort(rowMeans(ove_high,na.rm=T)))
# v    = v[mtdorder1]
# ove_low = readRDS(paste0('./realDE/result/cellbench/overlap/mast/bulk_sc_diffgene_overlaps_lessExprGene.rds'))
# u = apply(ove_low, 1, mean, na.rm=T)
# mtdorder2 = names(sort(rowMeans(ove_low,na.rm=T)))
# u = u[mtdorder2]


##################3
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
library(reshape2)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
coldf = readRDS('./resource/method_color.rds')
for (diffmtd in c('wilcox','mast')){
  rdir = paste0('./realDE/result/cellbench/overlap/',diffmtd,'/')
  ove_high =  readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
  ove_low = readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))
  rownames(ove_high) = coldf[match(rownames(ove_high),coldf[,'shortName']),'fullName']
  rownames(ove_low) = coldf[match(rownames(ove_low),coldf[,'shortName']),'fullName']
  ove_high = ove_high[complete.cases(ove_high),]
  ove_low = ove_low[complete.cases(ove_low),]
  mtdorder1 = names(sort(rowMeans(ove_high,na.rm=T)))
  mtdorder2 = names(sort(rowMeans(ove_low,na.rm=T)))
  pd_high = melt(ove_high)
  colnames(pd_high)=c('method','data','overlap')
  pd_high$method = factor(as.character(pd_high$method),levels=mtdorder1)
  
  napd = getNapd(coldf,pd_high,analysis.type='realDE',analysis.method=diffmtd,dataset='cellbench')
  if (!is.null(napd)) colnames(napd)[3] <- 'overlap'
  pd_high <- rbind(pd_high,napd[,1:3])
  pd_high$method = factor(as.character(pd_high$method),levels=c(levels(napd$method),mtdorder1))
  if (!is.null(napd)) levels(napd$method)  <- levels(pd_high$method)
  nav <- c('grey','black')
  names(nav) <- c('DifferentialFail','ImputationFail')
  
  if (is.null(napd)){
    p1 <- ggplot() + geom_tile(data=pd_high,aes(x=data,y=method,fill=overlap))+ 
      theme_minimal() + xlab('') + ylab('')+
      labs(fill='overlap') + ggtitle(paste0(ifelse(diffmtd=='mast','MAST','Wilcox'),',high-foldchange genes')) +
      theme_hm(pd_high$method)+
      scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
      guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))+
      scale_color_manual(values=nav)+
      theme(legend.key.width=unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'))
  } else {
    p1 <- ggplot() + geom_tile(data=pd_high,aes(x=data,y=method,fill=overlap))+ 
      geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
      theme_minimal() + xlab('') + ylab('')+
      labs(fill='overlap') + ggtitle(paste0(ifelse(diffmtd=='mast','MAST','Wilcox'),',high-foldchange genes')) +
      theme_hm(pd_high$method)+
      scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
      guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))+
      scale_color_manual(values=nav)+
      theme(legend.key.width=unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'))
  }
    
    
  
  pd_low = melt(ove_low)
  colnames(pd_low)=c('method','data','overlap')
  pd_low$method = factor(as.character(pd_low$method),levels=mtdorder2)
  
  napd = getNapd(coldf,pd_low,analysis.type='realDE',analysis.method=diffmtd,dataset='cellbench')
  if (!is.null(napd)) colnames(napd)[3] <- 'overlap'
  pd_low <- rbind(pd_low,napd[,1:3])
  pd_low$method = factor(as.character(pd_low$method),levels=c(levels(napd$method),mtdorder2))
  if (!is.null(napd)) levels(napd$method)  <- levels(pd_low$method)
  
  
  if (is.null(napd)){
    p2 <- ggplot() + geom_tile(data=pd_low,aes(x=data,y=method,fill=overlap))+ 
      theme_minimal() + xlab('') + ylab('')+
      labs(fill='overlap') + ggtitle(paste0(ifelse(diffmtd=='mast','MAST','Wilcox'),',low-foldchange genes')) +
      theme_hm(pd_low$method)+
      scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
      guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))+
      scale_color_manual(values=nav)+
      theme(legend.key.width=unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'))
    
  } else {
    p2 <- ggplot() + geom_tile(data=pd_low,aes(x=data,y=method,fill=overlap))+ 
      geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
      theme_minimal() + xlab('') + ylab('')+
      labs(fill='overlap') + ggtitle(paste0(ifelse(diffmtd=='mast','MAST','Wilcox'),',low-foldchange genes')) +
      theme_hm(pd_low$method)+
      scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
      guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))+
      scale_color_manual(values=nav)+
      theme(legend.key.width=unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'))
    
  }
    
  
  
  v_high = sort(rowMeans(ove_high,na.rm=T))
  v_low = sort(rowMeans(ove_low,na.rm=T)) 
  
  df_mean = data.frame(mean = v_high, method=names(v_high))
  pd_mean = melt(df_mean)
  colnames(pd_mean) = c('method','variable','mean')
  pd_mean$method = factor(as.character(pd_mean$method),levels=names(v_high))
  p3 = ggplot() + geom_tile(data=pd_mean,aes(x=variable,y=method,fill=mean))+
    scale_fill_gradientn(colors=brewer.pal(9,'BrBG')) +
    theme_classic()+
    guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))
  df_mean = data.frame(mean = v_low, method=names(v_low))
  pd_mean2 = melt(df_mean)
  colnames(pd_mean2) = c('method','variable','mean')
  pd_mean2$method = factor(as.character(pd_mean2$method),levels=names(v_low))
  p4 = ggplot() + geom_tile(data=pd_mean2,aes(x=variable,y=method,fill=mean))+
    scale_fill_gradientn(colors=brewer.pal(9,'BrBG')) +
    theme_classic()+
    guides(fill = guide_colourbar(barwidth = 0.6, barheight = 5,title.position = 'top',title.hjust=0.5,vjust=10))
  
  pdf(paste0('./realDE/plot/cellbench/plot/',diffmtd,'/overlap_hm_moreLessExprGene.pdf'),height=4.8,width=15)
  print(grid.arrange(p1,p2,p3,p4,nrow=1, layout_matrix = matrix(c(rep(1,2),rep(2,2),3,4),nrow=1)))
  dev.off()
}

