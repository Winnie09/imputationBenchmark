library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
source('./resource/function.R')
coldf = readRDS('./resource/method_color.rds')

plotfunc <- function(pd,napd=NULL,title='NULl',legend.lab = 'NULL') {
  fullpd <- rbind(pd,napd[,1:ncol(pd)])
  fullpd$method <- factor(as.character(fullpd$method),levels=c(levels(napd$method),levels(pd$method)))
  if (!is.null(napd)){
    # if (length(unique(napd$NA.reason)) > 1){
      nav <- c('black','grey')
      names(nav) <- c('ImputationFail',setdiff(unique(napd$NA.reason),'ImputationFail'))  
    # } else {
    #   nav = 'black'
    # }
    if (dataset!='hca'| length(levels(pd$data))<30){
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) + geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1),na.value = 'white') +
              theme_minimal() + xlab('') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=nav)+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              labs(fill=legend.lab)+
              ggtitle(title)+
              guides(color=guide_legend(order=1)))
    } else  {
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) + geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1),na.value = 'white') +
              theme_minimal() + xlab('466 pairs of cell types') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=nav)+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              labs(fill=legend.lab)+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              theme(axis.text.x=element_blank())+
              ggtitle(title))  
    }
    
  } else {
    if (dataset!='hca'){
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) +
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
              theme_minimal() + xlab('') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=c('grey','black'))+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              labs(fill=legend.lab)+
              ggtitle(title))  
    } else {
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) +
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1)) +
              theme_minimal() + xlab('466 pairs of cell types') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=c('grey','black'))+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              theme(axis.text.x=element_blank())+
              labs(fill=legend.lab)+
              ggtitle(title))
    }
    
  }
}

short2full <- function(pd){
  pdlev = levels(pd$method)
  pd[,'method'] = coldf[match(pd[,'method'],coldf[,'shortName']),'fullName']
  pd$method = factor(pd$method,levels=coldf[match(pdlev,coldf$shortName),'fullName'])
  pd
}

getNapd <- function(coldf,pd,analysis.type,analysis.method,dataset){
  namat <- setdiff(coldf[,2],pd[,1])
  namat = as.character(coldf[match(namat,coldf$fullName),'shortName'])
  if (length(namat) != 0){
    napd <- cbind(expand.grid(namat,unique(pd[,2])),NA)
    colnames(napd) <- c('method','data','eval')
    napd$NA.reason <- 'ImputationFail'
    for (i in unique(napd$method)){
      if (analysis.type == 'realDE'){
        reason.fail = 'DifferentialFail'
        
        result.path = paste0('./',analysis.type,'/result/',dataset,'/diff/',analysis.method,'/',i,'/res/')
        allf = sub('.rds','',list.files(paste0('./',analysis.type,'/result/',ifelse(dataset=='10x','cellbench',dataset),'/diff/',analysis.method,'/saver/res/')))
      } else if (analysis.type == 'nullDE') {
        reason.fail = 'DifferentialFail'
        result.path = paste0('./nullDE/result/',ifelse(dataset=='cellbench','10x',dataset),'/diff/',analysis.method,'/', i,'/res/')
        allf = sub('.rds','',list.files(paste0('./nullDE/result/',dataset,'/diff/',analysis.method,'/saver/res/')))
      } else if (analysis.type == 'clustering'){
        reason.fail = 'ClusteringFail'
      } else if (analysis.type == 'trajectory'){
        reason.fail  = 'TrajectoryFail'
      }  
      existf <- sub('.rds','',list.files(result.path))
      dataset = ifelse(dataset=='GSE81861','10xcellline',dataset)
      dataset = ifelse(dataset=='10x','10xcellline',dataset)
      if (dataset=='cellbench'){
        impute.exist <- intersect(list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/',dataset,'/',i)),'sc_10x_5cl.rds')
      } else {
        impute.exist <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/',dataset,'/',i))  
      }
      
      if (length(existf) > 0 ){
        napd[intersect(which(napd$method==i), match(setdiff(allf,existf),napd$data)),'NA.reason'] <- reason.fail
      } else if (length(impute.exist)>0){
        napd[which(napd$method==i),'NA.reason'] <- reason.fail
      }
    }  
    napd$method  = factor(coldf[match(napd$method, coldf$shortName),'fullName'],levels=coldf[match(levels(napd$method), coldf$shortName),'fullName'])
  } else {
    napd = NULL
  }
  
  napd
}

### realDE
## cellbench
for (dataset in c('cellbench','GSE81861','hca')){
  print(dataset)
  plotwidth = ifelse(dataset=='hca',6,4.8)
  pdf(paste0('./realDE/plot/',dataset,'/plot/mast/overlap_hm.pdf'),height=4.8,width=plotwidth)
  pd = short2full(readRDS(paste0('./realDE/plot/',dataset,'/plot/mast/overlap_hm_pd.rds')))
  colnames(pd)[3]='eval'
  plotfunc(pd,napd = getNapd(coldf,pd,'realDE','mast',dataset=dataset),'MAST','proportion')
  dev.off()
  
  pdf(paste0('./realDE/plot/',dataset,'/plot/wilcox/overlap_hm.pdf'),height=4.8,width=plotwidth)
  pd = short2full(readRDS(paste0('./realDE/plot/',dataset,'/plot/wilcox/overlap_hm_pd.rds')))
  colnames(pd)[3]='eval'
  plotfunc(pd,napd=getNapd(coldf,pd,'realDE','wilcox',dataset=dataset),'Wilcoxon','proportion')
  dev.off()
  
  # pdf(paste0('./realDE/plot/',dataset,'/plot/foldchange/overlap_hm.pdf'),height=5,width=4.8)
  # pd = short2full(readRDS(paste0('./realDE/plot/',dataset,'/plot/foldchange/overlap_hm_pd.rds')))
  # colnames(pd)[3]='eval'
  # plotfunc(pd,napd=getNapd(coldf,pd,'realDE','foldchange',dataset=dataset),'fold change','proportion')
  # dev.off()
}



### nullDE
for (dataset in c('10x','GSE81861','hca')){
  print(dataset)
  pdf(paste0('./nullDE/plot/plot/',dataset,'_wilcox_hm.pdf'),height=4.3,width=4.6)
  pd =short2full(readRDS(paste0('./nullDE/plot/plot/',dataset,'_wilcox_hm.rds')))
  colnames(pd)[3]='eval'
  pd = pd[,c('method','data','eval')]
  plotfunc(pd,napd=getNapd(coldf,pd,'nullDE','wilcox',dataset=dataset),'Wilcoxon','num.DEGs')
  dev.off()
  
  pdf(paste0('./nullDE/plot/plot/',dataset,'_mast_hm.pdf'),height=4.3,width=4.6)
  pd =short2full(readRDS(paste0('./nullDE/plot/plot/',dataset,'_mast_hm.rds')))
  colnames(pd)[3]='eval'
  pd = pd[,c('method','data','eval')]
  plotfunc(pd,napd=getNapd(coldf,pd,'nullDE','mast',dataset=dataset),'MAST','num.DEGs')
  dev.off()
}


