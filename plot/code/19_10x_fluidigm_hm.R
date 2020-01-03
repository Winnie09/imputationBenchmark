rm(list=ls())
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
colordf = readRDS('./resource/method_latent_color.rds')
source('./resource/function.R')
tb = read.csv('./doc/evalution_table.csv',stringsAsFactors = F)
tb = tb[!tb$evaluation%in%c('time','memory','scalability'),]
tb$category <- factor(tb$category,levels=c('imp_eval', 'differential', 'clustering', 'trajectory'))
tb$subcategory <- factor(tb$subcategory,levels=c('imp_eval', 'differential', 'null_differential','clustering', 'trajectory'))

mat = readRDS('./plot/plot/all_eval_score_noNA.rds')
rownames(mat) <- colordf[match(rownames(mat), colordf$shortName),'fullName']

eff = mat[,c('time','memory','scalability')]
mat = mat[,!colnames(mat)%in%c('time','memory','scalability')]

forestPlot <- function(mat) {
  library(ggplot2)
  m <- apply(mat,1,mean,na.rm=T)
  sd <- apply(mat,1,sd,na.rm=T)
  ld <- data.frame(Method=row.names(mat),low=m-2*sd,high=m+2*sd,stringsAsFactors = F)
  pd <- data.frame(Method=row.names(mat),p=m,stringsAsFactors = F)
  ld$Method <- factor(ld$Method,levels=pd[order(pd[,2]),1])
  pd$Method <- factor(pd$Method,levels=pd[order(pd[,2]),1])
  ggplot() + geom_linerange(data=ld,aes(x=Method,ymin=low,ymax=high)) + geom_point(data=pd,aes(x=Method,y=p)) + coord_flip()
}

str(mat)
library(ggplot2)
library(reshape2)
library(viridis)
UMI <- tb[tb$platform!='Fluidigm','evaluation']
Fluidigm <- tb[tb$platform!='UMI','evaluation']
Fluidigm_cellline <- tb[tb$platform!='UMI'& tb$type=='cellline' ,'evaluation']
UMI_celline <- tb[tb$platform!='Fuidigm'& tb$type=='cellline' ,'evaluation']
UMI_tissue <- tb[tb$platform!='Fuidigm'& tb$type=='tissue' ,'evaluation']

for (ii in 1:6){
    select = list(UMI=UMI,Fluidigm=Fluidigm, Fluidigm_cellline=Fluidigm_cellline,UMI_celline=UMI_celline,UMI_tissue=UMI_tissue, all=colnames(mat))[ii]
    df = mat[,select[[1]]]
    ## forest plot
    pdf(paste0('./plot/plot/all_eval_score_detail_',names(select),'_forest.pdf'),width=4,height=6)
    print(forestPlot(df) + ggtitle(names(select)))
    dev.off()
    ### hm data
    pd = melt(df)
    colnames(pd) = c('method','evaluation','score') 
    pd = cbind(pd, dataset=tb[match(pd$evaluation,tb$evaluation),'dataset'], subcategory=tb[match(pd$evaluation,tb$evaluation),'subcategory'])
    ## callapse dataset
    res1 = tapply(pd$score, list(pd$method,pd$subcategory,pd$dataset), mean, na.rm=T)
    tmpmat = matrix(NA,nrow=dim(res1)[1],ncol=dim(res1)[2])
    for (i in 1:nrow(tmpmat)){
      for (j in 1:ncol(tmpmat)){
        tmpmat[i,j] <- mean(res1[i,j,],na.rm=T)
      }
    }
    dimnames(tmpmat) <- dimnames(res1)[1:2]
    if ('null_differential'%in%colnames(tmpmat) & 'differential'%in%colnames(tmpmat)){
      tmpmat[, 'differential'] = rowMeans(tmpmat[,c('differential', 'null_differential')])
      tmpmat = tmpmat[,!colnames(tmpmat)%in%'null_differential']  
    }# else if (!'differential'%in%colnames(tmpmat)){
    #   colnames(tmpmat)[colnames(tmpmat)=='null_differential'] <- 'differential'
    # }
    
    ## callapse category
    
    pd2 <- melt(tmpmat)
    colnames(pd2) = c('method','evaluation','score') 
    mtdorder <- names(sort(tapply(pd2$score,list(pd2$method),mean,na.rm=T)))
    pd2$method = factor(as.character(pd2$method), levels=mtdorder)
    pd$method = factor(as.character(pd$method),levels=mtdorder)
    ### hm in details
    if (names(select)=='all'){
      pdf(paste0('./plot/plot/all_eval_score_detail_',names(select),'.pdf'),width=7,height=7)
    } else {
      pdf(paste0('./plot/plot/all_eval_score_detail_',names(select),'.pdf'),width=5,height=7)
    }
    print(ggplot(data=pd, aes(x=evaluation,y=method,fill=score)) + geom_tile()+
      scale_fill_viridis(discrete=FALSE, na.value='white') +
      theme_hm(pd$method) + xlab('') + ylab('')+
      labs(fill='score') + 
      ggtitle(names(select)))
    dev.off()
    
    ### hm by evalution category
    pdf(paste0('./plot/plot/all_eval_score_',names(select),'.pdf'),width=3.8,height=6)  
    print(ggplot(data=pd2, aes(x=evaluation,y=method,fill=score)) + geom_tile()+
      scale_fill_viridis(discrete=FALSE, na.value='white') +
      theme_hm(pd2$method) + xlab('') + ylab('')+
      labs(fill='score') + 
      ggtitle(names(select)))
    dev.off()
    
    if (names(select)=='all'){
      eff[rownames(eff)=='no_imp'] = NA
      effpd = melt(eff)
      colnames(effpd) = c('method','evaluation','score') 
      effpd$method = factor(as.character(effpd$method), levels=mtdorder)
      pdf(paste0('./plot/plot/efficiency.pdf'),width=0.8,height=6)  
      print(ggplot(data=effpd, aes(x=evaluation,y=method,fill=score)) + geom_tile()+
              scale_fill_viridis(discrete=FALSE, na.value='white') +
              theme_void() + xlab('') + ylab('') +
              theme(axis.text.y=element_blank(), axis.ticks = element_blank(),axis.text.x=element_text(face='bold',angle=90,hjust=1),legend.position = 'none')+
              ggtitle(''))
      dev.off()
    }
    # ## with sidebars    
    # library(pheatmap)
    # pheatmap(pd2, show_rownames=T, cluster_cols=F, cluster_rows=F, scale="none"),
    #          annotation_colors=my_color[2:14],
    #          cex=1, border_color=FALSE,annotation_row=eff,
    #          color = rev(morecols(20)),breaks=breaksList)
    # eff[1,] = 0
    # 
    # library(gplots)
    # library(RColorBrewer)
    # mypalette <- brewer.pal(11,"RdYlBu")
    # morecols <- colorRampPalette(mypalette)
    # morecols2 <- colorRampPalette(brewer.pal(6,"Set3"))
    # breaksList = seq(0, 1, by = 0.05)
    # eff[1,] <- 0
    # anno_color = lapply(colnames(eff), function(i) {
    #     n = length(eff[,i])
    #     col = morecols2(n)
    #     names(col) = eff[,i]
    #     col
    # })
    # eff <- list(time=eff[,1],memory=eff[,1],scalability=eff[,1])
    # tmp <- anno_color[[1]]
    # names(tmp) <- row.names(eff)
    # for (i in 1:3) anno_color[[i]] <- tmp
    # for (i in 1:ncol(eff)) eff[,i] <- row.names(eff)
    # names(anno_color) = colnames(eff)
    # pheatmap(tmpmat,scale = "none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=eff,
    #      treeheight_row = 0, treeheight_col = 0,annotation_colors=anno_color,
    #      legend=FALSE,show_rownames=TRUE,show_colnames=TRUE,
    #      # color=colorRampPalette(c(RColorBrewer::brewer.pal(3,"Spectral")[3],
    #      #                          RColorBrewer::brewer.pal(3,"Spectral")[2],
    #      #                          RColorBrewer::brewer.pal(3,"Spectral")[1]))(50),
    #      gaps_row=c(7,12),
    #      annotation_legend=FALSE,
    #      border_color="white",
    #      fontsize_row=10)
}


