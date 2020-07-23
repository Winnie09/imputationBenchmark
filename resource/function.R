library(Matrix)
process10x_rmDupGenes <- function(genebycellmat){
  tb = genebycellmat
  tb <- tb[rowSums(tb) > 0,]
  gn = rownames(tb)
  rs <- rowSums(tb)
  kid <- sapply(unique(gn),function(sid) {
    tmp <- which(gn==sid)
    if (length(tmp)==1) {
      tmp
    } else {
      tmp[which.max(rs[tmp])]
    }
  })
  if (is.list(kid)) kid = unlist(kid)
  tb <- tb[kid,]
  row.names(tb) <- gn[kid]
  tb <- tb[!grepl('^MT-',row.names(tb)),]
  tb = round(tb)
}


theme_hm <- function(method_vec) {
  theme_minimal() + 
    theme(axis.text.x = element_text(size=12,color='black',angle=90,hjust=1), 
          axis.text.y = element_text(size=12,color=ifelse(levels(method_vec)=='no_imp','red','black')),
          plot.title=element_text(size=12),
          legend.title=element_text(size=12))
}

theme_hm2 <- function(method_vec) {
  theme_classic() + 
    theme(axis.text.x = element_text(size=12,color='black',angle=90,hjust=1), 
          axis.text.y = element_text(size=12,color=ifelse(levels(method_vec)=='no_imp','red','black')),
          plot.title=element_text(size=12),
          legend.title=element_text(size=12))
}

plotfunc <- function(pd,napd=NULL,title='NULl',legend.lab = 'NULL') {
  fullpd <- rbind(pd,napd[,1:ncol(pd)])
  fullpd$method <- factor(as.character(fullpd$method),levels=c(levels(napd$method),levels(pd$method)))
  if (!is.null(napd)){
    if (length(unique(napd$NA.reason)) > 1){
      nav <- c('grey','black')
      names(nav) <- c('ImputationFail',setdiff(unique(napd$NA.reason),'ImputationFail'))  
    } else {
      nav = 'black'
    }
    if (dataset!='hca'){
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) + geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1),na.value = 'white') +
              theme_minimal() + xlab('') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=nav)+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              labs(fill=legend.lab)+
              ggtitle(title))    
    } else {
      print(ggplot() + geom_tile(data=fullpd,aes(x=data, y=method, fill=eval)) + geom_tile(data=napd,aes(x=data,y=method,color=NA.reason),fill='white',size=0.2) + 
              scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.95,1),na.value = 'white') +
              theme_minimal() + xlab('466 pairs of cell types') + ylab('') + theme_hm(fullpd$method) + 
              scale_color_manual(values=nav)+
              guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5,title.position = 'top'))+
              theme(legend.key.width=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'))+
              theme(axis.text.x=element_blank())+
              labs(fill=legend.lab)+
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
  namat <- setdiff(coldf[,2],unique(pd[,1]))
  if (length(namat) != 0){
    napd <- cbind(expand.grid(namat,unique(pd[,2])),NA)
    colnames(napd) <- c('method','data','eval')
    napd$NA.reason <- 'ImputationFail'
    for (i in unique(napd$method)){
      if (analysis.type == 'realDE'){
        reason.fail = 'DifferentialFail'
        result.path = paste0('./',analysis.type,'/result/',dataset,'/diff/',analysis.method,'/',i,'/res/')
        allf = sub('.rds','',list.files(paste0('./',analysis.type,'/result/',ifelse(dataset=='10x','cellbench',dataset),'/diff/',analysis.method,'/saver/res/')))
        existf <- sub('.rds','',list.files(result.path))
      } else if (analysis.type == 'nullDE') {
        reason.fail = 'DifferentialFail'
        result.path = paste0('./nullDE/result/',ifelse(dataset=='cellbench','10x',dataset),'/diff/',analysis.method,'/', i,'/res/')
        allf = sub('.rds','',list.files(paste0('./nullDE/result/',dataset,'/diff/',analysis.method,'/saver/res/')))
        existf <- sub('.rds','',list.files(result.path))
      } else if (analysis.type == 'clustering'){
        reason.fail = 'ClusteringFail'
        result.path = paste0('./clustering/perf/',dataset,'/',analysis.method,'/medianSil/') ## cellbench
        existf = NULL
        
      } else if (analysis.type == 'trajectory'){
        reason.fail  = 'TrajectoryFail'
      }  
      # existf <- sub('.rds','',list.files(result.path))
      dataset = ifelse(dataset=='cellbench','10xcellline',dataset)
      dataset = ifelse(dataset=='10x','10xcellline',dataset)
      impute.exist <- list.files(paste0('./result/procimpute/',dataset,'/',coldf[match(i,coldf$fullName),'shortName']))
      if (length(existf) > 0 ){
        napd[intersect(which(napd$method==i), match(setdiff(allf,existf),napd$data)),'NA.reason'] <- reason.fail
      } else if (length(impute.exist)>0){
        napd[which(napd$method==i),'NA.reason'] <- reason.fail
      }
    }  
  } else {
    napd = NULL
  }
  napd
}


