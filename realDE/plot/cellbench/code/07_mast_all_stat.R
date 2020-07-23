library(ggplot2)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
allmtd = list.files('./realDE/result/cellbench/diff/mast/')
mtd = allmtd[1]
allf = list.files(paste0('./realDE/result/cellbench/diff/mast/', mtd,'/res_allstat/'))
f = allf[1]
for (f in allf){
  print(f)
  res <- lapply(allmtd, function(mtd){
    if(file.exists(paste0('./realDE/result/cellbench/diff/mast/', mtd,'/res_allstat/',f))){
      tmp = readRDS(paste0('./realDE/result/cellbench/diff/mast/', mtd,'/res_allstat/',f))
      df = data.frame(tmp[,-1], data=sub('.rds','',f), method=mtd)  
      set.seed(12345)
      df = df[sample(1:nrow(df), 1e3), ]
    } else {
      print(mtd)
    }
  })
  
  pd = do.call(rbind,res)
  pd[,1] = as.numeric(pd[,1])
  pd[,2] = as.numeric(pd[,2])
  pd[,3] = as.numeric(pd[,3])
  pd = cbind(pd, coef_se = pd$stat/pd$z)
  pd = pd[!is.na(pd$method),]
  coldf=readRDS('./resource/method_latent_color.rds')
  pd$method = as.factor(coldf[match(as.character(pd$method), coldf[,'shortName']), 'fullName'])
  pd$pvalue = -log10(pd$pvalue+10^(-150))
  library(gridExtra)
  pdf(paste0('./realDE/plot/cellbench/plot/mast/',sub('.rds','',f),'_all_stat.pdf'),width=10, height=5)
  p1 <- ggplot() + geom_violin(data=pd, aes(x=method, y = pvalue, fill=method),alpha=0.3, scale='width') +
    theme_classic()+theme(legend.position = 'none')+
    theme(axis.text=element_text(color='black',angle=90,hjust=1), text=element_text(color='black'))+
    ggtitle(expression(paste(-log[10],'(pvalue + ', 10^{-150}, ')' ))) + xlab('')+ylab('')
    # scale_y_continuous(trans='log10')
  # dev.off()
  
  p2 <- ggplot() + geom_violin(data=pd, aes(x=method, y = stat, fill=method),alpha=0.3, scale='width') +
    theme_classic()+theme(legend.position = 'none')+
    theme(axis.text=element_text(color='black',angle=90,hjust=1), text=element_text(color='black'))+
    ggtitle('coefficient')+ xlab('') + ylab('')
  
  p3 <- ggplot() + geom_violin(data=pd, aes(x=method, y = coef_se, fill=method),alpha=0.3, scale='width') +
    theme_classic()+theme(legend.position = 'none')+
    theme(axis.text=element_text(color='black',angle=90,hjust=1), text=element_text(color='black'))+
    ggtitle('standard error of coefficient')+ xlab('')+ylab('')
  
  
  p4 <- ggplot() + geom_violin(data=pd, aes(x=method, y = z, fill=method),alpha=0.3, scale='width') +
    theme_classic()+theme(legend.position = 'none')+
    theme(axis.text=element_text(color='black',angle=90,hjust=1), text=element_text(color='black'))+
    ggtitle('z-score(coefficient divided by coefficient standard error)')+ xlab('')+ylab('')
  
  grid.arrange(p2,p3,p4,p1,nrow=2)
  dev.off()  
}

