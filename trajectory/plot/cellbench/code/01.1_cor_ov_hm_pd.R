trajmtd = commandArgs(trailingOnly = T)[1]
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
ddir = paste0('./trajectory/result/cellbench/',trajmtd,'/cor_ov/')

allf = list.files(ddir)
allf = setdiff(allf,'zinbwave.rds')
res <- lapply(allf,function(f){
  data = readRDS(paste0(ddir,f))
  data = matrix(unlist(data),ncol=2,byrow=T)
  data.frame(method = sub('.rds','',f), data = c("cellmix1.rds","cellmix2.rds","cellmix3.rds","cellmix4.rds","RNAmix_celseq2.rds","RNAmix_sortseq.rds"), 
             cor = data[,1], ov = data[,2], stringsAsFactors = F)
})
df = do.call(rbind,res)
df[,'data'] = sub('.rds','',df[,'data'] )


mtdorder1 = names(sort(tapply(df[,'cor'], list(df[,'method']), mean, na.rm=T), decreasing = F))
stat = tapply(df[,'cor'], list(df[,'method']), mean, na.rm=T)
mtdorder1 = c(setdiff(unique(df$method), mtdorder1), mtdorder1)
saveRDS(mtdorder1,paste0('./result/perf/rank/trajectory_cor_',trajmtd,'.rds'))
saveRDS(stat,paste0('./result/perf/assess/trajectory_cor_',trajmtd,'.rds'))
mtdorder2 = names(sort(tapply(df[,'ov'], list(df[,'method']), mean, na.rm=T), decreasing = F))
stat = tapply(df[,'ov'], list(df[,'method']), mean, na.rm=T)
mtdorder2 = c(setdiff(unique(df$method), mtdorder2), mtdorder2)
saveRDS(mtdorder2,paste0('./result/perf/rank/trajectory_ov_',trajmtd,'.rds'))
saveRDS(stat,paste0('./result/perf/assess/trajectory_ov_',trajmtd,'.rds'))

### use methods fullName
# df = df[!df[,'method']%in%names(which(tapply(is.na(df[,'cor']),list(df[,'method']),sum)==6)),]
colordf = readRDS('./resource/method_latent_color.rds')
df$method = colordf[match(df$method,colordf$shortName),'fullName']
df1 = df
df1[,'method'] = factor(df1$method, levels = colordf[match(mtdorder1,colordf$shortName),'fullName'])
df2 = df
df2[,'method'] = factor(df2$method, levels =colordf[match(mtdorder2,colordf$shortName),'fullName'])

### add incomplte methods to heatmap plotdata
addIncomMtd <- function(df){
  mtdorder = levels(df$method)
  add <- lapply(setdiff(colordf$fullName,unique(df$method)), function(i){
    data.frame(method = i, data=unique(df$data), cor=NA, ov = NA, stringsAsFactors = F)
  })
  add = do.call(rbind, add)
  df = rbind(df,add)
  df$method = factor(as.character(df$method), levels = c(setdiff(colordf$fullName, mtdorder), mtdorder))
  df  
}

pd1 <- addIncomMtd(df1)
pd2 <- addIncomMtd(df2)
saveRDS(list(pd1=pd1,pd2=pd2), paste0('./trajectory/plot/cellbench/plot/',trajmtd,'/cor_ov_hm_pd.rds'))
