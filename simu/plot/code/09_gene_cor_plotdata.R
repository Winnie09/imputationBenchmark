library(reshape2)
library(ggplot2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/')
allmtd = setdiff(allmtd,'pblr')
pd = NULL
# allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/')
f = as.character(commandArgs(trailingOnly = T)[1])
  print(f)
  dg = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',sub('.rds','',f),'/diffgn.rds'))
  data2 = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f,'/genebycell.rds'))
  dg = dg[dg%in%row.names(data2)]
  pd2 = data.frame(melt( corfunc(t(data2[dg,]), t(data2[setdiff(row.names(data2),dg),])) ), type= 'raw', f = f)
  set.seed(12345)
  tmppd = NULL
  tmppd = rbind(tmppd, pd2[sample(1e3),])
  res <- lapply(allmtd, function(mtd){
    print(mtd)
    existf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/', mtd)))
    if (f %in% existf){
      data = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/', mtd, '/', f,'.rds'))
      dg = dg[dg%in%row.names(data)]
      pd1 = data.frame(melt( corfunc(t(data[dg,]), t(data[setdiff(row.names(data),dg),])) ),type=mtd, f = f)
      data.frame(pd1[sample(1:nrow(pd1),1e3),,drop=F])
    }else{
      NULL
    }
  })
  names(res) = allmtd
  res = do.call(rbind,res)
  if (!is.null(res)){
    tmppd = rbind(tmppd, res)
  } 
saveRDS(tmppd,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/gene_cor_plotdata/',f,'.rds'))
# 
# 
# pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/diffgene_nondiffgene_cor.pdf'),width=20,height=12)
# ggplot(data=pd,aes(x=value,color=type,fill=type))+geom_density(alpha=0.1) + xlab('PCC between Diffgene') 
# + theme_classic() + facet_wrap(~type)
# dev.off()  
# 

