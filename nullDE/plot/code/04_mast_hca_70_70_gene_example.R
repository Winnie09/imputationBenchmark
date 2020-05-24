dataset = 'hca'
allmtd = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/wilcox/'))
af = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/wilcox/raw/res/'))
f = "70_70.rds"
rawres = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/mast/raw/res/',f))
fdr = p.adjust(rawres[,'pvalue'],method='fdr')
g = rawres[which.min(fdr), 'Gene']

  mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
  cluct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')
  cluct = cluct[complete.cases(cluct),]
  cluct_s = paste0(cluct[,1],'_',cluct[,2])
  mtch_s = paste0(mtch$cluster, '_', mtch$celltype)
  select_cluct = names(sort(table(mtch_s[mtch_s%in%cluct_s]),decreasing=T)[1])
  select_cell = mtch[mtch$cluster == sub('_.*','',select_cluct) & mtch$celltype == sub('.*_','',select_cluct), 'cell']
  
pd <- lapply(allmtd,function(mtd){
  print(mtd)
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))){
    imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))
    if (grepl('A.1',colnames(imp)[1])){
      colnames(imp) = sub('\\.','-',colnames(imp))
    }
    n1 = as.numeric(sub('_.*','',f))
    n2 = as.numeric(sub('.rds','',sub('.*_','',f)))
    imp = imp[,select_cell]
    
    set.seed(12345)
    id = sample(1:ncol(imp), n1+n2)
    x = imp[g,id[1:n1]]
    y = imp[g,id[(n1+1):(n1+n2)]]
    data.frame(method=mtd,x=x,y=y)  
  }
})

pd = do.call(rbind,pd)
fdr = sapply(allmtd,function(mtd){
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/mast/',mtd,'/res/',f))){
    tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/',dataset,'/diff/mast/',mtd,'/res/',f))
    fdr = p.adjust(tmp[,'pvalue'],method='fdr')
    fdr[tmp[,'Gene']==g]
  }
  
})

if (is.list(fdr)) fdr = unlist(fdr)
levels(pd$method) = paste0(levels(pd$method), ';', round(as.numeric(fdr[levels(pd$method)]),2))
library(ggplot2)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/plot/plot/hca_mast_',g,'_70_70.pdf'),height=12,width=12)
ggplot(data=pd,aes(x=x,y=y)) + geom_point(size=0.5)+ 
  geom_abline(slope=1,intercept = 0, color='steelblue', linetype='dashed')+
  facet_wrap(~method,scales = 'free') + theme_classic()
dev.off()
