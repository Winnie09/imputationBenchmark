clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')
ct <- ct[!is.na(ct[,2]),]

mtd = commandArgs(trailingOnly = T)[1]
suppressMessages(library(MAST))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/', mtd,'/res/'), showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))
if (grepl('A.1',colnames(imp)[1])){
  colnames(imp) = sub('\\.','-',colnames(imp))
}
cid <- expand.grid(ct[,1],ct[,1])
cid <- cid[cid[,1] > cid[,2],]
cid <- cid[ct[match(cid[,1],ct[,1]),2]!=ct[match(cid[,2],ct[,1]),2],]
for (i in 1:nrow(cid)){
  clu1 <- cid[i,1]
  clu2 <- cid[i,2]
  pval <- sapply(1:nrow(imp), function(i) {
    v1 = imp[i, clu[clu[,2]==clu1,1]]
    v2 = imp[i, clu[clu[,2]==clu2,1]]
    if (length(unique(c(v1,v2)))==1){
      1
    } else {
      wilcox.test(v1, v2)$p.value  
    }
  })
  fdrv <- p.adjust(pval,method='fdr')
  res = data.frame(Gene=rownames(imp), fdr = fdrv)
  saveRDS(res, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/', mtd,'/res/',ct[ct[,1]==clu1,2],'-',clu1,'_',ct[ct[,1]==clu2,2],'-',clu2,'.rds'))  
}

