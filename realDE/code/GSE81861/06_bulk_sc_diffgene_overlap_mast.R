allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/',f))){
      res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/',f))
      res$pvalue = p.adjust(res$pvalue,method='fdr')      
      res = res[order(res[,'pvalue'], -abs(res[,'stat'])),]
      ct1 <- sub('_.*','',f)
      ct2 <- sub('.rds','',sub('.*_','',f))
      bfex <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/bulkdiff/')
      bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/bulkdiff/',bf))
      tmp <- mean(sapply(c(1:100)*10,function(i) {
        mean(res[1:i,1] %in% gs)  ## discuss
      }))
    } else {
      return(NA)
    }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/overlap/mast/'
dir.create(rdir,showWarnings = F,recursive = T)
saveRDS(ove, paste0(rdir,'bulk_sc_diffgene_overlaps.rds'))

