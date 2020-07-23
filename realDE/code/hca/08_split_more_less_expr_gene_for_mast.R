allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/',mtd,'/res/'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/raw/MantonBM6.rds')
ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/ct/MantonBM6.rds')
library(parallel)
ove <- mclapply(allmtd, function(mtd){
  print(mtd)
  t(sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/',mtd,'/res/',f))){
      res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/',mtd,'/res/',f))
      res$pvalue = p.adjust(res$pvalue,method='fdr')
      res = res[order(res[,'pvalue'], -abs(res[,'stat'])),]
      ct1 <- sub('-.*','',sub('_.*','',f))
      ct2 <- sub('-.*','',sub('.*_','',f))
      tmpmat = raw[,ct%in%c(ct1,ct2)]
      
      #highid = which(rowMeans(tmpmat>0)>quantile(rowMeans(tmpmat>0),0.9))## more cell express
      #lowid = which(rowMeans(tmpmat>0)<quantile(rowMeans(tmpmat>0),0.1))## less
      bfex <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/bulkdiff/')
      bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/bulkdiff/',bf))
      fc <- rowMeans(raw[,which(ct == ct1)]) - rowMeans(raw[,which(ct == ct2)])
      fc <- abs(fc)
      highid = names(which(fc>quantile(fc,0.9)))## express in more cells
      lowid = names(which(fc<quantile(fc,0.1)))## less
      tmpres <- res[res[,1] %in% highid,]
      tmp1 <- mean(sapply(c(1:100)*10,function(i) {
        mean(tmpres[1:i,1] %in% gs)  ## discuss
      }))
      tmpres <- res[res[,1] %in% lowid,]
      tmp2 <- mean(sapply(c(1:100)*10,function(i) {
        mean(tmpres[1:i,1] %in% gs)  ## discuss
      }),na.rm=T)
      c(tmp1,tmp2)
    } else {
      return(c(NA,NA))
    }
  }))
},mc.cores=10) ## first 10 high, last 10 low


rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/overlap/mast/'
dir.create(rdir,showWarnings = F, recursive = T)
saveRDS(ove, paste0(rdir,'bulk_sc_diffgene_overlaps_fcGene.rds'))
if (is.list(ove){
  ove_high <- t(sapply(ove,function(i){
    i[,1]
  }))
  ove_low <- t(sapply(ove,function(i){
    i[,2]
  }))
  colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
  rownames(ove_high) <- rownames(ove_low) <- allmtd
  
}) else {
  ove_high = t(ove[1:(nrow(ove)/2),])
  ove_low = t(ove[(nrow(ove)/2+1):nrow(ove),])
  colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
}
saveRDS(ove_high, paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
saveRDS(ove_low, paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))
