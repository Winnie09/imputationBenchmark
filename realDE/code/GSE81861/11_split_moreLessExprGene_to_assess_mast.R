allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.rds')
ct = sub('_.*','',colnames(raw))
ove <- sapply(allmtd, function(mtd){
  print(mtd)
  t(sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/',f))){
      res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/',f))
      res = res[order(res[,'pvalue'], res[,'stat']),]
      ct1 <- sub('_.*','',f)
      ct2 <- sub('.rds*','',sub('.*_','',f))
      tmpmat = raw[,ct%in%c(ct1,ct2)]
      bfex <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/bulkdiff/')
      bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/bulkdiff/',sub('.rds','',f),'_diffgene.rds'))
      fc <- rowMeans(raw[,ct == ct1]) - rowMeans(raw[,ct == ct2])
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
}) ## first 10 high, last 10 low
ove_high = t(ove[1:10,])
ove_low = t(ove[11:20,])
colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/overlap/mast/'
dir.create(rdir,showWarnings = F, recursive = T)
saveRDS(ove_high, paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
saveRDS(ove_low, paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))

