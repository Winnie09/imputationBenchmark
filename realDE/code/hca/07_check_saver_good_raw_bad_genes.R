allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/',mtd,'/res/'))

ranklist <- lapply(allf,function(f){
  ct1 <- sub('-.*','',f)
  ct2 <- sub('-.*','',sub('.*_','',f))
  bfex <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/bulkdiff/')
  bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
  bf <- intersect(bfex,bf)
  gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/bulkdiff/',bf))
  
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/raw/res/',f))
  res = res[order(res$fdr),]
  rawrank <- 1:nrow(res)
  names(rawrank) <- res[,"Gene"]
  rawrank <- rawrank[intersect(names(rawrank),gs)]
  
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/wilcox/saver/res/',f))
  res = res[order(res[,'fdr']),]
  srank <- 1:nrow(res)
  names(srank) <- res[,"Gene"]
  srank <- srank[intersect(names(srank),gs)]
  rankdiff <- rawrank-srank[names(rawrank)] # negative: raw good, positive: saver good
})
names(ranklist)=allf

res <- sapply(allf[1:20],function(f){
  print(f)  
  imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/raw/MantonBM6.rds'))
  i = ranklist[[f]]
  g=names(sort(i[i>0],decreasing=T))[1:100] ## saver good
  g2=names(sort(i[i<0]))[1:100] ## raw good
  c(mean(rowSums(imp[g,,drop=F])/rowSums(imp[g,]>0) ),mean(rowMeans(imp[g,,drop=F]>0)),
    mean(rowSums(imp[g2,,drop=F])/rowSums(imp[g2,]>0)),mean(rowMeans(imp[g2,,drop=F]>0)))
})

