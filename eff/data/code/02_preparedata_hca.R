dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/data/processed/'
# for (num.cell in c(5e4, 1e5)){
  num.cell = 1e5
  dir.create(paste0(dir,num.cell),showWarnings = F)
  d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/hca/rawcount.rds')
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
    tb <- tb[kid,]
    row.names(tb) <- gn[kid]
    tb <- tb[!grepl('^MT-',row.names(tb)),]
    tb = round(tb)
  }
  set.seed(12345)
  if (num.cell < ncol(d)){
    id = sample(ncol(d), num.cell)  
  } else {
    id = sample(ncol(d), num.cell, replace = T)  
  }
  
  
  mat = process10x_rmDupGenes(d[,id])
  
  ## remove low expr genes
  # hist(rowMeans(mat>0), breaks=100)
  
  kid = which(rowMeans(mat > 0) >=0.05)
  mat = mat[kid,]
  mat = round(mat)
  
  saveRDS(mat, paste0(dir, num.cell,'/genebycell.rds'))      
  write.table(mat, file=paste0(dir,num.cell,'/genebycell.csv'), sep=',', quote=F)
  tm = t(mat)
  saveRDS(tm, paste0(dir,num.cell,'/cellbygene.rds'))
  write.table(tm, paste0(dir, num.cell,'/cellbygene.csv'), sep=',', quote=F)
  
  suppressMessages(library(scran))
  sce <- SingleCellExperiment(list(counts=mat))
  if (ncol(mat) < 21){
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
  } else {
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))
  }
  sf <- sizeFactors(sce)
  normmat <- sweep(mat,2,sf,'/')
  normmat <- log2(normmat + 1)
  saveRDS(normmat, paste0(dir, num.cell,'/norm_genebycell.rds'))      
  write.table(normmat, paste0(dir,num.cell,'/norm_genebycell.csv'), sep=',', quote=F)
  tm_norm = t(normmat)
  saveRDS(tm_norm, paste0(dir,num.cell,'/norm_cellbygene.rds'))
  write.table(tm_norm, paste0(dir, num.cell, '/norm_cellbygene.csv'), sep=',', quote=F)      
# }

