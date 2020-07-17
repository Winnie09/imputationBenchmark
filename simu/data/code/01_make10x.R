#get fullmat
fullmat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
ct = sub('_.*','',colnames(fullmat))
fullmat = fullmat[,which(ct == '293T')]
for (cellnum in 1e2*(1:5)) {
  set.seed(12345)
  allc <- sample(colnames(fullmat),cellnum*2)
  c1 <- allc[1:cellnum]
  c2 <- allc[(cellnum+1):(cellnum*2)]
  
  mat <- fullmat[,c(c1,c2)]
  add <- fullmat[,setdiff(1:ncol(fullmat),allc)]
  gn <- which(rowSums(mat) > 0)
  mat <- mat[gn,]
  add <- add[gn,]
  add <- add[rowSums(add) > 0,]
  fulladd <- add
  
  for (per in c(0.3,0.2,0.1,0.05)) {
    sampgn <- sample(row.names(add),nrow(mat)*per)
    simumat <- mat
    add <- fulladd[,sample(1:ncol(fulladd),cellnum)]
    for (g in sampgn[1:round(length(sampgn)/2)]) {
      simumat[g,c1] <- simumat[g,c1] + add[g,]
    }
    for (g in sampgn[(round(length(sampgn)/2)+1):length(sampgn)]) {
      simumat[g,c2] <- simumat[g,c2] + add[g,]
    }
    system(paste0('mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per))
    saveRDS(simumat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/genebycell.rds')) ## count
    write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/genebycell.csv'), sep=",",quote=F)
    
    saveRDS(sampgn,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/diffgn.rds'))
    saveRDS(t(simumat),file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/cellbygene.rds'))
    write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/cellbygene.csv'), sep=",",quote=F)
    
    libsize = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/libszie.rds')
    libsize = libsize[match(colnames(simumat),names(libsize))]
    normmat <- sweep(simumat,2,libsize,'/')
    normmat <- log2(normmat + 1)
    saveRDS(normmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/norm_genebycell.rds'))
    write.table(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/norm_genebycell.csv'),sep=",",quote=F)
    tnormmat = t(normmat)
    saveRDS(tnormmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/norm_cellbygene.rds'))
    write.table(tnormmat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/',cellnum,"_",per,'/norm_cellbygene.csv'), sep=",",quote=F)
  }
}

