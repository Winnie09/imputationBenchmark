# the data is 'mat'
## /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/code/make293t_simudata.R
numcell1 = as.numeric(commandArgs(trailingOnly = T)[1])
numcell2 = as.numeric(commandArgs(trailingOnly = T)[2])
probGene = as.numeric(commandArgs(trailingOnly = T)[3])
probRead= as.numeric(commandArgs(trailingOnly = T)[4])

fullmat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')
set.seed(12345)
sid <- sample(1:ncol(fullmat))
id1 <- sid[1:numcell1]
id2 <- sid[(1+numcell1):(numcell1 + numcell2)]
id <- as.vector(c(id1,id2))
mat <- fullmat[,colnames(fullmat)[id]]

colnames(mat)[1:numcell1] <- paste0(colnames(mat)[1:numcell1],':1')
colnames(mat)[(numcell1+1):(numcell1+numcell2)] <- paste0(colnames(mat)[(numcell1+1):(numcell1+numcell2)],':2')
add <- fullmat[,setdiff(colnames(fullmat),colnames(mat))]
gn <- which(rowSums(mat) > 0)
mat <- mat[gn,]
oriadd <- add[gn,]

sampid <- sample(1:ncol(oriadd))
sampc1id <- sampid[1:length(id1)]
sampc2id <- sampid[(1+length(id1)):(length(id1)+length(id2))]

add1 <- sapply(sampc1id,function(i) {
  if (probRead< 1) {
    ssv <- oriadd[,i]
    tarlib <- round(sum(ssv) * probRead)
    read <- rep(names(ssv),ssv)
    read <- sample(read,tarlib)
    v <- rep(0,length(ssv))
    names(v) <- names(ssv)
    tab <- table(read)
    v[names(tab)] <- as.vector(tab)
    v
  } else if (probRead== 1) {
    oriadd[,i]
  } else {
    rowSums(oriadd[,sample(colnames(oriadd),probRead)])
  }
})
# add1 <- add1[rowSums(add1) > 0,] ## for probRead = 0

add2 <- sapply(sampc2id,function(i) {
  if (probRead< 1) {
    ssv <- oriadd[,i]
    tarlib <- round(sum(ssv) * probRead)
    read <- rep(names(ssv),ssv)
    read <- sample(read,tarlib)
    v <- rep(0,length(ssv))
    names(v) <- names(ssv)
    tab <- table(read)
    v[names(tab)] <- as.vector(tab)
    v
  } else if (probRead== 1) {
    oriadd[,i]
  } else {
    rowSums(oriadd[,sample(colnames(oriadd),probRead)])
  }
})
## add2 <- add2[rowSums(add2) > 0,] ## for probRead = 0
## remove low wxpr genes
sampgn <- sample(intersect(row.names(add1),row.names(add2)),nrow(mat)*probGene)     
simumat <- mat
g <- sampgn[1:round(length(sampgn)/2)]
dg1 <- g
simumat[g,1:numcell1] <- simumat[g,1:numcell1] + add1[g,]      
g <- sampgn[(round(length(sampgn)/2)+1):length(sampgn)]
dg2 <- g
simumat[g,((numcell1+1):(numcell1+numcell2))] <- simumat[g,((numcell1+1):(numcell1+numcell2))] + add2[g,]
simumat <- simumat[rowMeans(simumat > 0) >= 0.01,]
## save files
system(paste0('mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead), showWarnings = FALSE)
saveRDS(dg1,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/diffgn1.rds'))
saveRDS(dg2,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/diffgn2.rds'))
saveRDS(simumat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/genebycell.rds'))
saveRDS(sampgn,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/diffgn.rds'))
write.table(simumat, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/genebycell.csv'), sep=",",quote=F)
saveRDS(t(simumat),file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/cellbygene.rds'))
write.table(t(simumat), file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/cellbygene.csv'), sep=",",quote=F)
