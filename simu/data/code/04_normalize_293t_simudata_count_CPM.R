# the data is 'mat'
## /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/code/make293t_simudata.R
numcell1 = as.numeric(commandArgs(trailingOnly = T)[1])
numcell2 = as.numeric(commandArgs(trailingOnly = T)[2])
probGene = as.numeric(commandArgs(trailingOnly = T)[3])
probRead= as.numeric(commandArgs(trailingOnly = T)[4])

simumat = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/genebycell.rds'))
libsize = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/libszie.rds')
cn = sub(':1','',colnames(simumat))
cn = sub(':2','',cn)
libsize = libsize[match(cn,names(libsize))]
normmat <- sweep(simumat,2,libsize,'/')
normmat <- log2(normmat + 1)
saveRDS(normmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_genebycell.rds'))
write.table(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_genebycell.csv'),sep=",",quote=F)
tnormmat = t(normmat)
saveRDS(tnormmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_cellbygene.rds'))
write.table(tnormmat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_cellbygene.csv'), sep=",",quote=F)

