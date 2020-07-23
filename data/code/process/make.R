# /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/process/make.sh
library(data.table)
allf <- list.files("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/500more_dge")
for (f in allf) {
  data <- fread(paste0("zcat /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/500more_dge/",f),data.table = F)
  row.names(data) <- data[,1]
  data <- as.matrix(data[,-1])
  #      data <- data[!grepl('mt-',row.names(data)),]
  #      data <- data[,colSums(data > 0) >= 500]
  sf <- sub('_.*','',f)
  libsize = colSums(data)/1000
  normdata <- sweep(data,2,libsize,'/')
  normdata <- log2(normdata + 1)
  saveRDS(normdata, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,"/norm_genebycell.rds"))
  #      system(paste0('mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf))
  #      saveRDS(data, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,"/genebycell.rds"))
  #      write.table(data,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,'/genebycell.csv'),sep=",",quote=F)
  #      data = t(data)
  #      saveRDS(data, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,'/cellbygene.rds'))
  #      write.table(data, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,'/cellbygene.csv'),sep=",", quote=F) ## rows cells, cols genes
  
  #      write.table(data, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,'/cellbygene_noname.csv'),sep=",", quote=F,row.names=F,col.names=F) ## rows cells, cols genes
  
  #      ct = as.matrix(rep(0,nrow(data)))
  #      write.table(ct, file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/500more_dge',sf,'/ct.csv'),sep=",",  col.names=FALSE, quote=F,row.names = F)
}
