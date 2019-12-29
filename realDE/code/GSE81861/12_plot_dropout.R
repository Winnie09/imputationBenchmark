allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/mast/',mtd,'/res/'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.rds')
ct = sub('_.*','',colnames(raw))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/GSE81861/plot/dropout/prop_cell_expressed.pdf',height=25,width=25)
par(mfrow=c(10,9))
for (sct in unique(ct)){
  for (mtd in allmtd){
    d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
    hist(rowMeans(d[,ct==sct]>0), main=paste(mtd,sct), col='grey',braks=30,xlab='cell.expressed.in(%)',ylab='#gene')  
  }
}
  
dev.off()
