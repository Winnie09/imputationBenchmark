  source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
  bulk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
  sexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')
  cl = sub('_.*','',colnames(sexpr))
  intgene = intersect(rownames(bulk),rownames(sexpr))
  bulk = bulk[intgene,]
  sexpr = sexpr[intgene,]
  colnames(bulk)[which(colnames(bulk)=='H1-hESC')] = 'H1'
  colnames(bulk)[which(colnames(bulk)=='IMR-90')] = 'IMR90'
  
get_scCellType_bulkCellType_cor <- function(ct1, ct2, bulkDiff){
  imp1 = sexpr[, which(cl==ct1)]
  imp2 = sexpr[, which(cl==ct2)]
  corvec = NULL
  corvec <- sapply(1:ncol(imp1),function(i) {
    sapply(1:ncol(imp2), function(j) {
      cor((imp1[,i] - imp2[,j]), bulkDiff,method='spearman')
    })
  })
  median(as.vector(corvec))
}

v = sapply(1:(ncol(bulk)-1), function(i){
  sapply((i+1):ncol(bulk), function(j){
    print(paste0(colnames(bulk)[i],'_',colnames(bulk)[j]))
    get_scCellType_bulkCellType_cor(ct1=colnames(bulk)[i], ct2=colnames(bulk)[j], bulkDiff = bulk[,colnames(bulk)[i]]-bulk[,colnames(bulk)[j]])
  })
})

saveRDS(v,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/diff/result/raw.rds'))


