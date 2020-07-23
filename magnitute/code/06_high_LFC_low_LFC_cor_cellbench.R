# ----------------------------------------------------------------
# use high- and low- LFC genes seperately (in bulk and sc overlap)
# ----------------------------------------------------------------
setwd('/dcl02/hongkai/data/whou')
library(data.table)
source('./resource/function.R')
method = commandArgs(trailingOnly = T)[1]
# method='magic'
bulk = readRDS('./rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
g = fread('./resource/gencode.v19.annotation.gtf',data.table = F)
g <- g[g[,3]=='gene',]
gn <- gsub('\"','',sub(' gene_name ','',sapply(g[,9],function(i) strsplit(i,';')[[1]][5])))
gl <- g[,5]-g[,4]+1
names(gl) <- gn
gl <- gl/1000
bulk <- bulk[row.names(bulk) %in% names(gl),]
bulk <- bulk/gl[row.names(bulk)]
lib <- colSums(bulk)/1e6
bulk <- t(t(bulk)/lib)
bulk <- log2(bulk + 1) ## TPM
colnames(bulk) = sub('_.*','',colnames(bulk))
bulk <- sapply(unique(colnames(bulk)),function(i) rowMeans(bulk[,colnames(bulk)==i]))
sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/',method,'/sc_10x_5cl.rds'))
example = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds'))
colnames(sexpr) = colnames(example)
cl = sub('.*:','',colnames(sexpr))
intgene = intersect(rownames(bulk),rownames(sexpr))
bulk = bulk[intgene,]
sexpr = sexpr[intgene,]

get_scCellType_bulkCellType_cor_high <- function(ct1, ct2, bulkDiff, pairid){
  imp1 = sexpr[, which(cl==ct1)]
  imp2 = sexpr[, which(cl==ct2)]
  corvec = NULL
  cnt = 0
  for (i in 1:ncol(imp1)){
    for (j in 1:ncol(imp2)){
      cnt = cnt + 1
      if (cnt == pairid ){
        v <- abs(imp1[,i] - imp2[,j])
        id = names(v[v > quantile(v, 0.9)])
        return(cor((imp1[id,i] - imp2[id,j]), bulkDiff[id],method='spearman'))
        break
      }  
    }
  }
}

get_scCellType_bulkCellType_cor_low <- function(ct1, ct2, bulkDiff, pairid){
  imp1 = sexpr[, which(cl==ct1)]
  imp2 = sexpr[, which(cl==ct2)]
  corvec = NULL
  cnt = 0
  for (i in 1:ncol(imp1)){
    for (j in 1:ncol(imp2)){
      cnt = cnt + 1
      if (cnt == pairid ){
        v <- abs(imp1[,i] - imp2[,j])
        id = names(v[v < quantile(v, 0.1)])
        return(cor((imp1[id,i] - imp2[id,j]), bulkDiff[id],method='spearman'))
        break
      }  
    }
  }
}


resdiff <- readRDS(paste0('./rna_imputation/diff/result/cellbench_10x_5cl/', method, '.rds'))
v = sapply(1:(ncol(bulk)-1), function(i){
  sapply((i+1):ncol(bulk), function(j){
    cn = paste0(colnames(bulk)[i],'_',colnames(bulk)[j])
    rescn <- resdiff[[cn]]
    if (length(rescn)%%2){
      pairid <- which(rescn == median(rescn))
    } else {
      pairid <- which(rescn == sort(rescn)[ceiling(length(rescn)/2)])
    }
    tmp_high <- get_scCellType_bulkCellType_cor_high(ct1=colnames(bulk)[i], ct2=colnames(bulk)[j], bulkDiff = bulk[,colnames(bulk)[i]]-bulk[,colnames(bulk)[j]], pairid = pairid)
    tmp_low <- get_scCellType_bulkCellType_cor_low(ct1=colnames(bulk)[i], ct2=colnames(bulk)[j], bulkDiff = bulk[,colnames(bulk)[i]]-bulk[,colnames(bulk)[j]], pairid = pairid)
    c(cn, tmp_high, tmp_low)
  })
})

m = t(do.call(cbind, v))
colnames(m) <- c('celltype', 'high_LFC_SCC', 'low_LFC_SCC')
dir.create('./rna_imputation/magnitute/result/sc_10x_5cl/',recursive = T, showWarnings = F)
saveRDS(m, paste0('./rna_imputation/magnitute/result/sc_10x_5cl/',method,'.rds'))


