library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
method = commandArgs(trailingOnly = T)[1]
# method='magic'
bulk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
g = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/gencode.v19.annotation.gtf',data.table = F)
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
sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',method,'/sc_10x_5cl.rds'))
example = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds'))
colnames(sexpr) = colnames(example)
cl = sub('.*:','',colnames(sexpr))
intgene = intersect(rownames(bulk),rownames(sexpr))
bulk = bulk[intgene,]
sexpr = sexpr[intgene,]

get_scCellType_bulkCellType_cor <- function(ct1, ct2, bulkDiff){
  imp1 = sexpr[, which(cl==ct1)]
  imp2 = sexpr[, which(cl==ct2)]
  corvec = NULL
  corvec <- sapply(1:ncol(imp1),function(i) {
    print(i)
    sapply(1:ncol(imp2), function(j) {
      cor((imp1[,i] - imp2[,j]), bulkDiff,method='spearman')
    })
  })
  as.vector(corvec)
}

v = sapply(1:(ncol(bulk)-1), function(i){
  sapply((i+1):ncol(bulk), function(j){
    cn = paste0(colnames(bulk)[i],'_',colnames(bulk)[j])
    tmp <- get_scCellType_bulkCellType_cor(ct1=colnames(bulk)[i], ct2=colnames(bulk)[j], bulkDiff = bulk[,colnames(bulk)[i]]-bulk[,colnames(bulk)[j]])
    names(tmp)=cn
    tmp
  })
})

for (i in 1:length(v)){
  if (!is.list(v[[i]])) v[[i]] = list(as.vector(v[[i]]))
}
cn = NULL
for (i in 1:(ncol(bulk)-1)){
  for (j in (i+1):ncol(bulk)){
    cn = c(cn,paste0(colnames(bulk)[i],'_',colnames(bulk)[j]))
  }
}
for (i in 1:length(v)){
    names(v[[i]]) = cn[1:length(v[[i]])]
    if (i!=length(v)) cn = cn[(length(v[[i]])+1):length(cn)]  
}
tmp = c(v[[1]],v[[2]],v[[3]],v[[4]])

dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/diff/result/sc_10x_5cl/',recursive = T, showWarnings = F)
saveRDS(tmp,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/diff/result/sc_10x_5cl/',method,'.rds'))




