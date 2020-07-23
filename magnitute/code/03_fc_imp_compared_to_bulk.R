setwd('/dcl02/hongkai/data/whou/')
library(data.table)
source('./resource/function.R')
get_bulk_foldchange <- function(expr){
  expr <- expr[, sort(colnames(expr))]
  mat <- lapply(1:(ncol(expr)-1), function(i){
      tmp <- sapply(((i+1):ncol(expr)), function(j){
        expr[,i] - expr[,j]
      })
      if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
      colnames(tmp) <- paste0(colnames(expr)[i], '_', colnames(expr)[(i+1):ncol(expr)])
      tmp
  })
  return(do.call(cbind, mat))
}

get_sc_foldchange <- function(expr, ct){
  expr_mean <- aggregate(t(expr), list(ct), mean)
  rownames(expr_mean) <- expr_mean[,1]
  expr_mean = expr_mean[,-1]
  expr_mean = t(expr_mean)
  expr_mean = as.matrix((expr_mean))
  expr_mean = expr_mean[, sort(colnames(expr_mean))]
  mat <- lapply(1:(ncol(expr_mean)-1), function(i){
      tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
        expr_mean[,i] - expr_mean[,j]
      })
      if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
      colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
      tmp
  })
  return(do.call(cbind, mat))
}



g = fread('./resource/gencode.v19.annotation.gtf',data.table = F)
g <- g[g[,3]=='gene',]
gn <- gsub('\"','',sub(' gene_name ','',sapply(g[,9],function(i) strsplit(i,';')[[1]][5])))
gl <- g[,5]-g[,4]+1
names(gl) <- gn
gl <- gl/1000

bulk = readRDS('./rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
bulk <- bulk[row.names(bulk) %in% names(gl),]
bulk <- bulk/gl[row.names(bulk)]
lib <- colSums(bulk)/1e6
bulk <- t(t(bulk)/lib)
bulk <- log2(bulk + 1) ## TPM
colnames(bulk) = sub('_.*','',colnames(bulk))
bulk <- sapply(unique(colnames(bulk)),function(i) rowMeans(bulk[,colnames(bulk)==i]))


sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/magic/sc_10x_5cl.rds'))
cl = sub('.*:','',colnames(sexpr))
intgene = intersect(rownames(bulk),rownames(sexpr))

sexpr = sexpr[intgene,]
umi_magic = get_sc_foldchange(sexpr, cl)

bulk = bulk[intgene,]
bk <- get_bulk_foldchange(bulk)

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds'))
cl = sub('.*:','',colnames(sexpr))
sexpr = sexpr[intgene,]
umi_saver = get_sc_foldchange(sexpr, cl)

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds'))
cl = sub('.*:','',colnames(sexpr))
sexpr = sexpr[intgene,]
umi_raw = get_sc_foldchange(sexpr, cl)


library(reshape2)
pd_bk <- melt(bk)
pd_bk_umi <- cbind(pd_bk, method = 'bulk')
pd_saver <- melt(umi_saver)
pd_saver_umi <- cbind(pd_saver, method = 'saver')
pd_magic <- melt(umi_magic)
pd_magic_umi <- cbind(pd_magic, method = 'magic')
pd_raw <- melt(umi_raw)
pd_raw_umi <- cbind(pd_raw, method = 'no_Imp')

library(ggplot2)
library(RColorBrewer)
p_umi <- ggplot(data = rbind(pd_bk_umi, pd_magic_umi, pd_saver_umi, pd_raw_umi), aes(x = Var2, y = value, color = method)) + geom_boxplot(alpha = 0.3, outlier.colour = NA)  +  
  scale_color_brewer(palette="Set1") + 
  coord_flip() +
  ylim(c(0, 5)) +
  ggtitle('UMI') + 
  xlab('log fold change') + 
  ylab('')

# --------
# fluidigm
# --------

bulk = readRDS('./rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
colnames(bulk)[which(colnames(bulk)=='H1-hESC')] = 'H1'
colnames(bulk)[which(colnames(bulk)=='IMR-90')] = 'IMR90'

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds'))
int = intersect(rownames(bulk), rownames(sexpr))
bulk = bulk[int,]
flui_saver <- get_sc_foldchange(sexpr[int, ], ct = sub('_.*','',colnames(sexpr)))

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/magic/GSE81861_Cell_Line_COUNT.rds'))
flui_magic <- get_sc_foldchange(sexpr[int, ], ct = sub('_.*','',colnames(sexpr[int, ])))

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds'))
flui_raw <- get_sc_foldchange(sexpr[int, ], ct = sub('_.*','',colnames(sexpr[int, ])))


bk <- get_bulk_foldchange(bulk)
pd_bk <- melt(bk)
pd_bk <- cbind(pd_bk, method = 'bulk')
pd_saver <- melt(flui_saver)
pd_saver <- cbind(pd_saver, method = 'saver')
pd_magic <- melt(flui_magic)
pd_magic <- cbind(pd_magic, method = 'magic')
pd_raw <- melt(flui_raw)
pd_raw <- cbind(pd_raw, method = 'no_Imp')

p_flui <- ggplot(data = rbind(pd_bk, pd_magic, pd_saver, pd_raw), aes(x = Var2, y = value, color = method)) + geom_boxplot(alpha = 0.3, outlier.colour = NA)  +  
  scale_color_brewer(palette="Set1") + 
  coord_flip() + 
  ylim(c(0,5)) +
  ggtitle('Fluidigm') + 
  xlab('log fold change') + 
  ylab('')

library(gridExtra)
pdf('./rna_imputation/magnitute/plot/fc_imp_compared_to_bulk.pdf', height = 5, width = 10)
grid.arrange(p_umi, p_flui, nrow = 1)
dev.off()



png('./rna_imputation/magnitute/plot/fc_imp_compared_to_bulk_scatter_fluidigm.png', width = 600, height = 1000)
par(mfrow = c(5,3))
for (pair in unique(pd_bk[,2])){
  df = data.frame(bulk = pd_bk[pd_bk[,2] == pair, 3],
                  saver = pd_saver[pd_saver[,2] == pair, 3],
                  magic = pd_magic[pd_magic[,2] == pair, 3],
                  raw = pd_magic[pd_raw[,2] == pair, 3],
                  stringsAsFactors = FALSE)
  plot(df[,4] ~ df[,1], xlab = 'bulk LFC', ylab = 'no_Imp LFC', main = pair, pch = 20, cex = .5)
  plot(df[,2] ~ df[,1], xlab = 'bulk LFC', ylab = 'SAVER LFC', main = pair, pch = 20, cex = .5)
  plot(df[,3] ~ df[,1], xlab = 'bulk LFC', ylab = 'SAVER LFC', main = pair, pch = 20, cex = .5)
}
dev.off()
    

png('./rna_imputation/magnitute/plot/fc_imp_compared_to_bulk_scatter_umi.png', width = 600, height = 1000)
par(mfrow = c(5,3))
for (pair in unique(pd_bk_umi[,2])){
  df = data.frame(bulk = pd_bk_umi[pd_bk_umi[,2] == pair, 3],
                  saver = pd_saver_umi[pd_saver_umi[,2] == pair, 3],
                  bulk = pd_magic_umi[pd_magic_umi[,2] == pair, 3],
                  raw = pd_raw_umi[pd_raw_umi[,2] == pair, 3],
                  stringsAsFactors = FALSE)
  plot(df[,4] ~ df[,1], xlab = 'bulk LFC', ylab = 'no_Imp LFC', main = pair, pch = 20, cex = .5)
  plot(df[,2] ~ df[,1], xlab = 'bulk LFC', ylab = 'SAVER LFC', main = pair, pch = 20, cex = .5)
  plot(df[,3] ~ df[,1], xlab = 'bulk LFC', ylab = 'SAVER LFC', main = pair, pch = 20, cex = .5)
}
dev.off()














