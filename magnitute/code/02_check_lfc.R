rm(list=ls())
# -------------
# load all data
# -------------
flui <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds')
str(flui)
umi <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds')
flui_magic <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/magic/GSE81861_Cell_Line_COUNT.rds')
umi_magic <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/magic/sc_10x_5cl.rds')
flui_saver <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds')
umi_saver <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds')

# -------------
# no imputation
# -------------
expr = flui
ct = sub('_.*', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
m_flui = do.call(cbind, mat)

expr = umi
ct = sub('.*:', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
str(mat)
m_umi = do.call(cbind, mat)

library(reshape2)
pd_flui <- melt(m_flui)
pd_flui_noImp <- cbind(pd_flui, platform = 'fluidigm')
pd_umi <- melt(m_umi)
pd_umi_noImp <- cbind(pd_umi, platform = 'umi')
library(ggplot2)
p_noImp <- ggplot(data = rbind(pd_flui_noImp, pd_umi_noImp)) + 
  geom_boxplot(aes(x = Var2, y = value, color = platform), outlier.color = NA) + 
  coord_flip() +
  xlab('') + ylab('log fold change of each gene') +
  ggtitle('no_imp')

# -------------
# MAGIC imputed
# -------------

expr = flui_magic
ct = sub('_.*', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
m_flui = do.call(cbind, mat)

expr = umi_magic
ct = sub('.*:', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
str(mat)
m_umi = do.call(cbind, mat)


pd_flui <- melt(m_flui)
pd_flui_magic <- cbind(pd_flui, platform = 'fluidigm')
pd_umi <- melt(m_umi)
pd_umi_magic <- cbind(pd_umi, platform = 'umi')
p_magic <- ggplot(data = rbind(pd_flui_magic, pd_umi_magic)) + 
  geom_boxplot(aes(x = Var2, y = value, color = platform), outlier.color = NA) + 
  coord_flip() +
  xlab('') + ylab('log fold change of each gene') +
  ggtitle('MAGIC')
  
# ----------------
# SAVER imputation
# ----------------

expr = flui_saver
ct = sub('_.*', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
m_flui = do.call(cbind, mat)

expr = umi_saver
ct = sub('.*:', '', colnames(expr))
expr_mean <- aggregate(t(expr), list(ct), mean)
rownames(expr_mean) <- expr_mean[,1]
expr_mean = expr_mean[,-1]
expr_mean = t(expr_mean)
expr_mean = as.matrix((expr_mean))
mat <- lapply(1:(ncol(expr_mean)-1), function(i){
    tmp <- sapply(((i+1):ncol(expr_mean)), function(j){
      expr_mean[,i] - expr_mean[,j]
    })
    if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
    colnames(tmp) <- paste0(colnames(expr_mean)[i], '_', colnames(expr_mean)[(i+1):ncol(expr_mean)])
    tmp
})
str(mat)
m_umi = do.call(cbind, mat)


pd_flui <- melt(m_flui)
pd_flui_saver <- cbind(pd_flui, platform = 'fluidigm')
pd_umi <- melt(m_umi)
pd_umi_saver <- cbind(pd_umi, platform = 'umi')
p_saver <- ggplot(data = rbind(pd_flui_saver, pd_umi_saver)) + 
  geom_boxplot(aes(x = Var2, y = value, color = platform), outlier.color = NA) + 
  coord_flip() +
  xlab('') + ylab('log fold change of each gene') +
  ggtitle('SAVER')

pdf('/Users/wenpinhou/Dropbox/rna_imputation/magnitute/plot/fc.pdf', width = 12, heigh = 3)
grid.arrange(p_noImp, p_magic, p_saver, nrow = 1)
dev.off()

