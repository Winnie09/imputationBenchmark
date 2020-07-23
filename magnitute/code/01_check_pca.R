flui <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds')
str(flui)
cv = apply(flui, 1, sd)/rowMeans(flui)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(flui[cv > 1, ]), scale. = T)
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
library(ggplot2)
pd_flui = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('_.*','', rownames(pca)))
p_flui <- ggplot() + 
  geom_point(data = pd_flui, aes(x = pc1, y = pc2, color = ct)) + 
  xlab('PC1 (31.8%)') + 
  ylab('PC2 (28.6%)') + 
  theme_classic() +
  ggtitle('Fluidigm_5cellline, no_imp')

umi <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds')
str(umi)
cv = apply(umi, 1, sd)/rowMeans(umi)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(umi[cv > 1, ]), scale. = T) # slow
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
pd_umi = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('.*:','', rownames(pca)))
p_umi <- ggplot() + 
  geom_point(data = pd_umi, aes(x = pc1, y = pc2, color = ct), size = 0.5) + 
  xlab('PC1 (14.8%)') + 
  ylab('PC2 (13.5%)') + 
  theme_classic() +
  ggtitle('10x_5cellline, no_imp')

library(gridExtra)
grid.arrange(p_flui,p_umi,nrow=1)




#############3
flui <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/magic/GSE81861_Cell_Line_COUNT.rds')
str(flui)
cv = apply(flui, 1, sd)/rowMeans(flui)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(flui[cv > 1, ]), scale. = T)
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
library(ggplot2)
pd_flui = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('_.*','', rownames(pca)))
p_flui_imp <- ggplot() + 
  geom_point(data = pd_flui, aes(x = pc1, y = pc2, color = ct)) + 
  xlab('PC1 (27.3%)') + 
  ylab('PC2 (24.8%)') + 
  theme_classic() +
  ggtitle('Fluidigm_5cellline, MAGIC')

umi <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/magic/sc_10x_5cl.rds')
cv = apply(umi, 1, sd)/rowMeans(umi)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(umi[cv > 1, ]), scale. = T)
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
pd_umi_imp = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('.*:','', rownames(pca)))
p_umi_imp <- ggplot() + 
  geom_point(data = pd_umi_imp, aes(x = pc1, y = pc2, color = ct), size = 0.5) + 
  xlab('PC1 (14.5%)') + 
  ylab('PC2 (12.5%)') + 
  theme_classic() +
  ggtitle('10x_5cellline, MAGIC')

#####
flui <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds')
str(flui)
cv = apply(flui, 1, sd)/rowMeans(flui)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(flui[cv > 1, ]), scale. = T)
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
library(ggplot2)
pd_flui2 = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('_.*','', rownames(pca)))
p_flui_imp2 <- ggplot() + 
  geom_point(data = pd_flui2, aes(x = pc1, y = pc2, color = ct)) + 
  xlab('PC1 (31.1%)') + 
  ylab('PC2 (25.6%)') + 
  theme_classic() +
  ggtitle('Fluidigm_5cellline, SAVER')

umi <- readRDS('/Users/wenpinhou/Dropbox/rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds')
cv = apply(umi, 1, sd)/rowMeans(umi)
summary(cv)
hist(cv)
set.seed(12345)
pr = prcomp(t(umi[cv > 1, ]), scale. = T)
pca <- pr$x
str(pca)
str(pr)
head(pr$sdev)
pd_umi_imp2 = data.frame(pc1 = pca[,1], pc2 = pca[,2], ct = gsub('.*:','', rownames(pca)))
p_umi_imp2 <- ggplot() + 
  geom_point(data = pd_umi_imp2, aes(x = pc1, y = pc2, color = ct), size = 0.5) + 
  xlab('PC1 (12.7%)') + 
  ylab('PC2 (10.7%)') + 
  theme_classic() +
  ggtitle('10x_5cellline, SAVER')

library(gridExtra)
pdf('/Users/wenpinhou/Dropbox/rna_imputation/magnitute/plot/noimp_magic_pca.pdf', width = 7, height = 7)
grid.arrange(p_flui, p_umi, p_flui_imp, p_umi_imp, p_flui_imp2, p_umi_imp2, nrow = 3)
dev.off()


#########
library(data.table)
tb = fread('/Users/wenpinhou/Dropbox/resource/gencode.v19.annotation.gtf',data.table = F)
tb <- tb[tb[,3]=='exon',]
gn <- sub('".*','',sub('.*gene_name "','',tb[,9]))
len <- tb[,5]-tb[,4]+1
gl <- rowsum(len,gn)
gl <- gl[,1]
gl <- log2(gl + 1)
int = intersect(rownames(umi), names(gl))
scc <- apply(umi[int,], 2, cor, gl[int])
pdf('/Users/wenpinhou/Dropbox/rna_imputation/magnitute/plot/no_imp_gene_length_gene_expr_scc.pdf', width = 6, height = 3)
par(mfrow = c(1,2))
plot(umi[int,1] ~ gl[int ], pch = 20, main = paste0('one cell, xy-axis SCC=', round(cor(umi[int,1], gl[int], method = 'spearman'), 2)), xlab = 'log2(gene_length + 1)', ylab = 'gene expression') 
hist(scc, main = 'all cells SCC, no_imp')
dev.off()


int = intersect(rownames(flui), names(gl))
scc <- apply(flui[int,], 2, cor, gl[int])
par(mfrow = c(1,2))
plot(flui[int,1] ~ gl[int ], pch = 20, main = paste0('one cell, xy-axis SCC=', round(cor(flui[int,1], gl[int], method = 'spearman'), 2)), xlab = 'log2(gene_length + 1)', ylab = 'gene expression') 
hist(scc, main = 'all cells SCC, no_imp')
dev.off()



