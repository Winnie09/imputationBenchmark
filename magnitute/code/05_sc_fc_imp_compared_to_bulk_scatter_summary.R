set.seed(1234)
setwd('/dcl02/hongkai/data/whou/')
library(data.table)
source('./resource/function.R')

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
bulk <- bulk[intgene,]
magic = sexpr[intgene,]

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/saver/sc_10x_5cl.rds'))
cl = sub('.*:','',colnames(sexpr))
saver = sexpr[intgene,]

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds'))
cl = sub('.*:','',colnames(sexpr))
raw = sexpr[intgene,]


ct <- sub('.*:','',colnames(sexpr))
uct <- unique(ct)
comp <- expand.grid(1:length(uct),1:length(uct))
comp <- comp[comp[,1] < comp[,2],]
comp <- data.frame(uct[comp[,1]],uct[comp[,2]],stringsAsFactors = F)
xbfc <- xmfc <- xsfc <- xrfc <- NULL
for (i in 1:nrow(comp)) {
  xbfc <- cbind(xbfc,bulk[,comp[i,1]]-bulk[,comp[i,2]])
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  xmfc <- cbind(xmfc,magic[,c1]-magic[,c2])
  xsfc <- cbind(xsfc,saver[,c1]-saver[,c2])
  xrfc <- cbind(xrfc,raw[,c1]-raw[,c2])
}



# --------
# fluidigm
# --------

bulk = readRDS('./rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
colnames(bulk)[which(colnames(bulk)=='H1-hESC')] = 'H1'
colnames(bulk)[which(colnames(bulk)=='IMR-90')] = 'IMR90'

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds'))
int = intersect(rownames(bulk), rownames(sexpr))
bulk = bulk[int,]
saver <- sexpr[int,]

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/magic/GSE81861_Cell_Line_COUNT.rds'))
magic <- sexpr[int,]

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds'))
raw <- sexpr[int,]


ct <- sub('_.*','',colnames(sexpr))
uct <- unique(ct)
comp <- expand.grid(1:length(uct),1:length(uct))
comp <- comp[comp[,1] < comp[,2],]
comp <- data.frame(uct[comp[,1]],uct[comp[,2]],stringsAsFactors = F)
bfc <- mfc <- sfc <- rfc <- NULL
for (i in 1:nrow(comp)) {
  bfc <- cbind(bfc,bulk[,comp[i,1]]-bulk[,comp[i,2]])
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- cbind(mfc,magic[,c1]-magic[,c2])
  sfc <- cbind(sfc,saver[,c1]-saver[,c2])
  rfc <- cbind(rfc,raw[,c1]-raw[,c2])
}

summary(apply(abs(xbfc),2,median))
summary(apply(abs(bfc),2,median))

summary(apply(abs(xrfc),2,median))
summary(apply(abs(rfc),2,median))

summary(apply(abs(xmfc),2,median))
summary(apply(abs(mfc),2,median))

summary(apply(abs(xsfc),2,median))
summary(apply(abs(sfc),2,median))


summary(sapply(1:ncol(xbfc),function(i) cor(xbfc[,i],xsfc[,i])))
summary(sapply(1:ncol(bfc),function(i) cor(bfc[,i],sfc[,i])))

intxf <- intersect(row.names(bfc),row.names(xbfc))
summary(sapply(1:ncol(xbfc),function(i) cor(xbfc[intxf,i],xsfc[intxf,i])))
summary(sapply(1:ncol(bfc),function(i) cor(bfc[intxf,i],sfc[intxf,i])))


g2 <- names(sort(rowMeans(abs(bfc)),decreasing = T))[1:8000]
summary(sapply(1:ncol(bfc),function(i) cor(bfc[g2,i],sfc[g2,i])))




