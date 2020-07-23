set.seed(12345)
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
png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_umi.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[,comp[i,1]]-bulk[,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[,c1]-magic[,c2]
  sfc <- saver[,c1]-saver[,c2]
  rfc <- raw[,c1]-raw[,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.off()

umidata <- list(comp = comp,
             bulk = bulk,
             ct = ct,
             magic = magic,
             saver = saver,
             raw = raw)
             

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
png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_fluidigm.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[,comp[i,1]]-bulk[,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[,c1]-magic[,c2]
  sfc <- saver[,c1]-saver[,c2]
  rfc <- raw[,c1]-raw[,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.off()

# ------------------------------------
# use umi and fluidigm intersect genes
# ------------------------------------
str(umidata$magic)
str(umidata$bulk)
intgene <- intersect(rownames(umidata$bulk), rownames(bulk))
str(intgene)
umidata$bulk = umidata$bulk[intgene, ]
umidata$saver = umidata$saver[intgene, ]
umidata$magic = umidata$magic[intgene, ]
umidata$raw = umidata$raw[intgene, ]

png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_umi_intgene.png', width = 600, height = 2000)
par(mfrow=c(nrow(umidata$comp),3))
for (i in 1:nrow(umidata$comp)) {
  bfc <- umidata$bulk[,umidata$comp[i,1]]-umidata$bulk[,umidata$comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(umidata$ct==umidata$comp[i,1],1),1)
  c2 <- sample(which(umidata$ct==umidata$comp[i,2],1),1)
  mfc <- umidata$magic[,c1]-umidata$magic[,c2]
  sfc <- umidata$saver[,c1]-umidata$saver[,c2]
  rfc <- umidata$raw[,c1]-umidata$raw[,c2]
  plot(rfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
}
dev.off()


png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_fluidigm_intgene.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[intgene,comp[i,1]]-bulk[intgene,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[intgene,c1]-magic[intgene,c2]
  sfc <- saver[intgene,c1]-saver[intgene,c2]
  rfc <- raw[intgene,c1]-raw[intgene,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.offset.seed(12345)
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
png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_umi.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[,comp[i,1]]-bulk[,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[,c1]-magic[,c2]
  sfc <- saver[,c1]-saver[,c2]
  rfc <- raw[,c1]-raw[,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.off()

umidata <- list(comp = comp,
             bulk = bulk,
             ct = ct,
             magic = magic,
             saver = saver,
             raw = raw)
             

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
png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_fluidigm.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[,comp[i,1]]-bulk[,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[,c1]-magic[,c2]
  sfc <- saver[,c1]-saver[,c2]
  rfc <- raw[,c1]-raw[,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.off()

# ------------------------------------
# use umi and fluidigm intersect genes
# ------------------------------------
str(umidata$magic)
str(umidata$bulk)
intgene <- intersect(rownames(umidata$bulk), rownames(bulk))
str(intgene)
umidata$bulk = umidata$bulk[intgene, ]
umidata$saver = umidata$saver[intgene, ]
umidata$magic = umidata$magic[intgene, ]
umidata$raw = umidata$raw[intgene, ]

png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_umi_intgene.png', width = 600, height = 2000)
par(mfrow=c(nrow(umidata$comp),3))
for (i in 1:nrow(umidata$comp)) {
  bfc <- umidata$bulk[,umidata$comp[i,1]]-umidata$bulk[,umidata$comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(umidata$ct==umidata$comp[i,1],1),1)
  c2 <- sample(which(umidata$ct==umidata$comp[i,2],1),1)
  mfc <- umidata$magic[,c1]-umidata$magic[,c2]
  sfc <- umidata$saver[,c1]-umidata$saver[,c2]
  rfc <- umidata$raw[,c1]-umidata$raw[,c2]
  plot(rfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(umidata$comp[i,],collapse = '_'))
}
dev.off()


png('./rna_imputation/magnitute/plot/sc_fc_imp_compared_to_bulk_scatter_fluidigm_intgene.png', width = 600, height = 2000)
par(mfrow=c(nrow(comp),3))
for (i in 1:nrow(comp)) {
  bfc <- bulk[intgene,comp[i,1]]-bulk[intgene,comp[i,2]]
  set.seed(12345)
  c1 <- sample(which(ct==comp[i,1],1),1)
  c2 <- sample(which(ct==comp[i,2],1),1)
  mfc <- magic[intgene,c1]-magic[intgene,c2]
  sfc <- saver[intgene,c1]-saver[intgene,c2]
  rfc <- raw[intgene,c1]-raw[intgene,c2]
  plot(rfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(mfc~bfc,main=paste0(comp[i,],collapse = '_'))
  plot(sfc~bfc,main=paste0(comp[i,],collapse = '_'))
}
dev.off()


