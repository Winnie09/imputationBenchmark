setwd('/dcl02/hongkai/data/whou/')
library(data.table)
source('./resource/function.R')
get_bulk_foldchange <- function(expr){
  expr <- expr[, sort(colnames(expr))]
  mat <- lapply(1:(ncol(expr)-1), function(i){
      tmp <- sapply(((i+1):ncol(expr)), function(j){
        abs(expr[,i] - expr[,j])
      })
      if (is.vector(tmp)) tmp = data.frame(v1 = tmp)
      colnames(tmp) <- paste0(colnames(expr)[i], '_', colnames(expr)[(i+1):ncol(expr)])
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
bulk = bulk[intgene,]
bk_umi <- get_bulk_foldchange(bulk)

# --------
# fluidigm
# --------
bulk = readRDS('./rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
colnames(bulk)[which(colnames(bulk)=='H1-hESC')] = 'H1'
colnames(bulk)[which(colnames(bulk)=='IMR-90')] = 'IMR90'

sexpr = readRDS(paste0('./rna_imputation/result/procimpute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds'))
int = intersect(rownames(bulk), rownames(sexpr))
bulk = bulk[int,]
bk_flui <- get_bulk_foldchange(bulk)


pd4 = readRDS('./rna_imputation/plot/pd/15_imp_diff_eval_4subfig.pd4.rds')
coldf = readRDS('./rna_imputation/resource/method_color.rds')


df = data.frame(Dataset = c(colnames(bk_umi), colnames(bk_flui)), BulkLFC = c(apply(bk_umi, 2, median), apply(bk_flui, 2, median)), stringsAsFactors = FALSE)

lv = levels(pd4$Dataset)

dtb <- data.frame(full = lv, short = paste0(sub('\\(.*', '', lv), '_', sub('\\(.*', '', sub('.*_', '', lv))), stringsAsFactors = FALSE)

tf <- function(x) {
  s1 <- sub('_.*','',x)
  s2 <- sub('.*_','',x)
  sapply(1:length(s1),function(i) paste0(sort(c(s1[i],s2[i])),collapse = '_'))
}

df$Dataset = dtb[match(tf(df$Dataset), tf(dtb$short)),'full']


pd = cbind(pd4, BulkLFC = df[match(pd4$Dataset, df$Dataset),'BulkLFC'])

library(ggplot2)
p <- ggplot(data = pd, aes(x = BulkLFC, y = corMedian, color = platform)) + 
  geom_point() +
  theme_classic() + 
  facet_wrap(~Method) +
  xlab("Median of all genes' absolute LFC in bulk") +
  ylab('Median single-cell LFC (same as y-axis in Figure 2H)')
ggsave('./rna_imputation/magnitute/plot/performance_bulkLFC_scatter.pdf',
       p, width = 7, heigh = 7)




