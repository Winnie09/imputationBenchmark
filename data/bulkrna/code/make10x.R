library(data.table)

allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/10xtsv/')
meta = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/10xmeta.csv')
f = 'GSE129240_rsem_expected_counts.tsv'
tmp <- fread(paste0('/scratch/users/whou10@jhu.edu/Wenpin/rna_imputation/data/bulkrna/10xtsv/',f),data.table = F)
rn = tmp[,1]
tmp = as.matrix(tmp[,-1])
rownames(tmp) = sub('\\..*','',rn)
libsize = colSums(tmp)/1e6
a = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
a[,1] <- as.character(a[,1])
a[,2] <- as.character(a[,2])
a[,1] = sub('\\..*','',a[,1])
tmp <- tmp[row.names(tmp) %in% a[,1],]
rownames(tmp) = a[match(rownames(tmp),a[,'geneid']),'genename']
library(Matrix)
process10x_rmDupGenes <- function(genebycellmat){
  tb = genebycellmat
  tb <- tb[rowSums(tb) > 0,]
  gn = rownames(tb)  
  rs <- rowSums(tb)
  kid <- sapply(unique(gn),function(sid) {
    tmp <- which(gn==sid)
    if (length(tmp)==1) {
      tmp
    } else {
      tmp[which.max(rs[tmp])]
    }
  })
  tb <- tb[kid,]
  row.names(tb) <- gn[kid]
  tb <- tb[!grepl('^MT-',row.names(tb)),]
  tb = round(tb)
}

cnt = process10x_rmDupGenes(tmp)
# mat <- mat[,c(1,2,13,14)]
 mat = sweep(cnt, 2, libsize,'/')
mat = log2(mat+1)
mat = mat[,c(1,2,13,14)]
colnames(mat) = c('239T','239T','Jurkat','Jurkat')
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/jurkat_hex.rds')

library(limma)

design <- cbind(1,c(1,1,0,0))
fit <- lmFit(voom(mat,design),design=design)
fit <- eBayes(fit)
res <- topTable(fit,n=nrow(mat),coef=2)
diffgene <- row.names(res)[res[,5] < 0.05]

