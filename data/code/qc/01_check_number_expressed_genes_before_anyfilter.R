num.expressed.gene = list()
prop.mito.read = list()
## GSE
tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/GSE81861/GSE81861_Cell_Line_COUNT.csv',as.is = T)
rownames(tb) = tb[,1]
tb = as.matrix(tb[,-1])
tb = round(tb)
summary(colSums(tb>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1074    8784   13218   11885   15383   18728 
num.expressed.gene[['GSE81861']] =  colSums(tb>0)
g = sapply(rownames(tb),function(i)strsplit(i,'_')[[1]][2])
prop.mito.read[['GSE81861']] = colSums(tb[(grepl('^MT',g)), ])/colSums(tb)


## 10x
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

tb1 = as.matrix(readMM('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/matrix.mtx'))
g1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/genes.tsv')
barc1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/barcodes.tsv')
barc1 = paste0('293T_', barc1[,1])
rownames(tb1) = g1[,2]
colnames(tb1) = barc1

tb2 = as.matrix(readMM('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/matrix.mtx'))
g2 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/genes.tsv')
barc2 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/barcodes.tsv')
barc2 = paste0('jurkat_',barc2[,1])
rownames(tb2) = g2[,2]
colnames(tb2) = barc2
## cbind
allg = union(rownames(tb1),rownames(tb2))
mat = matrix(0,nrow=length(allg),ncol=sum(ncol(tb1),ncol(tb2)))
dimnames(mat) = list(allg, c(colnames(tb1),colnames(tb2)))
mat[rownames(tb1),colnames(tb1)] = tb1
mat[rownames(tb2),colnames(tb2)] = tb2
print(summary(colSums(mat>0)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  524    2927    3277    3293    3644    6088 
num.expressed.gene[['10xcellline']] <- colSums(mat>0)
prop.mito.read[['10xcellline']] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)


## cellmix
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/'
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/'
allf = c('cellmix1.count.csv','cellmix2.count.csv','cellmix3.count.csv','cellmix4.count.csv','cellmix5.count.csv')
for (f in allf){
  print(f)
  savef = sub('.count.csv','',f)
  meta = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.metadata.csv'),as.is=T))
  cnt = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.count.csv'), as.is=T))
  m =  meta[match(colnames(cnt),row.names(meta)),c('H1975', 'H2228', 'HCC827')]
  m =matrix(as.numeric(m),nrow=nrow(m))
  # colnames(cnt) = paste0('H1975_H2228_HCC827:', apply(t(round(apply(m,1,'/',9),2)),1,paste,collapse='_'))
  colnames(cnt) = paste0('H1975_H2228_HCC827:', apply(m,1,paste,collapse='_'))
  colnames(cnt) = paste0(colnames(cnt),':',1:length(colnames(cnt)))
  mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
  mtch = as.matrix(mtch)
  row.names(cnt) = mtch[match(row.names(cnt), sub('\\..*','',mtch[,'geneid'])),'genename']
  mat = cnt[!is.na(row.names(cnt)),]
  print(summary(colSums(mat>0)))
  num.expressed.gene[[sub('.csv','',f)]] <- colSums(mat>0)
  prop.mito.read[[sub('.csv','',f)]] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)
}

# [1] "cellmix1.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     937    1583    2100    2198    2612    6028 
# [1] "cellmix2.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1661    2636    3246    3396    3995    8197 
# [1] "cellmix3.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1583    2628    3330    3492    4148    8207 
# [1] "cellmix4.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2355    3971    4855    4955    5827   11597 
# [1] "cellmix5.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   15431   16788   17596   17432   18319   18953 

## rnamix
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/'
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/'
allf = c('RNAmix_sortseq.count.csv','RNAmix_celseq2.count.csv')
for (f in allf){
  print(f)
  savef = sub('.count.csv','',f)
  meta = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.metadata.csv'),as.is=T))
  cnt = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.count.csv'), as.is=T))
  colnames(cnt) = paste0('H1975_H2228_HCC827:', apply(meta[match(colnames(cnt),row.names(meta)),c( 'H1975_prop','H2228_prop', 'HCC827_prop')],1,paste,collapse='_'))
  colnames(cnt) = paste0(colnames(cnt),':',1:length(colnames(cnt)))  
  mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
  mtch = as.matrix(mtch)
  row.names(cnt) = mtch[match(row.names(cnt), sub('\\..*','',mtch[,'geneid'])),'genename']
  mat = cnt[!is.na(row.names(cnt)),]
  print(summary(colSums(mat>0)))
  num.expressed.gene[[sub('.csv','',f)]] <- colSums(mat>0)
  prop.mito.read[[sub('.csv','',f)]] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)
  
}

# [1] "RNAmix_sortseq.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2519    4878    6692    6676    8361   11346 
# [1] "RNAmix_celseq2.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    3212    5599    6940    7003    8429   10910 




## 3cl
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/'
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/'
allf = c('sc_celseq2.count.csv','sc_dropseq.count.csv','sc_10x.count.csv')

for (f in allf){
  print(f)
  savef = sub('.count.csv','',f)
  meta = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.metadata.csv'),as.is=T))
  cnt = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.count.csv'), as.is=T))
  colnames(cnt) = paste0(colnames(cnt),':',meta[match(colnames(cnt),row.names(meta)),'cell_line_demuxlet'])
  mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
  mtch = as.matrix(mtch)
  row.names(cnt) = mtch[match(row.names(cnt), sub('\\..*','',mtch[,'geneid'])),'genename']
  mat = cnt[!is.na(row.names(cnt)),]
  print(summary(colSums(mat>0)))
  num.expressed.gene[[sub('.csv','',f)]] <- colSums(mat>0)
  prop.mito.read[[sub('.csv','',f)]] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)
  
}
#   [1] "sc_celseq2.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2795    5782    7286    7017    8362   12668 
# [1] "sc_dropseq.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    3092    4699    5483    5678    6431    8897 
# [1] "sc_10x.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    6577    8439    8951    8930    9434   11127

## 5cl
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/'
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/'
allf = c('sc_celseq2_5cl_p1.count.csv','sc_celseq2_5cl_p2.count.csv','sc_celseq2_5cl_p3.count.csv','sc_10x_5cl.count.csv')
for (f in allf){
  print(f)
  savef = sub('.count.csv','',f)
  meta = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.metadata.csv'),as.is=T))
  cnt = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.count.csv'), as.is=T))
  colnames(cnt) = paste0(colnames(cnt),':',meta[match(colnames(cnt),row.names(meta)),'cell_line_demuxlet'])
  if (f != "sc_10x_5cl.count.csv"){
    mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
    mtch = as.matrix(mtch)
    row.names(cnt) = mtch[match(row.names(cnt), sub('\\..*','',mtch[,'geneid'])),'genename']
    mat = cnt[!is.na(row.names(cnt)),]
  } else {
    mat = cnt
  }
  print(summary(colSums(mat>0)))
  num.expressed.gene[[sub('.csv','',f)]] <- colSums(mat>0)
  prop.mito.read[[sub('.csv','',f)]] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)
}
# [1] "sc_celseq2_5cl_p1.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2938    5015    5724    5978    6881   11183 
# [1] "sc_celseq2_5cl_p2.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2926    4741    5407    5556    6378   10289 
# [1] "sc_celseq2_5cl_p3.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     800    3918    4597    4713    5405    8734 
# [1] "sc_10x_5cl.count.csv"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1328    3649    4437    4360    5135    7549 

## hca
filelist <- paste0(list.files('/home-4/zji4@jhu.edu/scratch/scdata/HCA/data/immune_census/align/',pattern = '^MantonBM6',full.names = T),'/outs/filtered_gene_bc_matrices/hg19')
smat <- sapply(1:length(filelist),function(id) {
  d <- filelist[id]
  genefile <- paste0(d,"/genes.tsv")
  barcodefile <- paste0(d,"/barcodes.tsv")
  matrixfile <- paste0(d,"/matrix.mtx")
  gn <- read.table(genefile,as.is=T)
  cn <- readLines(barcodefile)
  expr <- as.matrix(readMM(matrixfile))
  row.names(expr) <- paste0(gn[,1],':',gn[,2])
  colnames(expr) <- cn
  expr
})

gene <- unique(unlist(lapply(smat,row.names)))
cell <- unique(unlist(lapply(smat,colnames)))
mat <- matrix(0,nrow=length(gene),ncol=length(cell),dimnames = list(gene,cell))
for (i in smat) {
  mat[row.names(i),colnames(i)] <- i
}
print(summary(colSums(mat>0)))
num.expressed.gene[['hca']] <- colSums(mat>0)
prop.mito.read[['hca']] <- colSums(mat[(grepl('^MT',rownames(mat))), ])/colSums(mat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  365    1093    1480    2003    2482    7208 


saveRDS(num.expressed.gene,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/qc/num.expressed.gene.rds')
saveRDS(prop.mito.read,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/qc/prop.mito.read.rds')
