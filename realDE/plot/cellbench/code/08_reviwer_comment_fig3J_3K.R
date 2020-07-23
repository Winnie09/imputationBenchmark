raw <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds')
ct <- sub('.*:','',colnames(raw))

mtd <- 'raw'
res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/mast/',mtd,'/res/',f))
res$pvalue = p.adjust(res$pvalue,method='fdr')
res = res[order(res[,'pvalue'], -abs(res[,'stat'])),]
ct1 <- sub('_.*','',f)
ct2 <- sub('.rds*','',sub('.*_','',f))
tmpmat = raw[,ct%in%c(ct1,ct2)]
fc <- rowMeans(raw[,ct == ct1]) - rowMeans(raw[,ct == ct2])
fc <- abs(fc)
highid = names(which(fc>quantile(fc,0.9)))## express in more cells
lowid = names(which(fc<quantile(fc,0.1)))## less

setwd("/scratch/users/whou10@jhu.edu/Wenpin/rna_imputation/realDE/result/cellbench/diff/mast/raw/res")
af <- list.files(getwd())
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast_noImp_p_stat_distribution.pdf', width=6, height=9)
par(mfrow=c(5,4))
for (f in af){
  tmp <- readRDS(paste0('./', f,'.rds'))[tmp[,1] %in% lowid, ]
  hist(tmp[,3], main=f, col='grey', breaks=100, xlab='Test statistics', ylab='Frequency')
  hist(tmp[,2], main=f, col='grey', breaks=100, xlab='P-value', ylab='Frequency')
}
dev.off()
  
#########
## bulk, pusedobulk
raw <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds')
bk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds')
int <- intersect(rownames(bk), rownames(raw))
raw <- raw[int, ]
bk <- bk[int, ]
g1 <- intersect(int, highid)
g2 <- intersect(int, lowid)
af <- gsub('.rds','',af)
v1 <- v2 <- u1 <- u2 <-  NULL
for (f in af){
  print(f)
  ct <- sub('.*:','',colnames(raw))
  ct1 <- sub('_.*','',f)
  ct2 <- sub('.rds*','',sub('.*_','',f))
  fc <- abs(rowMeans(raw[,ct == ct1]) - rowMeans(raw[,ct == ct2]))
  t1 <- sapply(seq(1,nrow(raw)), function(i) abs(t.test(raw[i,ct == ct1], raw[i,ct == ct2])$statistic))
  names(t1) <- rownames(raw)
  
  bk_ct <- sub('_.*','',colnames(bk))
  bk_fc <- abs(bk[, bk_ct == ct1] - bk[, bk_ct == ct2])
  
  v1[f] <- cor(fc[g1], bk_fc[g1], method = 'spearman')
  v2[f] <- cor(fc[g2], bk_fc[g2], method = 'spearman')
  u1[f] <- cor(t1[g1], bk_fc[g1], method = 'spearman')
  u2[f] <- cor(t1[g2], bk_fc[g2], method = 'spearman')
}

#######
## diff perf
diffmtd <- 'mast'
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/',diffmtd,'/')
ove_high =  readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_moreExprGene.rds'))
ove_low = readRDS(paste0(rdir,'bulk_sc_diffgene_overlaps_lessExprGene.rds'))
ove_low <- ove_low[, colnames(ove_high)]
ove_high = ove_high[complete.cases(ove_high),]
ove_low = ove_low[complete.cases(ove_low),]
colSums(ove_high)
  
######
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast_pseudobulk_bulk_LFC.pdf', width=10.5, height=5)
par(mfrow=c(1,2))
plot(colSums(ove_high) ~ v1[colnames(ove_high)], main=paste0('xy-axes correlation:',round(cor(colSums(ove_high), v1[colnames(ove_high)]),2)), ylab='Fig.3J colSums', xlab='SCC between pseudobulk and bulk abs(LFC)', pch=20)
plot(colSums(ove_low) ~ v1[colnames(ove_low)], main=paste0('xy-axes correlation:',round(cor(colSums(ove_low), v1[colnames(ove_low)]),2)),  ylab='Fig.3K colSums', xlab='SCC between pseudobulk and bulk abs(LFC)', pch=20)
dev.off()

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast_pseudobulk_stat_bulk_LFC.pdf', width=10.5, height=5)
par(mfrow=c(1,2))
plot(colSums(ove_high) ~ u1[colnames(ove_high)], main=paste0('xy-axes correlation:',round(cor(colSums(ove_high), u1[colnames(ove_high)]),2)), ylab='Fig.3J colSums', xlab='SCC between pseudobulk abs(t.statistics) and bulk abs(LFC')
plot(colSums(ove_low) ~ u1[colnames(ove_low)], main=paste0('xy-axes correlation:',round(cor(colSums(ove_low), u1[colnames(ove_low)]),2)), ylab='Fig.3K colSums', xlab='SCC between pseudobulk abs(t.statistics) and bulk abs(LFC')
dev.off()

gs <- NULL
for (f in af){
   gs[f] = length(readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/bulkdiff/', f, '_diffgene.rds')))
}

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/mast_number_goldstardardDEG.pdf', width=10.5, height=5)
par(mfrow=c(1,2))
plot(colSums(ove_high)~gs[colnames(ove_high)], main = paste0('xy-axes correlation:',round(cor(colSums(ove_high), gs[colnames(ove_high)]),2)), ylab='Fig.3J colSums', xlab='Number of gold standard DEGs', pch=20)
plot(colSums(ove_low)~gs[colnames(ove_low)], main = paste0('xy-axes correlation:',round(cor(colSums(ove_low), gs[colnames(ove_low)]),2)),  ylab='Fig.3K colSums', xlab='Number of gold standard DEGs', pch=20)
dev.off()




