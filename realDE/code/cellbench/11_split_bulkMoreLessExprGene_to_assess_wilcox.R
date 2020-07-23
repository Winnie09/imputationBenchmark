allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/wilcox/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/wilcox/',mtd,'/res/'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/sc_10x_5cl.rds')
ct = sub('.*:','',colnames(raw))
bkcnt = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
bkct = sub('_.*','',colnames(bkcnt))
intid = intersect(rownames(raw),rownames(bkcnt))
raw = raw[intid,]
bkcnt = bkcnt[intid,]
bkexpr = sweep(bkcnt,2,(colSums(bkcnt)/1e6),'/')
bkexpr = log2(bkexpr+1)
ove <- sapply(allmtd, function(mtd){
  print(mtd)
  t(sapply(allf,function(f) {
    print(f)
    if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/wilcox/',mtd,'/res/',f))){
      res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/wilcox/',mtd,'/res/',f))
      res = res[order(res[,'fdr']),]
      ct1 <- sub('_.*','',f)
      ct2 <- sub('.rds*','',sub('.*_','',f))
      tmpmat = raw[,ct%in%c(ct1,ct2)]
      fc <- rowMeans(bkexpr[,bkct == ct1]) - rowMeans(bkexpr[,bkct == ct2]) ## bulk
      fc <- abs(fc)
      
      highid = names(which(fc>quantile(fc,0.9)))## express in more cells
      lowid = names(which(fc<quantile(fc,0.1)))## less
      bfex <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/bulkdiff/')
      bf <- c(paste0(ct1,'_',ct2,'_diffgene.rds'),paste0(ct2,'_',ct1,'_diffgene.rds'))
      bf <- intersect(bfex,bf)
      gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/bulkdiff/',sub('.rds','',f),'_diffgene.rds'))
      tmpres <- res[res[,1] %in% highid,]
      tmp1 <- mean(sapply(c(1:100)*10,function(i) {
        mean(tmpres[1:i,1] %in% gs)  ## discuss
      }))
      tmpres <- res[res[,1] %in% lowid,]
      tmp2 <- mean(sapply(c(1:100)*10,function(i) {
        mean(tmpres[1:i,1] %in% gs)  ## discuss
      }),na.rm=T)
      c(tmp1,tmp2)
    } else {
      return(c(NA,NA))
    }
  }))
}) ## first 10 high, last 10 low
ove_high = t(ove[1:10,])
ove_low = t(ove[11:20,])
colnames(ove_high) <- colnames(ove_low) <- sub('.rds','', allf)
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/wilcox/'
dir.create(rdir,showWarnings = F, recursive = T)
saveRDS(ove_high, paste0(rdir,'bulk_sc_diffgene_overlaps_bulkMoreExprGene.rds'))
saveRDS(ove_low, paste0(rdir,'bulk_sc_diffgene_overlaps_bulkLessExprGene.rds'))
