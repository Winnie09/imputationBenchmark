bk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/cellbench_quality/plot/plot/bulk_onecell_scatter.pdf', width=13,height=8)
par(mfrow=c(3,5))
for (f in af){
    print(f)
    sc = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/',f,'/genebycell.rds'))
    g = intersect(rownames(bk), rownames(sc))
    bk = bk[g,]
    sc =sc[g,]
    bkid = grep('HCC827', colnames(bk))
    mixf = c(af[grep('mix',af)])
    if (!f %in% mixf) {
      scid = grep('HCC827', colnames(sc))
      
    } else {
      ct = sapply(colnames(sc), function(i) strsplit(i,':')[[1]][2])
      if ('0.00_0.00_1.00' %in% ct){
        scid = which(ct == '0.00_0.00_1.00')  
      } else if('0_0_9' %in% ct){
        scid = which(ct == '0_0_9')  
      } else if('0_0_90' %in% ct){
        scid = which(ct == '0_0_90')  
      }
    }
    
    for (i in scid[1:1]){
      y = log2(bk[,bkid[1]] + 1)
      x = log2(sc[,i] + 1)
      print(plot(y~x, pch=20, xlab= paste0('cell:',colnames(sc)[i]), main=f, ylab='bulk (log2+1)', cex=.1, col=rgb(0,0,0,alpha=0.5)))
    }
}
dev.off()
