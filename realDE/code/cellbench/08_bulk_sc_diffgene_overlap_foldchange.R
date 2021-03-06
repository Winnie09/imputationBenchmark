allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/foldchange/')
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/foldchange/',allmtd[1],'/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/foldchange/',mtd,'/',f))){
        res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/diff/foldchange/',mtd,'/',f))
        res = res[order(res[,'fc'],decreasing=T),]
        gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/bulkdiff/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(res[1:i,1] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})
ove = t(ove)
colnames(ove) = sub('.rds','', colnames(ove))
saveRDS(ove, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/cellbench/overlap/foldchange/bulk_sc_diffgene_overlaps.rds')
