allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/wilcox/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/wilcox/',mtd,'/res/'))
ove <- sapply(allmtd, function(mtd){
  sapply(allf,function(f) {
      print(f)
      if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/wilcox/',mtd,'/res/',f))){
        res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/wilcox/',mtd,'/res/',f))
        res = res[order(res[,'fdr']),]
        gs = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/bulkdiff/', sub('.rds','',f),'_diffgene.rds'))
        tmp <- mean(sapply(c(1:100)*10,function(i) {
          mean(res[1:i,1] %in% gs)  ## discuss
        }))  
      } else {
        return(NA)
      }
  })
})

