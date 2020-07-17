library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saucie'))
getf <- sub('.rds','', list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saucie/'))
runf <- setdiff(allf, getf)
#runf <- allf
sink('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/code/procimpute/saucie_redo.txt')
res <- lapply(runf, function(f){
      d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
      sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saucie/',f,'.csv'),data.table = F)))
      if (nrow(d) != nrow(sexpr)){
            print(f)
            writeLines(f)
      } else {
            row.names(sexpr) <- row.names(d)
            colnames(sexpr) <- colnames(d)
            saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saucie/",f,'.rds'))  
      }
})
sink()

