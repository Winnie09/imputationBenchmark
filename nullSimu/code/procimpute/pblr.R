allf = sub('.csv', '',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/pblr'))
res <- sapply(allf, function(f){
            sexpr <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/pblr/',f,'.csv'),as.is=T, header=F)
            sexpr = as.matrix(sexpr)
            sexpr = log2(sexpr + 1)
            d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/',f, '/genebycell.rds'))
            row.names(sexpr) <- row.names(d)
            colnames(sexpr) <- colnames(d)
            saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/pblr/',f,'.rds'))    
})

