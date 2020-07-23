allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/')
allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell')
for (mtd in allmtd){
      tmpf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/',mtd))
      print(mtd)
      if (length(grep('csv',tmpf[1])) != 0){
            tmpf = gsub('.csv','',tmpf)
      } else if (length(grep('rds',tmpf[1])) != 0){
            tmpf = gsub('.rds','',tmpf)
      } else if (length(grep('mat',tmpf[1])) != 0){
            tmpf = gsub('.mat','',tmpf)
      } 
      tmpf = gsub('clust','', gsub('scimpute_count','', gsub('pars1','', gsub('totalCounts_by_cell','', tmpf))))
      tmpf = unique(tmpf)
      print(length(allf) - length(tmpf))
      sink(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/code/impute/',mtd,'/redof.txt'))
      writeLines(setdiff(allf,tmpf))
      sink()
}

