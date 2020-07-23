library(matrixStats)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/')
library(parallel)
res <- mclapply(allmtd, function(mtd){
      set.seed(12345)
      print(mtd)
      allf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/', mtd)))
      mdata <- sdata <- tdata <- mdata2 <- sdata2 <- tdata2  <-NULL
      for (f in allf){
            print(f)
            mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/', mtd,'/',f,'.rds'))
            gn =  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f,'/diffgn.rds'))
            cnum1 = as.numeric(strsplit(f,'_')[[1]][1])
            cnum2 = as.numeric(strsplit(f,'_')[[1]][2])
            p = as.numeric(strsplit(f,'_')[[1]][3])
            q = as.numeric(strsplit(f,'_')[[1]][4])
            ## diff genes
            mat1 = mat[intersect(row.names(mat),gn),1:cnum1]
            m1 = rowMeans(mat1)
            mat2 = mat[intersect(row.names(mat),gn),(cnum1+1):ncol(mat)]
            m2 = rowMeans(mat2)
            m = sample(abs(m1 - m2), min(length(m1), 200))
            
            s1 <- apply(mat1,1,sd)
            s1 = rowSds(mat1)
            s2 = rowSds(mat2)
            s = sample(sqrt((s1^2/cnum1) +  (s2^2/cnum2)), min(length(s1),200))
            
            t = m/s
            mdata = rbind(mdata, data.frame(method = mtd, m = m))    
            sdata = rbind(sdata, data.frame(method = mtd, s = s))
            tdata = rbind(tdata, data.frame(method = mtd, t = t))
            
            ###### non-diff genes
            mat21 = mat[setdiff(row.names(mat),gn),1:cnum1]
            m21 = rowMeans(mat21)
            mat22 = mat[setdiff(row.names(mat),gn),(cnum1+1):ncol(mat)]
            m22 = rowMeans(mat22)
            m2 = sample(abs(m21 - m22), 2000)
            
            s21 = rowSds(mat21)
            s22 = rowSds(mat22)
            s2 = sample(sqrt((s21^2/cnum1) +  (s22^2/cnum2)), 2000)
            
            t2 = m2/s2
            mdata2 = rbind(mdata2, data.frame(method = mtd, m = m2))    
            sdata2 = rbind(sdata2, data.frame(method = mtd, s = s2))
            tdata2 = rbind(tdata2, data.frame(method = mtd, t = t2))
      }
      saveRDS(list(mdata = mdata, tdata = tdata, sdata = sdata, mdata2 = mdata2, sdata2 = sdata2, tdata2=tdata2),
              paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/gene_t_statistic_distribution_plotdata/',mtd,'.rds'))
      return(0)
}, mc.cores = 12)

