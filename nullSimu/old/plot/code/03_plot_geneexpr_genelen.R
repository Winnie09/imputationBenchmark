dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/'
m = readRDS(paste0(dir,'genebycell.rds'))

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/geneexpr_genelen.pdf',width=11,height=5)
par(mfrow=c(2,5))
ct = sub('_.*','',colnames(m))
for (sct in unique(ct)){
      expr = rowMeans(m[,which(ct==sct)])
      len = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19.rds')
      df = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
      names(len) = df[match(names(len), df$geneid), 'genename']
      g = intersect(row.names(m), names(len))
      len = len[g]
      expr = expr[g]
      cor(len,expr)
      cor(log10(len),log10(expr))
      plot(log10(expr)~log10(len),pch='.',col=rgb(0,0,0,0.1),main=paste0('fluidigm,',sct,',SCC=',round(cor(log10(expr+1e-5),log10(len),method='pearson'),2)))      
}



## 10x
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/'
m = readRDS(paste0(dir,'genebycell.rds'))
ct = sub('_.*','',colnames(m))
for (sct in unique(ct)){
      expr = rowMeans(m[,which(ct==sct)])
      len = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19.rds')
      df = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
      names(len) = df[match(names(len), df$geneid), 'genename']
      g = intersect(row.names(m), names(len))
      len = len[g]
      expr = expr[g]
      cor(len,expr)
      cor(log10(len),log10(expr))
      plot(log10(expr)~log10(len),pch='.',col=rgb(0,0,0,0.1),main=paste0('10x,SCC=',round(cor(log10(expr+1e-5),log10(len),method='pearson'),2)))
}
dev.off()
