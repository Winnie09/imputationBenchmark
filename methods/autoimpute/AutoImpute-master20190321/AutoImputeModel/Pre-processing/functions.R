# ----------------------------------------------------------------------
# Normalization by library size
# ----------------------------------------------------------------------

normalize_by_umi_2 <-function(x, min.count=2, min.cell=3) {
  print("[!normalization] Onset...")
  mat  = x
  gene_symbols = colnames(x)
  cs <- colSums(mat>min.count)
  x_use_genes <- which(cs > min.cell)
  
  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]

  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  print("[!normalization] Complete.")
  list(m=x_norm,use_genes=gene_symbols)
}

# ----------------------------------------------------------------------
# Dispersion Genes and umi_2bsetting (Gene-selection and log transform)
# ----------------------------------------------------------------------
matrix.subset<-function(normalized_data, ngenes_keep = 1000){
  print("[!geneselection] Onset...")
  df<-get_variable_gene(normalized_data$m)
  gc()
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  system("mkdir ./tmp/",ignore.stderr = T)
  system("rm ./tmp/genes",ignore.stderr = T)
  write.csv(features, file = "./tmp/genes", quote = F,row.names = F)
  write.csv(normalized_data$use_genes[features], file = "./tmp/genes_used_all", quote = F,row.names = F)
  print("[!normalization] Complete.")
  print("[!log-transform] Onset...")
  genes = read.csv(file = "./tmp/genes")
  features = genes$x
  
  # Log transformation
  m_n_whole<-normalized_data$m[,features]
  m_filt<-Matrix(log2(m_n_whole+1),sparse = T)
  print(paste("[!info] Log normalized matrix dimensions :",dim(m_filt)[1], " x ", dim(m_filt)[2]))
  print("[!log-transform] Complete.")
  return(m_filt)
  
}

get_variable_gene<-function(m) {
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion), bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}

