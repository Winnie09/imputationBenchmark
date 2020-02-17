## input: d(matrix), k(number of clusters)
## output: clu(vector)
suppressMessages(library(scran))
suppressMessages(library(igraph))
mtd = as.character(commandArgs(trailingOnly = T)[[1]])
f = as.character(commandArgs(trailingOnly = T)[[2]])
print(mtd)
print(f)
allf = c("cellmix1","cellmix2","cellmix3","cellmix4","sc_10x","sc_celseq2","sc_dropseq","sc_10x_5cl", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2","sc_celseq2_5cl_p3")
df = data.frame(allf = allf, k = c(rep(3,7), rep(5,4)))
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f,'.rds'))
d = t(d) ## cell by gene
k = as.numeric(df[which(allf==f),'k'])
get_graph <- function(rowData, num.clu){
  if (nrow(rowData)>3e3){
    start <- 50
  } else if (nrow(rowData)>2e3){
    start <- 20
  } else {
    start <- 1
  }
  for (num.nn in start:100){
    print(num.nn)
    graph = buildSNNGraph(rowData, transposed=T,k=num.nn)
    res = cluster_louvain(graph)$membership
    print(length(table(res)))
    if (length(table(res)) == num.clu){
      return(graph)
      break
    } else if (length(table(res)) < num.clu){
      break
    }
  }
}
graph <- get_graph(d,k)
if (is.null(graph)){
  hclu = hclust(dist(d))
  clu = cutree(hclu,k=k)
} else {
  res = cluster_louvain(graph)$membership
  cc <- aggregate(d,list(res),mean)
  cc <- as.matrix(cc[,-1])
  hclu <- hclust(dist(cc))
  clu <- cutree(hclu,k)
  clu <- clu[res]
  names(clu) = row.names(d)
}
saveRDS(clu,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/',mtd,'/',f,'.rds'))    


