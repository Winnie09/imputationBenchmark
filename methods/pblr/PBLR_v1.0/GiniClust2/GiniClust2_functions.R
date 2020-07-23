# R functions used

calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind = is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x = x[!na.ind]
  n = length(x)
  mu = mean(x)
  N = if (unbiased) n * (n - 1) else n * n
  ox = x[order(x)]
  dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
  dsum / (mu * N)
}

jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}


Mean.in.log2space=function(x) {   #input is a vetor in log2 space, the return is their mean in log2 space.
  return(log2(mean(2^(x)-1)+1))
}

calcul.pca = function(x,n){
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(ppk=ppk,pca=pca, m=m,use_genes=use_genes)
}


jaccard_dist_large_matrix <- function(m) { #use when matrix larger than 65,536 (2^16) cells
  A = tcrossprod(m)
  
  r<-35000 #65,536 was max possible before, now 70000 is- can increase r to increase this
  # approach is to split up matrix into 4, then recombine
  im = which(A[1:r,1:r] > 0,arr.ind=TRUE)
  im2 = which(A[(r+1):dim(A)[1],(r+1):dim(A)[1]] > 0,arr.ind=TRUE)
  im2[,1] = im2[,1]+r
  im2[,2] = im2[,2]+r
  im3 = which(A[1:r,(r+1):dim(A)[1]] > 0,arr.ind=TRUE)
  im3[,2]=im3[,2]+r
  im4 = which(A[(r+1):dim(A)[1],1:r] > 0,arr.ind=TRUE)
  im4[,1]=im4[,1]+r
  im_all<-rbind(im,im2,im3,im4)

  b = rowSums(m) 
  
  ## only non-zero values of common
  Aim = A[im_all]
  #print("line 4")
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im_all[,1],
    j = im_all[,2],
    x = Aim / (b[im_all[,1]] + b[im_all[,2]] - Aim),
    dims = dim(A)
  )
  
  J1<-J[1:r,1:r]
  oneMinusJ1<-1-J1
  index.cell.zero.1<-index.cell.zero[which(index.cell.zero<(r+1))]
  oneMinusJ1[index.cell.zero.1,index.cell.zero.1]=0
  
  J2<-J[(r+1):dim(A)[1],(r+1):dim(A)[1]]
  oneMinusJ2<-1-J2
  index.cell.zero.2<-index.cell.zero[which(index.cell.zero>r)]
  oneMinusJ2[index.cell.zero.2-r,index.cell.zero.2-r]=0
  
  J3<-J[1:r,(r+1):dim(A)[1]]
  oneMinusJ3<-1-J3
  oneMinusJ3[index.cell.zero.1,index.cell.zero.2-r]=0
  
  J4<-J[(r+1):dim(A)[1],1:r]
  oneMinusJ4<-1-J4
  oneMinusJ4[index.cell.zero.2-r,index.cell.zero.1]=0
  
  #combine:
  
  cell.cell.jaccard.distance<-matrix(rep(0,dim(A)[1]**2),ncol=dim(A)[1],nrow=dim(A)[1])
  cell.cell.jaccard.distance[1:r,1:r]<-as.vector(oneMinusJ1)
  cell.cell.jaccard.distance[(r+1):dim(A)[1],(r+1):dim(A)[1]]<-as.vector(oneMinusJ2)
  cell.cell.jaccard.distance[1:r,(r+1):dim(A)[1]]<-as.vector(oneMinusJ3)
  cell.cell.jaccard.distance[(r+1):dim(A)[1],1:r]<-as.vector(oneMinusJ4)
  
  return(cell.cell.jaccard.distance)
}

# #for MAST DE
# #additional functions
# Mean.in.log2space=function(x,pseudo.count) {
#   return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
# }
# 
# stat.log2=function(data.m, group.v, pseudo.count){
#   #data.m=data.used.log2
#   log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
#   log2.mean.r <- t(log2.mean.r)
#   colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
#   log2.mean.r = log2.mean.r[-1,]
#   log2.mean.r = as.data.frame(log2.mean.r)
#   log2.mean.r = unfactor(log2.mean.r)  #from varhandle
#   log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
#   log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
#   log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
#   results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
#   colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
#   rownames(results) = rownames(log2.mean.r)
#   return(results)
# }


####### function m.auc  ######
#install.packages("ROCR")

library("ROCR")
v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)
  
}  
####### function m.auc END ######

