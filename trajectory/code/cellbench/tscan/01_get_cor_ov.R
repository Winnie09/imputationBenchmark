prep_traj_order <- function(group){
  group <- gsub('_',' ',group)
  cds <- list()
  cds$H2228_to_H1975 = NA
  cds$H2228_to_H1975[group=="0 9 0"] = 0
  cds$H2228_to_H1975[group=="1 7 1"] = 1
  cds$H2228_to_H1975[group=="2 5 2"] = 2
  cds$H2228_to_H1975[group=="3 3 3"] = 3
  cds$H2228_to_H1975[group=="5 2 2"] = 4
  cds$H2228_to_H1975[group=="7 1 1"] = 5
  cds$H2228_to_H1975[group=="9 0 0"] = 6
  
  cds$H2228_to_HCC827 = NA
  cds$H2228_to_HCC827[group=="0 9 0"] = 0
  cds$H2228_to_HCC827[group=="1 7 1"] = 1
  cds$H2228_to_HCC827[group=="2 5 2"] = 2
  cds$H2228_to_HCC827[group=="3 3 3"] = 3
  cds$H2228_to_HCC827[group=="2 2 5"] = 4
  cds$H2228_to_HCC827[group=="1 1 7"] = 5
  cds$H2228_to_HCC827[group=="0 0 9"] = 6
  
  cds$H1975_to_HCC827 = NA
  cds$H1975_to_HCC827[group=="9 0 0"] = 0
  cds$H1975_to_HCC827[group=="7 1 1"] = 1
  cds$H1975_to_HCC827[group=="5 2 2"] = 2
  cds$H1975_to_HCC827[group=="3 3 3"] = 3
  cds$H1975_to_HCC827[group=="2 2 5"] = 4
  cds$H1975_to_HCC827[group=="1 1 7"] = 5
  cds$H1975_to_HCC827[group=="0 0 9"] = 6
  cds$cell <- names(group)
  return(cds)
}

prep_RNA_traj_order <- function(group){
  group <- gsub('_',' ',group)
  cds <- list()
  cds$H2228_to_H1975 = NA
  cds$H2228_to_H1975[group=="1 0 0"] = 0
  cds$H2228_to_H1975[group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_H1975[group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_H1975[group=="0.16 0.68 0.16"] = 3
  cds$H2228_to_H1975[group=="0 1 0"] = 4
  
  cds$H2228_to_HCC827 = NA
  cds$H2228_to_HCC827[group=="1 0 0"] = 0
  cds$H2228_to_HCC827[group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_HCC827[group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_HCC827[group=="0.16 0.16 0.68"] = 3
  cds$H2228_to_HCC827[group=="0 0 1"] = 4
  
  cds$H1975_to_HCC827 = NA
  cds$H1975_to_HCC827[group=="0 1 0"] = 0
  cds$H1975_to_HCC827[group=="0.16 0.68 0.16"] = 1
  cds$H1975_to_HCC827[group=="0.33 0.33 0.33"] = 2
  cds$H1975_to_HCC827[group=="0.16 0.16 0.68"] = 3
  cds$H1975_to_HCC827[group=="0 0 1"] = 4
  cds$cell <- names(group)
  return(cds)
}

prep_traj_wrapper <- function(group){
  if ('0_0_9' %in% group){
    return(prep_traj_order(group))
  } else {
    return(prep_RNA_traj_order(group))
  }
}

get_cds <- function(expression_matrix){
  set.seed(12345)
  group <- sub('.*:','',colnames(expression_matrix))
  names(group) <- colnames(expression_matrix)
  prep <- prep_traj_wrapper(group)
  rsd <- apply(expression_matrix,1,sd)
  rm <- rowMeans(expression_matrix)
  cv <- rsd/rm
  cv[is.na(cv)] <- 0
  expression_matrix <- expression_matrix[cv > median(cv),]
  d <- prcomp(t(expression_matrix), scale = T)
  sdev <- d$sdev[1:20]
  x <- 1:20
  optpoint <- which.min(sapply(2:10, function(i) {
    x2 <- pmax(0, x - i)
    sum(lm(sdev ~ x + x2)$residuals^2)
  }))
  pcadim = optpoint + 1
  pr <- d$x[,1:pcadim]
  for (clun in 4:nrow(pr)) {
    print(clun)
    clu <- kmeans(pr,clun)$cluster
    if ('0_9_0' %in% group) {
      tmp <- table(clu[group=='0_9_0'])
      tmp = names(sort(tmp,decreasing = T))[1]  ## cluster with the most H2228 cells
    } else { ## add <
      tmp <- table(clu[group=='0_1_0'])
      tmp = names(sort(tmp,decreasing = T))[1]
    }  ## add >
    mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
    path <- all_simple_paths(mcl$MSTtree,tmp)
    path <- lapply(path,as.vector)
    ord <- list()
    for (i in 1:length(path)) {
      target <- path[[i]]
      over <- sapply(setdiff(1:length(path),i),function(j) {
        mean(target %in% path[[j]])
      })
      if (max(over) < 1) {
        ord[[length(ord)+1]] <- TSCANorder(mcl,MSTorder=target,orderonly = T)  
      }
    }
    if (length(ord) > 1) {break()}
  }
  list(ord=ord,clu=clu,mcl=mcl,group=prep)
}

#################################### add begin
get_cds_latent <- function(expression_matrix){
  set.seed(12345)
  group <- sub('.*:','',colnames(expression_matrix))
  names(group) <- colnames(expression_matrix)
  prep <- prep_traj_wrapper(group)
  library(TSCAN)
  library(igraph)
  pr <- t(expression_matrix)
  
  for (clun in 4:nrow(pr)) {
    clu <- kmeans(pr,clun)$cluster
    if ('0_9_0' %in% group) {
      tmp <- table(clu[group=='0_9_0'])
      tmp = names(sort(tmp,decreasing = T))[1]  
    } else { ## add <
      tmp <- table(clu[group=='0_1_0'])
      tmp = names(sort(tmp,decreasing = T))[1]
    }  ## add >
    mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
    path <- all_simple_paths(mcl$MSTtree,tmp)
    path <- lapply(path,as.vector)
    ord <- list()
    for (i in 1:length(path)) {
      target <- path[[i]]
      over <- sapply(setdiff(1:length(path),i),function(j) {
        mean(target %in% path[[j]])
      })
      if (max(over) < 1) {
        ord[[length(ord)+1]] <- TSCANorder(mcl,MSTorder=target,orderonly = T)  
      }
    }
    if (length(ord) > 1) {break()}
  }
  list(ord=ord,clu=clu,mcl=mcl,group=prep)
}
#################################### end
get_max_score = function(ord,group){
  pt <- 1:length(ord)
  
  cor1 = cor(pt, group$H2228_to_H1975[match(ord,group$cell)], use = "pairwise.complete.obs", method='spearman')
  cor2 = cor(pt, group$H2228_to_HCC827[match(ord,group$cell)], use = "pairwise.complete.obs", method='spearman')
  
  ov1 = sum(!is.na(group$H2228_to_H1975[match(ord,group$cell)]))/sum(!is.na(group$H2228_to_H1975))
  ov2 = sum(!is.na(group$H2228_to_HCC827[match(ord,group$cell)]))/sum(!is.na(group$H2228_to_HCC827))
  
  res_df = data.frame(corr=abs(c(cor1,cor2)),overlap=c(ov1,ov2))
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  return(res_df[1,])
}

get_cor_ov <- function(cds){
  max_score_df = list()
  corr = c()
  ov = c()
  if (length(cds$ord) > 1) {
    for(i in 1:length(cds$ord)){
      max_j = get_max_score(cds$ord[[i]],cds$group)
      corr = c(corr, max_j$corr)
      ov = c(ov, max_j$overlap)
    }
    
    res_df = data.frame(corr=corr,overlap=ov)
    res_df = res_df[complete.cases(res_df),]
    res_df = res_df[order(res_df$corr,decreasing = T),]
    if (nrow(res_df)>=2){
      res_df = res_df[1:2,]
      return(colMeans(res_df))  
    } else {
      res_df
    }
  } else {
    c(NA,NA)
  }
}

mtd = as.character(commandArgs(trailingOnly = T)[[1]])
library(TSCAN)
library(igraph)
print(mtd)
rdir1 = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/cellbench/tscan/cds/',mtd,'/')
dir.create(rdir1, showWarnings = FALSE, recursive = T)
rdir2 = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/cellbench/tscan/cor_ov/'
dir.create(rdir2, showWarnings = FALSE, recursive = T)
getf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd))
allf = c("cellmix1.rds","cellmix2.rds","cellmix3.rds","cellmix4.rds","RNAmix_celseq2.rds",'RNAmix_sortseq.rds')
runf = intersect(allf, getf)
res_df <- lapply(runf, function(f){
  print(f)
  expression_matrix = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f))
  prefix = unique(sub(':.*','',colnames(expression_matrix)))  
  if (f %in% c('cellmix1.rds','cellmix2.rds',"cellmix3.rds","cellmix4.rds")){
    raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/',sub('.rds','',f),'/genebycell.rds'))
    colnames(expression_matrix) = colnames(raw)
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
  } else {
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
    cn <- sapply(cn,function(i) {
      j <- as.numeric(strsplit(i,'_')[[1]])
      paste0(j,collapse='_')
    },USE.NAMES = F)
  }
  colnames(expression_matrix) = paste0(prefix,':',cn)
  str(expression_matrix)
  if (grepl('latent',mtd)){
    cds <- get_cds_latent(expression_matrix)
  } else {
    cds <- get_cds(expression_matrix)  
  }
  saveRDS(cds,paste0(rdir1,f))
  res <- get_cor_ov(cds)
  print(res)
  return(res)  
})
saveRDS(res_df, paste0(rdir2, mtd,'.rds'))


