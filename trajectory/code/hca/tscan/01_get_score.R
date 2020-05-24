ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
ctlevel <- data.frame(ct=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
row.names(ctlevel) <- ctlevel[,1]

correctorder <- wrongorder <- NULL
for(pid in c('immunepath','monopath','erypath')) {
  evct <- ctlevel[ctlevel[,pid]==1,1]
  pair <- expand.grid(evct,evct)
  pair[,1] <- as.character(pair[,1])
  pair[,2] <- as.character(pair[,2])
  pair <- pair[pair[,1]!=pair[,2],]
  corid <- which(ctlevel[pair[,1],'level'] < ctlevel[pair[,2],'level'])
  wroid <- which(ctlevel[pair[,1],'level'] > ctlevel[pair[,2],'level'])
  correctorder <- c(correctorder,sapply(corid,function(si) paste0(pair[si,],collapse = '_')))
  wrongorder <- c(wrongorder,sapply(wroid,function(si) paste0(pair[si,],collapse = '_')))
}
correctorder <- unique(correctorder)
wrongorder <- unique(wrongorder)

get_cds <- function(expression_matrix,latent=F){
  set.seed(12345)
  if (latent) {
    pr <- t(expression_matrix)
  } else {
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
  }
  for (clun in 4:nrow(pr)) {
    clu <- kmeans(pr,clun)$cluster
    tmp <- which.min(sapply(1:clun,function(scn) mean(ctlevel[match(ct[match(names(clu)[clu==scn],ct[,1]),3],ctlevel[,1]),2],na.rm=T)))
    
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
  list(ord=ord,clu=clu,mcl=mcl)
}

scorefunc <- function(ord) {
  rowSums(sapply(ord,function(so) {
    so <- sub('\\.','-',so)
    soct <- ct[match(so,ct[,1]),3]
    eid <- expand.grid(1:length(soct),1:length(soct))
    eid <- eid[eid[,1]<eid[,2],]
    eid <- sapply(1:nrow(eid),function(i) paste0(soct[eid[i,1]],'_',soct[eid[i,2]]))
    c(sum(eid %in% correctorder),sum(eid %in% wrongorder))
  }))
}

mtd = as.character(commandArgs(trailingOnly = T)[[1]])
library(TSCAN)
library(igraph)
print(mtd)
rdir1 = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/tscan/cds/',mtd,'/')
dir.create(rdir1, showWarnings = FALSE, recursive = T)
rdir2 = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/tscan/cor_ov/'
dir.create(rdir2, showWarnings = FALSE, recursive = T)
f = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd))
expression_matrix = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/',f))

str(expression_matrix)
if (grepl('latent',mtd)){
  cds <- get_cds(expression_matrix,latent=T)
} else {
  cds <- get_cds(expression_matrix,latent=F)
}
saveRDS(cds,paste0(rdir1,f))

score <- scorefunc(cds$ord)
saveRDS(score, paste0(rdir2, mtd,'.rds'))












