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

get_cds <- function(expression_matrix){
  set.seed(12345)
  
  cell_metadata <- data.frame(cell=colnames(expression_matrix))
  row.names(cell_metadata) <- colnames(expression_matrix)
  gene_annotation <- data.frame(gene_short_name=row.names(expression_matrix))
  row.names(gene_annotation) <- row.names(expression_matrix)
  
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  cds <- newCellDataSet(as.matrix(expression_matrix),phenoData = pd, featureData = fd,expressionFamily=uninormal())
  #cds = reduceDimension(cds, reduction_method = "ICA",norm_method="none",pseudo_expr=0)            #####
  #cds = orderCells(cds,num_paths=3)
  flag <- 0
  tryCatch({cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0);flag <- 1},warning=function(w){},error=function(e){})
  
  if (flag == 0) {
    cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0,auto_param_selection=F)
  }
  cds = orderCells(cds)
  
  tmp <- as.numeric(as.character(pData(cds)$State))
  names(tmp) <- colnames(expression_matrix)
  utmp <- unique(tmp)
  
  checkroot <- sapply(utmp,function(ss) {
    cds = orderCells(cds,root_state=ss)
    length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell)
  })
  utmp <- utmp[checkroot > 0]
  tmp <- utmp[which.min(sapply(utmp,function(scn) mean(ctlevel[match(ct[match(names(tmp)[tmp==scn],ct[,1]),3],ctlevel[,1]),2],na.rm=T)))]
  cds = orderCells(cds,root_state=tmp) # reorder the cells
  return(cds)
}

scorefunc <- function(cds) {
  if (length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0) {
    sl <- NULL
    for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
      tryCatch({cds_reduced <- buildBranchCellDataSet(cds,branch_point=i)},error=function(e) {})
      df = data.frame(pData(cds_reduced),stringsAsFactors = F)
      df <- df[order(df$Pseudotime),]
      sl <- rbind(sl,rowSums(sapply(unique(df$Branch),function(ub) {
        so <- as.character(df[df[,'Branch']==ub,1])
        soct <- ct[match(so,ct[,1]),3]
        eid <- expand.grid(1:length(soct),1:length(soct))
        eid <- eid[eid[,1]<eid[,2],]
        eid <- sprintf('%s_%s',soct[eid[,1]],soct[eid[,2]])
        c(sum(eid %in% correctorder),sum(eid %in% wrongorder))
      })))
    }
    sl[which.max(sl[,1]/rowSums(sl)),]
  } else {
    NA
  }
}

mtd = as.character(commandArgs(trailingOnly = T)[[1]])
suppressMessages(library(monocle))
print(mtd)
rdir1 = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/monocle2/cds/',mtd,'/')
dir.create(rdir1, showWarnings = FALSE, recursive = T)
rdir2 = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/monocle2/cor_ov/'
dir.create(rdir2, showWarnings = FALSE, recursive = T)
f = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd))
expression_matrix = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/',f))
str(expression_matrix)
cds <- get_cds(expression_matrix)
saveRDS(cds,paste0(rdir1,f))

score <- scorefunc(cds)
saveRDS(score, paste0(rdir2, mtd,'.rds'))
