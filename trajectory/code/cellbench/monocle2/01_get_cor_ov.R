library(parallel)
prep_traj_order <- function(cds,group){
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
  return(cds)
}

prep_RNA_traj_order <- function(cds,group){
  cds$H2228_to_H1975 = NA
  cds$H2228_to_H1975[cds$group=="1 0 0"] = 0
  cds$H2228_to_H1975[cds$group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_H1975[cds$group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_H1975[cds$group=="0.16 0.68 0.16"] = 3
  cds$H2228_to_H1975[cds$group=="0 1 0"] = 4
  
  cds$H2228_to_HCC827 = NA
  cds$H2228_to_HCC827[cds$group=="1 0 0"] = 0
  cds$H2228_to_HCC827[cds$group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_HCC827[cds$group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_HCC827[cds$group=="0.16 0.16 0.68"] = 3
  cds$H2228_to_HCC827[cds$group=="0 0 1"] = 4
  
  cds$H1975_to_HCC827 = NA
  cds$H1975_to_HCC827[cds$group=="0 1 0"] = 0
  cds$H1975_to_HCC827[cds$group=="0.16 0.68 0.16"] = 1
  cds$H1975_to_HCC827[cds$group=="0.33 0.33 0.33"] = 2
  cds$H1975_to_HCC827[cds$group=="0.16 0.16 0.68"] = 3
  cds$H1975_to_HCC827[cds$group=="0 0 1"] = 4
  return(cds)
}

prep_traj_wrapper <- function(cds,group){
  if ('0 0 9' %in% cds$group){
    return(prep_traj_order(cds,group))
  } else {
    return(prep_RNA_traj_order(cds,group))
  }
}

get_cds <- function(expression_matrix){
  set.seed(12345)
  group <- sub('.*:','',colnames(expression_matrix))
  colnames(expression_matrix) <- paste0(colnames(expression_matrix),'_',1:ncol(expression_matrix))
  print(str(group))
  prop <- t(sapply(group,function(i) {
    as.numeric(strsplit(sub('.*:','',i),'_')[[1]])
  }))
  # prop <- prop[,-4]
  cell_metadata <- data.frame(cell=colnames(expression_matrix),p1=prop[,1],p2=prop[,2],p3=prop[,3])
  row.names(cell_metadata) <- colnames(expression_matrix)
  gene_annotation <- data.frame(gene_short_name=row.names(expression_matrix))
  row.names(gene_annotation) <- row.names(expression_matrix)
  
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  cds <- newCellDataSet(as.matrix(expression_matrix),phenoData = pd, featureData = fd,expressionFamily=uninormal())
  
  cds$group = gsub('_',' ',group)
  cds <- prep_traj_wrapper(cds,cds$group)
  #cds = reduceDimension(cds, reduction_method = "ICA",norm_method="none",pseudo_expr=0)            #####
  #cds = orderCells(cds,num_paths=3)
  flag <- 0
  tryCatch({cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0);flag <- 1},warning=function(w){},error=function(e){})
  if (flag == 0) {
    cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0,auto_param_selection=F)
  }
  
  print(length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points))
  cds = orderCells(cds)
  tmp = table(pData(cds)$State[pData(cds)$group=='0_1_0'])  # specify H2228 as root state.
  tmp = tmp[order(tmp,decreasing = T)]
  cds = orderCells(cds,root_state=names(tmp)[1],num_paths=3) # reorder the cells
  return(cds)
}

get_max_score = function(col_anno, states){
  col_st = col_anno[col_anno$State %in% states,]
  cor1 = cor(col_st$Pseudotime, col_st$H2228_to_H1975, use = "pairwise.complete.obs",method='spearman')
  cor2 = cor(col_st$Pseudotime, col_st$H2228_to_HCC827, use = "pairwise.complete.obs",method='spearman')
  cor3 = cor(col_st$Pseudotime, col_st$H1975_to_HCC827, use = "pairwise.complete.obs",method='spearman')
  ov1 = sum(col_anno[!is.na(col_anno$H2228_to_H1975),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_H1975))
  ov2 = sum(col_anno[!is.na(col_anno$H2228_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_HCC827))
  ov3 = sum(col_anno[!is.na(col_anno$H1975_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H1975_to_HCC827))
  #sp1 = 1-sum(col_anno[is.na(col_anno$H2228_to_H1975),"State"] %in% states)/sum(is.na(col_anno$H2228_to_H1975))
  #sp2 = 1-sum(col_anno[is.na(col_anno$H2228_to_HCC827),"State"] %in% states)/sum(is.na(col_anno$H2228_to_HCC827))
  #sp3 = 1-sum(col_anno[is.na(col_anno$H1975_to_HCC827),"State"] %in% states)/sum(is.na(col_anno$H1975_to_HCC827))
  #res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3),sp=c(sp1,sp2,sp3))
  res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3))
  res_df = res_df[order(res_df$overlap,decreasing = T),]
  return(res_df[1,])
}

get_cor_ov <- function(cds){
  max_score_df = list()
  if (length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0) {
    for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
      print(i)
      rm(cds_reduced)
      tryCatch({cds_reduced <- buildBranchCellDataSet(cds, branch_point=i)},error=function(e) {}) # get branch info
      if (exists('cds_reduced')) {
        col_data_df = as.data.frame(pData(cds))
        max_score_df[[i]] = lapply(1:length(unique(pData(cds_reduced)$Branch)),function(j){
          state1 = table(pData(cds_reduced)$State[pData(cds_reduced)$Branch == unique(pData(cds_reduced)$Branch)[j]])
          state1 = state1[state1>=1] ## Aug30,19. change > to >=
          print(state1)
          max_score_df1 = get_max_score(col_data_df, names(state1))
        })
        max_score_df[[i]] = Reduce(rbind,max_score_df[[i]])      
      }
    }      
    if (length(max_score_df) ==0) {
      c(NA,NA)
    } else {
      max_score_df <- max_score_df[!sapply(max_score_df,is.null)]
      tmp = unlist(lapply(max_score_df,function(x){colMeans(x)[2]}))         ### autoimpute
      finalres = max_score_df[[which(tmp==max(tmp))]]
      return(colMeans(finalres))
    }
  } else {
    c(NA,NA)
  }
  
}

mtd = as.character(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
suppressMessages(library('monocle'))
print(mtd)
getf = list.files(paste0('./result/procimpute/cellbench/',mtd))
allf = c("cellmix1.rds","cellmix2.rds","cellmix3.rds","cellmix4.rds","RNAmix_celseq2.rds",'RNAmix_sortseq.rds')
runf = intersect(allf, getf)
res_df <- lapply(runf, function(f){
  print(f)
  expression_matrix = readRDS(paste0('./result/procimpute/cellbench/',mtd,'/',f))
  prefix = unique(sub(':.*','',colnames(expression_matrix)))  
  if (f %in% c('cellmix1.rds','cellmix2.rds',"cellmix3.rds","cellmix4.rds")){
    raw = readRDS(paste0('./data/processed/cellbench/',sub('.rds','',f),'/genebycell.rds'))
    colnames(expression_matrix) = colnames(raw)
    # cn <- sapply(cn,function(i) {
    #   paste0(round(as.numeric(strsplit(i,'_')[[1]]) * 9),collapse = '_')
    # },USE.NAMES = F)
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
  } else {
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
    # cn <- sapply(cn,function(i) {
    #   j <- as.numeric(strsplit(i,'_')[[1]][1:2])
    #   paste0(c(j,round((9-round(sum(j)*9))/9,2)),collapse = '_')
    # },USE.NAMES = F)
    cn <- sapply(cn,function(i) {
      j <- as.numeric(strsplit(i,'_')[[1]])
      paste0(j,collapse='_')
    },USE.NAMES = F)
  } 
  colnames(expression_matrix) = paste0(prefix,':',cn)
  str(expression_matrix)
  tryCatch({cds <- get_cds(expression_matrix)},warning=function(w){},error=function(e) {})
  if (exists('cds')) {
    dir.create(file.path(paste0('./trajectory/result/cellbench/monocle2/cds/',mtd)), showWarnings = FALSE, recursive = T)
    saveRDS(cds,paste0('./trajectory/result/cellbench/monocle2/cds/',mtd,'/',f))
    # max(as.numeric(cds$State))
    # plot_cell_trajectory(cds)
    res <- get_cor_ov(cds)
    print(res)
    return(res)  
  } else {c(NA,NA)}
})
saveRDS(res_df, paste0('./trajectory/result/cellbench/monocle2/cor_ov/', mtd,'.rds'))
