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


mtd = 'magic'
f = 'cellmix2.rds'
suppressMessages(library('monocle'))
## load expr
expression_matrix = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f))
id <- which(!duplicated(t(expression_matrix)))
expression_matrix <- expression_matrix[,id]

prefix = unique(sub(':.*','',colnames(expression_matrix)))
# cn <- sub('.*:','',colnames(expression_matrix))
cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
cn <- sapply(cn,function(i) {
      j <- as.numeric(strsplit(i,'_')[[1]][1:2])
      paste0(c(j,round((9-round(sum(j)*9))/9,2)),collapse = '_')
},USE.NAMES = F)
if (f %in% c('cellmix1.rds','cellmix2.rds',"cellmix3.rds","cellmix4.rds")){
      cn <- sapply(cn,function(i) {
            paste0(round(as.numeric(strsplit(i,'_')[[1]]) * 9),collapse = '_')
      },USE.NAMES = F)
} 
colnames(expression_matrix) = paste0(prefix,':',cn)
str(expression_matrix)

## get cds
for (seed in 1:100){
      set.seed(seed)
      print(seed)
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
      cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0, max_components = 5)    ## monocel2
      dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/code/check',mtd,'/'))
      saveRDS(cds,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/code/check',mtd,'/',seed,'_cds.rds'))
      print(length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points))
}

