tb = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.rds')
cl = unique(sub('_.*','',colnames(tb)))
cl[which(cl=='IMR90')] = 'IMR-90'
cl[which(cl=='H1')] = 'H1-hESC'
load('/home-4/whou10@jhu.edu/scratch/Wenpin/geneexpr_geneopenness/data/metadata/RNAseq/raw/rawmatrix_meta.rda')
names(hg19_file) = gsub(' ','.',names(hg19_file))
intcl = intersect(hg19_file$Biosample.term.name,cl)

load('/home-4/whou10@jhu.edu/scratch/Wenpin/geneexpr_geneopenness/data/raw/hg19/RNAseq/matrix/FPKM/rawmatrix.rda')
colnames(FPKM) = sub('.tsv','',colnames(FPKM))
FPKM = log10(FPKM + 1) # [1] 56082   766
numexprgene = colSums(FPKM > 0)
FPKM = FPKM[,numexprgene > 10000] # [1] 56082   763
expr = FPKM[,hg19_file[hg19_file$Biosample.term.name %in% intcl,'File.accession']]
encff = hg19_file[hg19_file$Biosample.term.name %in% intcl,'File.accession']
encsr = hg19_file[match(encff,hg19_file$File.accession), 'Experiment.accession']
library(Matrix.utils)
expr = t(as.matrix(aggregate.Matrix(t(expr),encsr,fun = 'mean')))
gid_name = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
expr = expr[rownames(expr)%in%gid_name$geneid,]
rownames(expr) = gid_name[match(rownames(expr), gid_name[,'geneid']),'genename']

res = sapply(unique(rownames(expr)),function(gn){
  id = which(rownames(expr)==gn)
  if (length(id)>1){
    id[which.max(rowSums(expr[id,]))]
  } else {
    id
  }
})
expr <- expr[res,]
colnames(expr) = paste0(colnames(expr),'_', hg19_file[match(colnames(expr),hg19_file$Experiment.accession),'Biosample.term.name'])
saveRDS(expr,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_encsr.rds')

bcl = sub('.*_','',colnames(expr))
library(Matrix.utils)
expr = t(as.matrix(aggregate.Matrix(t(expr),bcl,fun = 'mean')))
saveRDS(expr,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
