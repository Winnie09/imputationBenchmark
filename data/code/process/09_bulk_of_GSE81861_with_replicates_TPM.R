meta <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/geneexpr_geneopenness/data/May19/metadata/match/match_DNase_RNA_UW_CSHL_together_RNASample.csv')
load('/home-4/whou10@jhu.edu/scratch/Wenpin/geneexpr_geneopenness/data/May19/metadata/RNAseq/processed/processed.rda')
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/geneexpr_geneopenness/data/May19/hg19/RNAseq/matrix/TPM.rds')
meta <- meta[meta[,2] %in% c('GM12878','H1-hESC','A549','IMR-90','K562'),]
hg19_file <- hg19_file[hg19_file[,4] %in% meta[,1],]
hg19_file <- hg19_file[hg19_file[,'Library made from']=='polyadenylated mRNA',]
d <- d[,paste0(hg19_file[,'File accession'],'.tsv')]
v = hg19_file[,'Biosample term name'] 
v[v=='H1-hESC'] = 'H1'
v[v=='IMR-90'] = 'IMR90'
colnames(d) <- paste0(v,'_',1:ncol(d))
tb = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
rownames(d) = tb[match(rownames(d), tb[,'geneid']), 'genename']
d = d[!duplicated(rownames(d)),]
saveRDS(d, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/bulk_of_GSE81861_with_replicates_TPM.rds')
