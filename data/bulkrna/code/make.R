# download ENCODE metadata, and then pair mca data with its bulk.
tb = read.csv('/home-4/whou10@jhu.edu/data/ENCODE/mouse/bulkrna/meta/metadata.tsv',sep='\t',as.is=T)
tb = tb[tb$Output.type == 'gene quantifications' & tb$Biosample.treatments == "" & tb$Assembly == "mm10",]

mcadata = list.files('/home-4/whou10@jhu.edu/data/mca/raw/500more_dge/')
ct = gsub("\\..*","",gsub("_.*","",mcadata))

ct = sub("1","",ct)
ct = sub("2","",ct)
ct = sub("3","",ct)
ct = sub("4","",ct)
ct = sub("5","",ct)
ct = sub("6","",ct)
ct = tolower(ct)
ct[which(ct == 'bonemarrow')] = 'bone marrow'
ct[which(ct == 'mammarygland')] = 'mammary gland'
ct[which(ct == 'placentae')] = 'placenta'

tmp = list()
for (i in 1:length(ct)){
  print(i)
  if (sum(which(tb$Biosample.term.name == ct[i])) != 0) tmp[[i]] = data.frame(mca_filename = mcadata[i], ENCODE_celltype = tb$Biosample.term.name[which(tb$Biosample.term.name == ct[i])], ENCODE_accession = tb$File.accession[which(tb$Biosample.term.name == ct[i])])
} 
df = do.call(rbind,tmp)
saveRDS(df, '/home-4/whou10@jhu.edu/data/mca/meta/correspond_ENCODE_bulk.rds')

setwd('/home-4/whou10@jhu.edu/data/ENCODE/mouse/bulkrna/mca/')
for (i in unique(df$ENCODE_accession)){
  system(paste('wget ',tb$File.download.URL[which(tb$File.accession == i)]))
}
