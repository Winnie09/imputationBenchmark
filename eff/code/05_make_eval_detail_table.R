setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')

t = readRDS('./rna_imputation/eff/result/time.rds')
rownames(t) = t[,1]
t = t[,-1]
colnames(t) = paste0('time_',colnames(t),'(min)')
t_scaled = readRDS('./rna_imputation/eff/result/time_scaled.rds')
colnames(t_scaled) = paste0('time_',colnames(t_scaled),'_scaled')
time = readRDS('./rna_imputation/result/perf/assess/time.rds')

m = readRDS('./rna_imputation/eff/result/memory.rds')
rownames(m) = m[,1]
m = m[,-1]
dn = dimnames(m)
m = matrix(as.numeric(sub('K','',m)),ncol=4)
dimnames(m) = dn
m = m/(1024^2)
colnames(m) = paste0('memory_',colnames(m),'(GB)')

m_scaled = readRDS('./rna_imputation/eff/result/memory_scaled.rds')
colnames(m_scaled) = paste0('memory_',colnames(m_scaled),'_scaled')
memory = readRDS('./rna_imputation/result/perf/assess/memory.rds')
scalability = readRDS('./rna_imputation/result/perf/assess/scalability.rds')

o = rownames(t)
df = cbind(t, t_scaled[o,], time=time[o], m[o,], m_scaled[o,], memory=memory[o], scalability=scalability[o])

tb = readRDS('./rna_imputation/resource/method_color.rds')
rownames(df) = tb[match(rownames(df), tb[,'shortName']),'fullName'] 
write.csv(df,'./rna_imputation/eff/result/efficiency_detail_scores.csv')



