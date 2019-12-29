cd impute
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/code/impute/`
do
cd /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/code/impute/$i/
sh `ls /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/code/impute/$i | grep 'runall' `   
cd ..
done
