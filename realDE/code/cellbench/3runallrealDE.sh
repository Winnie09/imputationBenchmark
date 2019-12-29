for i in `cat /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/impute_method.txt`
do
qsub 3runrealDE.sh $i
done 

