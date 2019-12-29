for fd in `ls impute`
do
rm impute/$fd/out
for i in `ls impute/$fd | grep .sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> impute/$fd/out
done
done
