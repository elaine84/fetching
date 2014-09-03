#!/bin/dash
#$ -S /bin/dash
#$ -o $JOB_ID.output
#$ -e $JOB_ID.error
#$ -cwd
#$ -pe orte 2
#$ -l h_rt=16:00:00

export PYTHONPATH=~/path/to/fetching:$PYTHONPATH

n=$NSLOTS
py=pymixture
data_size=1000000
r=1.1
l=50000
b=1
d=64
c=0.999
u=0
y=1
s=2

REV=`git rev-parse --short HEAD`
TIME=`date "+%s"`
META=py=$py-data_size=$data_size-seed=$seed-jump=$jump-n=$n-r=$r-l=$l-b=$b-d=$d-u=$u-y=$y-s=$s-c=$c-$REV-$TIME
DATE=`date +%Y%m%d`

mpirun -np $n --mca btl self,tcp fetching -r $r -l $l -b $b -d $d -u $u -y $y -s $s -c $c -p py $py data_size=$data_size  seed=$seed jump=$jump

REV=`git rev-parse --short HEAD`
TIME=`date "+%s"`

python scripts/summary.py o:$JOB_ID.output e:$JOB_ID.error
python scripts/plot.py 
