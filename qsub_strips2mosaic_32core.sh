#!/bin/bash

#PBS -l walltime=160:00:00,nodes=1:ppn=32
#PBS -m n
#PBS -k oe
#PBS -j oe

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
echo $PBS_JOBID
echo $PBS_O_HOST
echo $PBS_NODEFILE
echo $a1

module load gdal/2.1.3
module load matlab/2019a

echo $p1
echo $p2
echo $p3
echo $p4
echo $p5
echo $p6
echo $p7

if [ -z $p7 ]; then
    cmd="try; addpath('${p1}'); addpath('${p6}'); addpath('${p6}/intersections'); ${p2}('${p3}','${p4}',${p5}); catch e; disp(getReport(e)); exit(1); end; exit(0)"
    echo $cmd
    time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
else
    cmd="try; addpath('${p1}'); addpath('${p6}'); addpath('${p6}/intersections'); ${p2}('${p3}','${p4}',${p5},'${p7}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"
    echo $cmd
    time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
fi
