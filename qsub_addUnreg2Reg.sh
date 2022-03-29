#!/bin/bash

#PBS -l walltime=200:00:00,nodes=1:ppn=4,mem=16gb
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

cmd="try; addpath('${p1}'); addpath('${p4}'); ${p2}('${p3}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"
echo $cmd
time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
