#!/bin/bash

#PBS -l walltime=40:00:00,nodes=1:ppn=8,mem=48gb
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

cmd="addpath('${p1}'); addpath('${p4}'); writeTileToTif('${p2}',${p5},'${p3}'); exit"
echo $cmd
time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
