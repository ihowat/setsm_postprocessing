#!/bin/bash

#PBS -l walltime=40:00:00,nodes=1:ppn=2,mem=8gb
#PBS -m n
#PBS -k oe
#PBS -j oe

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
echo $PBS_JOBID
echo $PBS_O_HOST
echo $PBS_NODEFILE
echo $a1

module load gdal/2.1.1
module load matlab/2019a

echo $p1
echo $p2
echo $p3
echo $p4

tiles="${p3//;/','}"

cmd="addpath('${p1}'); addpath('${p4}'); batch_batchMergeQuadTileBuffer('${p2}',{'${tiles}'}); exit"
echo $cmd
time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
