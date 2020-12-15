#!/bin/bash

#PBS -l walltime=40:00:00,nodes=1:ppn=4,mem=24gb
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

if [ "${p5}" == "false" ]; then
    echo "Building tifs and meta files"
    cmd="addpath('${p1}'); addpath('${p4}'); writeTileToTifv4('${p2}','${p3}'); tileMetav4('${p2}'); exit"

else
    echo "Building meta files only"
    cmd="addpath('${p1}'); addpath('${p4}'); tileMetav4('${p2}'); exit"

fi

echo $cmd
time matlab -nojvm -nodisplay -nosplash -r "${cmd}"
