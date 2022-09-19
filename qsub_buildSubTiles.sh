#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=2
#PBS -m n
#PBS -k oe
#PBS -j oe

echo "PBS SUBMISSION DIR: ${PBS_O_WORKDIR}"
echo "JOB ID: ${PBS_JOBID}"
echo "NODE: ${PBS_O_HOST}"
echo "NODE FILE: ${PBS_NODEFILE}"

# move to temp folder to reduce startup time
mkdir -p ~/matlab_temp
cd ~/matlab_temp

cwd=$(pwd)
echo "CWD: ${cwd}"
echo

module load gdal/2.1.3
module load matlab/2019a


echo "var1: ${p1}"
echo "var2: ${p2}"
echo "var3: ${p3}"
echo "var4: ${p4}"
echo

echo matlab -nojvm -nodisplay -nosplash -r "addpath('${p1}'); addpath('${p4}'); ${p2}('${p3}'); exit"
time matlab -nojvm -nodisplay -nosplash -r "addpath('${p1}'); addpath('${p4}'); ${p2}('${p3}'); exit"

echo "Done"
