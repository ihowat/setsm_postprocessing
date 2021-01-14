#!/bin/bash

#PBS -l walltime=200:00:00,nodes=1:ppn=16,mem=64gb
#PBS -m n
#PBS -k oe
#PBS -j oe
#PBS -q batch

echo ________________________________________________________
echo
echo PBS Job Log
echo Start time: $(date)
echo
echo Job name: $PBS_JOBNAME
echo Job ID: $PBS_JOBID
echo Submitted by user: $USER
echo User effective group ID: $(id -ng)
echo
echo Hostname of submission: $PBS_O_HOST
echo Submitted to cluster: $PBS_SERVER
echo Submitted to queue: $PBS_QUEUE
echo Requested nodes per job: $PBS_NUM_NODES
echo Requested cores per node: $PBS_NUM_PPN
echo Requested cores per job: $PBS_NP
echo Node list file: $PBS_NODEFILE
echo Nodes assigned to job: $(cat $PBS_NODEFILE)
echo Running node index: $PBS_O_NODENUM
echo
echo Running on hostname: $HOSTNAME
echo Parent PID: $PPID
echo Process PID: $$
echo
echo Working directory: $PBS_O_WORKDIR
echo ________________________________________________________
echo

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
echo "var5: ${p5}"
echo "var6: ${p6}"
echo "var7: ${p7}"
echo "var8: ${p8}"
echo "var9: ${p9}"
echo "var10: ${p10}"
echo "var11: ${p11}"
echo

echo matlab -nojvm -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); ${p3}('${p4}','${p5}','${p6}','${p7}','landTile','${p8}','refDemFile','${p9}','make2m',${p10}); catch e; disp(getReport(e)); exit(1); end; exit(0);"
time matlab -nojvm -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); ${p3}('${p4}','${p5}','${p6}','${p7}','landTile','${p8}','refDemFile','${p9}','make2m',${p10}); catch e; disp(getReport(e)); exit(1); end; exit(0);"
status=$?

echo "Done"

# create finfile if matlab command exited without error
if (( status == 0 )); then
    touch "${p11}"
fi
