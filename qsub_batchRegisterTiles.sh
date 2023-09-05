#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=4,mem=30gb
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
#echo Node list file: $PBS_NODEFILE
#echo Nodes assigned to job: $(cat $PBS_NODEFILE)
#echo Running node index: $PBS_O_NODENUM
echo
echo Running on hostname: $HOSTNAME
echo Parent PID: $PPID
echo Process PID: $$
echo
working_dir="$PBS_O_WORKDIR"

echo Number of physical cores: $(grep "^core id" /proc/cpuinfo | sort -u | wc -l)
echo Number of virtual cores: $(nproc --all)
mem_report_cmd="free -g"
echo "Memory report [GB] ('${mem_report_cmd}'):"
echo -------------------------------------------------------------------------
${mem_report_cmd}
echo -------------------------------------------------------------------------
echo

echo "Changing to working directory: ${working_dir}"
cd "$working_dir" || exit 1
echo

## Arguments/Options
tiledir="$ARG_TILEDIR"
resolution="$ARG_RESOLUTION"
is2dir="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile/";

set -uo pipefail

if [ -n "$tiledir" ]; then
    # Make sure this is an absolute path
    tiledir=$(readlink -f "$tiledir")
else
    tiledir="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles/"
#    tiledir="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing/"
fi
if [ -z "$resolution" ]; then
    resolution='10m';
#    resolution='2m';
fi


## Validate arguments
if [ ! -d "$tiledir" ]; then
    echo "Tiledir does not exist: ${tiledir}"
    exit 1
fi
if [ ! -d "$is2dir" ]; then
    echo "is2dir does not exist: ${is2dir}"
    exit 1
fi

matlab_cmd="try; addpath('/mnt/pgc/data/common/repos/setsm_postprocessing4'); batchRegisterTiles('${tiledir}', '${is2dir}', 'resolution','${resolution}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"
#matlab_cmd="try; addpath('/mnt/pgc/data/scratch/erik/repos/setsm_postprocessing4'); batchRegisterTiles('${tiledir}', '${is2dir}', 'resolution','${resolution}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"

echo "Argument tile directory: ${tiledir}"
echo "Matlab command: \"${matlab_cmd}\""

time matlab -nojvm -nodisplay -nosplash -r "$matlab_cmd"
