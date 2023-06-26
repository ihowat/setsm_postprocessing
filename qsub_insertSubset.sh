#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=2,mem=15gb
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

cd "$PBS_O_WORKDIR"

#module load gdal/2.1.3
module load matlab/2019a

## Arguments/Options
tiledir="$ARG_TILEDIR"
resolution="$ARG_RESOLUTION"
annualRootDir="$ARG_ANNUALROOTDIR"
boxShp="$ARG_BOXSHP"
pShp="$ARG_PSHP"

set -uo pipefail

if [ -n "$tiledir" ]; then
    # Make sure this is an absolute path
    tiledir=$(readlink -f "$tiledir")
else
    tiledir=""
fi
if [ -z "$resolution" ]; then
#    resolution='10m';
    resolution='2m';
fi
if [ -z "$annualRootDir" ]; then
#    annualRootDir="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/results/output_tiles_annual"
    annualRootDir="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/results/output_tiles_annual_greenland_ian"
fi
if [ -z "$boxShp" ]; then
    boxShp="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/glacier_patching/arcticdem_v41_patches_box.shp"
fi
if [ -z "$pShp" ]; then
    pShp="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/glacier_patching/arcticdem_v41_patches_p.shp"
fi


## Validate arguments
if [ ! -d "$tiledir" ]; then
    echo "Tiledir does not exist: ${tiledir}"
    exit 1
fi

#matlab_cmd="try; addpath('/mnt/pgc/data/common/repos/setsm_postprocessing_pgc'); batch_insertSubset('${tiledir}', '${resolution}', '${annualRootDir}', '${boxShp}', '${pShp}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"
matlab_cmd="try; addpath('/mnt/pgc/data/scratch/erik/repos/setsm_postprocessing_pgc'); addpath('/mnt/pgc/data/scratch/erik/repos/setsm_postprocessing4'); batch_insertSubset('${tiledir}', '${resolution}', '${annualRootDir}', '${boxShp}', '${pShp}'); catch e; disp(getReport(e)); exit(1); end; exit(0)"

echo "Argument tile directory: ${tiledir}"
echo "Matlab command: \"${matlab_cmd}\""

time matlab -nojvm -nodisplay -nosplash -r "$matlab_cmd"
