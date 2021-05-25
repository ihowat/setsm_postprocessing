#!/bin/bash

## PGC settings
##PBS -l walltime=200:00:00,nodes=1:ppn=16,mem=64gb
##PBS -m n
##PBS -k oe
##PBS -j oe
##PBS -q batch

## BW settings
##PBS -l nodes=1:ppn=32:xe,gres=shifter
##PBS -l walltime=96:00:00
##PBS -v CRAY_ROOTFS=SHIFTER,UDI="ubuntu:xenial"
##PBS -e $PBS_JOBID.err
##PBS -o $PBS_JOBID.out
##PBS -m n
##PBS -q high


set -uo pipefail

set +u
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
echo Working directory: $PBS_O_WORKDIR
echo ________________________________________________________
echo
CORES_PER_NODE=$PBS_NUM_PPN

#echo ________________________________________
#echo
#echo SLURM Job Log
#echo Start time: $(date)
#echo
#echo Job name: $SLURM_JOB_NAME
#echo Job ID: $SLURM_JOBID
#echo Submitted by user: $USER
#echo User effective group ID: $(id -ng)
#echo
#echo SLURM account used: $SLURM_ACCOUNT
#echo Hostname of submission: $SLURM_SUBMIT_HOST
#echo Submitted to cluster: $SLURM_CLUSTER_NAME
#echo Submitted to node: $SLURMD_NODENAME
#echo Cores on node: $SLURM_CPUS_ON_NODE
#echo Requested cores per task: $SLURM_CPUS_PER_TASK
#echo Requested cores per job: $SLURM_NTASKS
#echo Requested walltime: $SBATCH_TIMELIMIT
##echo Nodes assigned to job: $SLURM_JOB_NODELIST
##echo Running node index: $SLURM_NODEID
#echo
#echo Running on hostname: $HOSTNAME
#echo Parent PID: $PPID
#echo Process PID: $$
#echo
#echo Working directory: $SLURM_SUBMIT_DIR
#echo ________________________________________________________
#echo
#CORES_PER_NODE=$SLURM_CPUS_ON_NODE
set -u

set +u
system="$ARG_SYSTEM"
scriptdir="$ARG_SCRIPTDIR"
libdir="$ARG_LIBDIR"
tileName="$ARG_TILENAME"
outDir="$ARG_OUTDIR"
projection="$ARG_PROJECTION"
tileDefFile="$ARG_TILEDEFFILE"
stripDatabaseFile="$ARG_STRIPDATABASEFILE"
stripsDirectory="$ARG_STRIPSDIRECTORY"
waterTileDir="$ARG_WATERTILEDIR"
refDemFile="$ARG_REFDEMFILE"
tileqcDir="$ARG_TILEQCDIR"
tileParamListFile="$ARG_TILEPARAMLISTFILE"
make2m="$ARG_MAKE2M"
finfile="$ARG_FINFILE"
logfile="$ARG_LOGFILE"
set -u

if [ -z "$tileName" ]; then
    tileName="$1"
fi
if [ -z "$tileName" ]; then
    echo "argument 'tileName' not supplied, exiting"
    exit 0
fi
outDir="${outDir/<tilename>/${tileName}}"

utm_tileprefix=$(echo "$tileName" | grep -Eo '^utm[0-9]{2}[ns]')
if [ -n "$utm_tileprefix" ]; then

    if [ -z "$projection" ]; then
        projection="$utm_tileprefix"
    fi

    hemi_abbrev="${utm_tileprefix: -1}"
    if [ "$hemi_abbrev" = 'n' ]; then
        hemisphere='North'
    elif [ "$hemi_abbrev" = 's' ]; then
        hemisphere='South'
    else
        hemisphere=''
    fi
    if [ -n "$hemisphere" ]; then
        tileDefFile="${tileDefFile/<hemisphere>/${hemisphere}}"
    fi

    refDemFile="${refDemFile/<tileprefix>/${utm_tileprefix}}"
fi
if [ -z "$projection" ]; then
    echo "argument 'projection' not supplied, exiting"
    exit 0
fi

finfile="${finfile/<tilename>/${tileName}}"
logfile="${logfile/<tilename>/${tileName}}"

matlab_cmd="\
addpath('${scriptdir}'); addpath('${libdir}'); \
parpool($CORES_PER_NODE); \
run_buildSubTiles(\
'${tileName}','${outDir}',\
'${projection}','${tileDefFile}',\
'${stripDatabaseFile}','${stripsDirectory}',\
'${waterTileDir}','${refDemFile}',\
'${tileqcDir}','${tileParamListFile}',\
${make2m})"


# System-specific settings
if [ "$system" = 'pgc' ]; then
    MATLAB_WORKING_DIR="${HOME}/matlab_working_dir"
    MATLAB_ENV="module load matlab/2019a"
    MATLAB_PROGRAM="matlab"
    MATLAB_SETTINGS="-nodisplay -nosplash"
    GDAL_ENV="module load gdal/2.1.3"
    export BWPY_PREFIX=""
    APRUN_PREFIX=""

elif [ "$system" = 'bw' ]; then
    MATLAB_WORKING_DIR="/scratch/sciteam/GS_bazu/mosaic_data/matlab_working_dir"
    MATLAB_ENV=""
    MATLAB_PROGRAM="/projects/sciteam/bazu/matlab/bin/matlab"
    MATLAB_SETTINGS="-nodisplay -nodesktop -nosplash"
    export LM_LICENSE_FILE="1711@bwlm1.ncsa.illinois.edu:1711@bwlm2.ncsa.illinois.edu"
    GDAL_ENV="module load bwpy/2.0.2"
    export BWPY_PREFIX="bwpy-environ -- "
    APRUN_PREFIX="aprun -b -N 1 -d 32"
    #export CRAY_ROOTFS=SHIFTER
    #export UDI="ubuntu:xenial"
    #echo
    #echo "CRAY_ROOTFS=$CRAY_ROOTFS"
    #echo "UDI=$UDI"
    #echo
fi


task_cmd="${MATLAB_PROGRAM} ${MATLAB_SETTINGS} -r \"${matlab_cmd}\""

if [ -n "$(env | grep '^SWIFT_WORKER_PID=')" ]; then
    echo "In a Swift job"
elif [ -n "$(env | grep '^PBS_JOBID=')" ]; then
    echo "In a PBS job"
    task_cmd="${APRUN_PREFIX} ${task_cmd}"
else
    echo "Not in a Swift or PBS job"
fi

# Load environments
if [ -n "$MATLAB_ENV" ]; then
    eval "$MATLAB_ENV"
fi
if [ -n "$GDAL_ENV" ]; then
    eval "$GDAL_ENV"
fi


# Move to empty folder to reduce startup time
if [ ! -d "$MATLAB_WORKING_DIR" ]; then
    mkdir -p "$MATLAB_WORKING_DIR"
fi
cd "$MATLAB_WORKING_DIR"
cd_status=$?
if (( cd_status != 0 )); then
    echo "Failed to change to Matlab working dir: ${MATLAB_WORKING_DIR}"
    exit 0
fi

cwd=$(pwd)
echo "CWD: ${cwd}"
echo

echo "logfile: ${logfile}"
echo "finfile: ${finfile}"
echo

echo "Task CMD: ${task_cmd}"
echo
time eval "$task_cmd" 2>&1 | tee "$logfile"
task_return_code=$?

echo
echo "Task return code: ${task_return_code}"
if [ ! -f "$logfile" ]; then
    log_error=''
    echo "WARNING! Logfile does not exist: ${logfile}"
    echo "Cannot check logfile for potential errmsgs from task"
else
    log_error=$(grep -i -m1 'error' "$logfile")
    if [ -n "$log_error" ]; then
        echo "Found errmsg in logfile (${logfile}):"
        echo "$log_error"
    else
        echo "Found no errmsg in logfile (${logfile})"
    fi
fi
echo

# Create finfile if Matlab command exited without error
if (( task_return_code == 0 )) && [ -z "$log_error" ]; then
    echo "Considering run successful due to [task return code of zero] AND [no errmsgs found in logfile]"
    echo "Creating finfile: ${finfile}"
    touch "$finfile"
else
    echo "Considering run unsuccessful due to either [non-zero task return code] OR [errmsgs found in logfile]"
    echo "Will not create finfile"
fi

echo
echo "Done"
