#!/bin/bash

#PBS -l nodes=1:ppn=32:xe,gres=shifter
#PBS -l walltime=96:00:00
#PBS -v CRAY_ROOTFS=SHIFTER,UDI="ubuntu:xenial"
#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out
#PBS -m n
#PBS -q high


# Compute node parameters
JOB_ID="$PBS_JOBID"
CORES_PER_NODE="$PBS_NUM_PPN"

# Standard BW wrappers
APRUN_PREFIX="aprun -b -N 1 -d ${CORES_PER_NODE} -cc none --"
BWPY_PREFIX="bwpy-environ --"  # bwpy doesn't work within the Matlab Shifter image...
                               # It's possible to use bwpy outside of the Matlab process aprun,
                               # but a certain procedure must be followed to load/unload
                               # the Shifter image in the jobscript for running Matlab.

# Matlab configuration for BW
MATLAB_PROGRAM="/projects/sciteam/bazu/matlab/R2020a/bin/matlab"
export MATLABHOST=$(printf 'nid%05d' "$(head -n1 "$PBS_NODEFILE")")
export LM_LICENSE_FILE="27000@matlab-pgc.cse.umn.edu"

# Core Matlab run settings
MATLAB_SETTINGS="-nodisplay -nodesktop -nosplash"
MATLAB_USE_PARPOOL=true

# Load Python/GDAL env if needed
set +u
source /projects/sciteam/bazu/tools/miniconda3/bin/activate /projects/sciteam/bazu/tools/miniconda3/envs/gdal2
export PATH="${PATH}:/projects/sciteam/bazu/tools/miniconda3/envs/gdal2/bin/"
set -u

### EDIT THIS BLOCK ###
# Matlab settings specific to the script you want to run
MATLAB_SCRIPTS_DIR="/projects/sciteam/bazu/tools/setsm_postprocessing_pgc"
MATLAB_SCRIPTS_DIR2="/projects/sciteam/bazu/tools/setsm_postprocessing4"
MATLAB_WORKING_DIR="/scratch/sciteam/GS_bazu/mosaic_data/matlab_working_dir"
MATLAB_TEMP_DIR="/scratch/sciteam/GS_bazu/mosaic_data/matlab_temp_dir"
#MATLAB_LOGFILE="/scratch/sciteam/GS_bazu/mosaic_data/logs/example.log"
MATLAB_LOGFILE=''


if [ "$MATLAB_USE_PARPOOL" = true ]; then
    job_working_dir="${MATLAB_WORKING_DIR}"
    job_temp_dir="${MATLAB_TEMP_DIR}/${JOB_ID}"
    matlab_parpool_init="\
pc = parcluster('local'); \
pc.JobStorageLocation = '${job_temp_dir}'; \
pc.NumWorkers = ${CORES_PER_NODE}; \
parpool(pc, ${CORES_PER_NODE});"
else
    job_working_dir="$MATLAB_TEMP_DIR"
    job_temp_dir="$MATLAB_TEMP_DIR"
    matlab_parpool_init="\
ps = parallel.Settings; \
ps.Pool.AutoCreate = false;"
fi


### ENABLE AND EDIT ONE OF THE FOLLOWING TWO BLOCKS ###

# For interactive MATLAB job
matlab_cmd=''

## For running MATLAB script
#matlab_cmd="\
#addpath('${MATLAB_SCRIPTS_DIR}'); \
#addpath('${MATLAB_SCRIPTS_DIR2}'); \
#${matlab_parpool_init} \
#try; matlab_example; catch e; disp(getReport(e)); exit(1); end; exit(0);"\


task_cmd="${MATLAB_PROGRAM} ${MATLAB_SETTINGS}"
if [ -n "$matlab_cmd" ]; then
    task_cmd="${task_cmd} -r \"${matlab_cmd}\""
fi
task_cmd="${APRUN_PREFIX} bash -c '\
export LD_LIBRARY_PATH=\"/projects/sciteam/bazu/matlab/lib-GLIBC2.12:\${LD_LIBRARY_PATH}\"; \
$(echo "$task_cmd" | sed "s|'|'\"'\"'|g");'"


# Move to working folder.
# Working folder should be empty to reduce Matlab startup time.
if [ ! -d "$job_working_dir" ]; then
    mkdir -p "$job_working_dir"
fi
if [ ! -d "$job_temp_dir" ]; then
    mkdir -p "$job_temp_dir"
fi
cd "$job_working_dir"
cd_status=$?
if (( cd_status != 0 )); then
    echo "Failed to change to working dir: ${job_working_dir}"
    exit 0
fi


echo "Running command: ${task_cmd}"
if [ -n "$MATLAB_LOGFILE" ]; then
    time eval "$task_cmd" 2>&1 | tee "$MATLAB_LOGFILE"
else
    time eval "$task_cmd"
fi


if [ "$MATLAB_USE_PARPOOL" = true ] && [ -d "$job_temp_dir" ]; then
    echo "Removing Matlab temp dir for parallel tasks: ${job_temp_dir}"
    sleep 10s
    rm -rf "$job_temp_dir"
fi

echo
echo "Done"
