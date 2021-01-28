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

task_cmd="${task_cmd//|COMMA|/,}"
echo "Task CMD: ${task_cmd}"
echo

job_logfile="/home/${USER}/${PBS_JOBNAME}.o$(echo "$PBS_JOBID" | cut -d'.' -f1)"
echo "job logfile: ${job_logfile}"
echo "finfile: ${finfile}"
echo

echo "$task_cmd"
time eval "$task_cmd"
task_return_code=$?

echo

echo "Task return code: ${task_return_code}"
if [ ! -f "$job_logfile" ]; then
    log_error=''
    echo "WARNING! Job logfile does not exist: ${job_logfile}"
    echo "Cannot check job logfile for potential errmsgs from task"
else
    log_error=$(grep -i -m1 'error' "$job_logfile")
    if [ -n "$log_error" ]; then
        echo "Found errmsg in job logfile (${job_logfile}):"
        echo "$log_error"
    else
        echo "Found no errmsg in job logfile (${job_logfile})"
    fi
fi
echo

# create finfile if matlab command exited without error
if (( task_return_code == 0 )) && [ -z "$log_error" ]; then
    echo "Considering run successful due to [task return code of zero] AND [no errmsgs found in job logfile]"
    echo "Creating finfile: ${finfile}"
    touch "$finfile"
else
    echo "Considering run unsuccessful due to either [non-zero task return code] OR [errmsgs found in job logfile]"
    echo "Will not create finfile"
fi

echo
echo "Done"
