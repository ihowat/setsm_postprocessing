#!/bin/bash

#PBS -l walltime=100:00:00,nodes=1:ppn=8,mem=48gb
#PBS -m n
#PBS -k oe
#PBS -j oe
#PBS -q old

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
CORES_PER_NODE=$PBS_NUM_PPN

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
echo "var12: ${p12}"

dstfile="$p6"
finfile="$p10"
echo "dstfile: var6: ${dstfile}"
echo "finfile: var10: ${finfile}"
job_logfile="/home/${USER}/${PBS_JOBNAME}.o$(echo "$PBS_JOBID" | cut -d'.' -f1)"
echo "job logfile: ${job_logfile}"

echo

export BWPY_PREFIX=""

# Check if quad arg is present
if [ "${p9}" == 'null' ]; then
    echo "Quad arg not present. Running full tile"

    echo matlab -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); warning('off','all'); parpool($CORES_PER_NODE); [x0,x1,y0,y1]=getTileExtents('${p7}','${p8}'); disp(x0); ${p3}('${p4}',${p5},'${p6}','projection','${p11}','version','${p12}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);"
    time matlab -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); warning('off','all'); parpool($CORES_PER_NODE); [x0,x1,y0,y1]=getTileExtents('${p7}','${p8}'); disp(x0); ${p3}('${p4}',${p5},'${p6}','projection','${p11}','version','${p12}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);"
    task_return_code=$?

else
    echo "Running quadrant ${p9}"

    echo matlab -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); warning('off','all'); parpool($CORES_PER_NODE); [x0,x1,y0,y1]=getTileExtents('${p7}','${p8}','quadrant','${p9}'); ${p3}('${p4}',${p5},'${p6}','projection','${p11}','version','${p12}','quadrant','${p9}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);"
    time matlab -nodisplay -nosplash -r "try; addpath('${p1}'); addpath('${p2}'); warning('off','all'); parpool($CORES_PER_NODE); [x0,x1,y0,y1]=getTileExtents('${p7}','${p8}','quadrant','${p9}'); ${p3}('${p4}',${p5},'${p6}','projection','${p11}','version','${p12}','quadrant','${p9}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);"
    task_return_code=$?

fi

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
#    if [ -f "$dstfile" ]; then
#        echo "Removing dstfile: ${dstfile}"
#        rm "$dstfile"
#    fi
fi

echo
echo "Done"
