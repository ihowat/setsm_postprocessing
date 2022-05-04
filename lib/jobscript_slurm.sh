#!/bin/bash

#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 2
#SBATCH --mem=8G
#SBATCH -o %x.o%j

echo ________________________________________
echo
echo SLURM Job Log
echo Start time: $(date)
echo
echo Job name: $SLURM_JOB_NAME
echo Job ID: $SLURM_JOBID
echo Submitted by user: $USER
echo User effective group ID: $(id -ng)
echo
echo SLURM account used: $SLURM_ACCOUNT
echo Hostname of submission: $SLURM_SUBMIT_HOST
echo Submitted to cluster: $SLURM_CLUSTER_NAME
echo Submitted to node: $SLURMD_NODENAME
echo Cores on node: $SLURM_CPUS_ON_NODE
echo Requested cores per task: $SLURM_CPUS_PER_TASK
echo Requested cores per job: $SLURM_NTASKS
echo Requested walltime: $SBATCH_TIMELIMIT
#echo Nodes assigned to job: $SLURM_JOB_NODELIST
#echo Running node index: $SLURM_NODEID
echo
echo Running on hostname: $HOSTNAME
echo Parent PID: $PPID
echo Process PID: $$
echo
working_dir="$SLURM_SUBMIT_DIR"

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


## Environment variables passed in through job submit command
init_env_script_path="$JOBSCRIPT_INIT_ENV_SCRIPT_PATH"
init_env_requests="$JOBSCRIPT_INIT_ENV_REQUESTS"
init_env_requests="${init_env_requests//\@COMMA\@/,}"; init_env_requests="${task_cmd//\@SPACE\@/ }";

set -uo pipefail

## The whole task command must be passed in as an environment variable named 'TASK_CMD'.
# With PBS, use the qsub '-v' option.
# With SLURM, use the sbatch '--export' option.
# In your Python script, replace commas with '@COMMA@' and spaces with '@SPACE@' (no quotes)
# in the string you send in with the TASK_CMD variable.
task_cmd="$JOBSCRIPT_TASK_CMD";
task_cmd="${task_cmd//\@COMMA\@/,}"; task_cmd="${task_cmd//\@SPACE\@/ }";


## Configure environment loading below:

#envload_cmd="module load matlab/2019a"
##envload_cmd=""
#if [ -n "$envload_cmd" ]; then
#    echo "Loading Matlab environment with the following command:"
#    echo -- "$envload_cmd"
#    eval "$envload_cmd"
#fi
#
##envload_cmd="source ~/.bashrc ; conda activate pgc"
#envload_cmd=""
#if [ -n "$envload_cmd" ]; then
#    echo "Loading Python/GDAL environment with the following command:"
#    echo -- "$envload_cmd"
#    eval "$envload_cmd"
#fi

## Or run requested environment(s) through init script
if [ -n "$init_env_script_path" ] && [ -n "$init_env_requests" ]; then
    bash "$init_env_script_path" ${init_env_requests}
    rc=$?
    if (( rc != 0 )); then
        echo
        echo "Non-zero return code (${rc}) from environment init script"
        echo "Exiting without executing task"
        exit 1
    fi
fi

echo
echo "Executing task command:"
echo -- "$task_cmd"
echo ________________________________________________________
echo

time eval "$task_cmd"
