#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=2,mem=8gb
#PBS -m n
#PBS -k oe
#PBS -j oe

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

# Environment variables passed in through job submit command
init_env_script_path="$JOBSCRIPT_INIT_ENV_SCRIPT_PATH"
init_env_requests="$JOBSCRIPT_INIT_ENV_REQUESTS"
init_env_requests="${init_env_requests//\@COMMA\@/,}"; init_env_requests="${init_env_requests//\@SPACE\@/ }";

if [ -n "$init_env_script_path" ] && [ -n "$init_env_requests" ]; then
    source "$init_env_script_path" ${init_env_requests}
    rc=$?
    if (( rc != 0 )); then
        echo "Non-zero return code (${rc}) from environment init script"
        echo "Exiting without executing task"
        exit 1
    fi
fi


## Run task command(s)

if [ -n "$JOBSCRIPT_TASK_BUNDLE_FILE" ]; then
    task_run_type='batch'
else
    task_run_type='single'
fi

if [ "$task_run_type" = 'single' ]; then
    ## The whole task command must be passed in as an environment variable named 'TASK_CMD'.
    # With PBS, use the qsub '-v' option.
    # With SLURM, use the sbatch '--export' option.
    # In your Python script, replace commas with '@COMMA@' and spaces with '@SPACE@' (no quotes)
    # in the string you send in with the TASK_CMD variable.
    task_cmd="$JOBSCRIPT_TASK_CMD";
    task_cmd="${task_cmd//\@COMMA\@/,}"; task_cmd="${task_cmd//\@SPACE\@/ }";

    echo ________________________________________________________
    echo
    echo $(date) -- "Executing task command:"
    printf '%s\n' "$task_cmd"
    echo ________________________________________________________
    echo

    time eval "$task_cmd"

elif [ "$task_run_type" = 'batch' ]; then
    task_bundle_file="$JOBSCRIPT_TASK_BUNDLE_FILE"
    if [ -n "$JOBSCRIPT_TASK_BUNDLE_MODE" ]; then
        batch_run_mode="$JOBSCRIPT_TASK_BUNDLE_MODE"
    else
        batch_run_mode='serial'
    fi

    num_tasks=$(wc -l "$task_bundle_file" | cut -d' ' -f1)

    echo "Will run batch of ${num_tasks} tasks in ${batch_run_mode} mode from bundle file: ${task_bundle_file}"

    task_num=0
    while IFS= read -r task_cmd; do
        ((task_num++))

        echo ________________________________________________________
        echo
        echo $(date) -- "Executing task command (${task_num}/${num_tasks}):"
        printf '%s\n' "$task_cmd"
        echo ________________________________________________________
        echo

        if [ "$batch_run_mode" = 'serial' ]; then
            time eval "$task_cmd"
        elif [ "$batch_run_mode" = 'parallel' ]; then
            task_cmd="${task_cmd}; echo; echo; Finished task (${task_num}/${num_tasks})"
            time eval "$task_cmd" &
        fi
    done < "$task_bundle_file"
fi
