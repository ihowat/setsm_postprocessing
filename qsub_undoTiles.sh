#!/bin/bash

##PBS -l walltime=100:00:00,nodes=1:ppn=4,mem=30gb
##PBS -m n
##PBS -k oe
##PBS -j oe
##PBS -q batch

#SBATCH --time 10:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 6
#SBATCH --mem=30G
#SBATCH --partition=batch
#SBATCH -o %x.o%j

#echo ________________________________________________________
#echo
#echo PBS Job Log
#echo Start time: $(date)
#echo
#echo Job name: $PBS_JOBNAME
#echo Job ID: $PBS_JOBID
#echo Submitted by user: $USER
#echo User effective group ID: $(id -ng)
#echo
#echo Hostname of submission: $PBS_O_HOST
#echo Submitted to cluster: $PBS_SERVER
#echo Submitted to queue: $PBS_QUEUE
#echo Requested nodes per job: $PBS_NUM_NODES
#echo Requested cores per node: $PBS_NUM_PPN
#echo Requested cores per job: $PBS_NP
##echo Node list file: $PBS_NODEFILE
##echo Nodes assigned to job: $(cat $PBS_NODEFILE)
##echo Running node index: $PBS_O_NODENUM
#echo
#echo Running on hostname: $HOSTNAME
#echo Parent PID: $PPID
#echo Process PID: $$
#echo
#working_dir="$PBS_O_WORKDIR"
#
#echo Number of physical cores: $(grep "^core id" /proc/cpuinfo | sort -u | wc -l)
#echo Number of virtual cores: $(nproc --all)
#mem_report_cmd="free -g"
#echo "Memory report [GB] ('${mem_report_cmd}'):"
#echo -------------------------------------------------------------------------
#${mem_report_cmd}
#echo -------------------------------------------------------------------------
#echo

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

set -uo pipefail

#module load gdal/2.1.3
module load matlab/2019a

## Arguments/Options
tiledir="$ARG_TILEDIR"

## Validate arguments
if [ ! -d "$tiledir" ]; then
    echo "Tiledir does not exist: ${tiledir}"
    exit 1
fi

matlab_cmd="try; addpath('/mnt/pgc/data/common/repos/setsm_postprocessing4'); matfile_list=dir(['${tiledir}','/*.mat']); matfile_paths=fullfile({matfile_list.folder},{matfile_list.name}); undoTile(matfile_paths); catch e; disp(getReport(e)); exit(1); end; exit(0)"
#matlab_cmd="try; addpath('/mnt/pgc/data/scratch/erik/repos/setsm_postprocessing4'); matfile_list=dir(['${tiledir}','/*.mat']); matfile_paths=fullfile({matfile_list.folder},{matfile_list.name}); undoTile(matfile_paths); catch e; disp(getReport(e)); exit(1); end; exit(0)"

echo "Argument tile directory: ${tiledir}"
echo "Matlab command: \"${matlab_cmd}\""

time matlab -nojvm -nodisplay -nosplash -r "$matlab_cmd"
