#!/bin/bash

##PBS -l nodes=1:ppn=2,walltime=200:00:00
##PBS -m n
##PBS -k oe
##PBS -j oe
##PBS -q batch

#SBATCH --time 200:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 2
#SBATCH --mem=10G
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

source /home/husby036/installed/build/miniconda3/bin/activate /home/husby036/installed/build/miniconda3/envs/pgc

## Bash settings
set -uo pipefail


tile_math_script="$ARG_JOBSCRIPT"
tile_row="$ARG_ROWNUM"

#for tile_col in {1..80}; do  # Arctic has 80 cols, Antarctic has 60 cols
#for tile_col in {1..60}; do  # Arctic has 80 cols, Antarctic has 60 cols
for dem_file in /mnt/pgc/data/elev/dem/copernicus-dem-30m/COP30_hh/Copernicus_DSM_COG_10_${tile_row}_*_DEM.tif; do

#    tilename=$(printf '%02d_%02d' "$tile_row" "$tile_col")
#    echo "Running tile math for tile: ${tilename}"

    echo "Running tile math for tile: ${dem_file}"

#    export ARG_OUT_FILE="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles/${tilename}_10m_cop30_wgs84.tif"
#    export ARG_RASTER_A="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles/${tilename}_10m_cop30.tif"
##    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/egm2008/mosaic/arctic_tiles/${tilename}_10m_egm08.tif"
#    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/egm2008/us_nga_egm2008_1_3413_arctic-10m.vrt"
#    export ARG_ADD_OR_SUBTRACT=add
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_ADJUST_RES=false
#    export ARG_ADJUST_RES_DEGREES=false
#    export ARRG_OVERWRITE=false

#    export ARG_OUT_FILE="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/results/output_tiles/${tilename}/${tilename}_10m_dem_diff-cop30.tif"
#    export ARG_RASTER_A="/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/results/output_tiles/${tilename}/${tilename}_10m_dem.tif"
#    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles/${tilename}_10m_cop30_wgs84.tif"
#    export ARG_ADD_OR_SUBTRACT=subtract
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_ADJUST_RES=false
#    export ARG_ADJUST_RES_DEGREES=false
#    export ARRG_OVERWRITE=false

#    export ARG_OUT_FILE="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles/${tilename}_10m_cop30_wgs84.tif"
#    export ARG_RASTER_A="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles/${tilename}_10m_cop30.tif"
##    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/egm2008/mosaic/rema_tiles/${tilename}_10m_egm08.tif"
#    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/egm2008/us_nga_egm2008_1_3031_rema-10m.vrt"
#    export ARG_ADD_OR_SUBTRACT=add
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_ADJUST_RES=false
#    export ARG_ADJUST_RES_DEGREES=false
#    export ARRG_OVERWRITE=false

#    export ARG_OUT_FILE="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles/${tilename}/${tilename}_10m_dem_diff-cop30.tif"
#    export ARG_RASTER_A="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles/${tilename}/${tilename}_10m_dem.tif"
#    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles/${tilename}_10m_cop30_wgs84.tif"
#    export ARG_ADD_OR_SUBTRACT=subtract
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_ADJUST_RES=false
#    export ARG_ADJUST_RES_DEGREES=false
#    export ARRG_OVERWRITE=false

    out_fname=$(basename "${dem_file%.*}_wgs84.tif")
    export ARG_OUT_FILE="/mnt/pgc/data/elev/dem/copernicus-dem-30m/COP30_hh_wgs84-height/${out_fname}"
    export ARG_RASTER_A="$dem_file"
    export ARG_RASTER_B="/mnt/pgc/data/elev/dem/egm2008/us_nga_egm2008_1_4326_1-arcsec.vrt"
    export ARG_ADD_OR_SUBTRACT=add
    export ARG_PREDICTOR=3
    export ARG_ROUND_FLOAT=yes
    export ARG_RESAMPLING=bilinear
    export ARG_RESOLUTION=
    export ARG_ADJUST_RES=true
    export ARG_ADJUST_RES_DEGREES=true
    export ARG_OVERWRITE=false

    bash "$tile_math_script"
    echo

done
