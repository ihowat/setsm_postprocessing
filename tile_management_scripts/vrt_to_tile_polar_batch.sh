#!/bin/bash

#PBS -l nodes=1:ppn=2,walltime=200:00:00
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

source /home/husby036/installed/build/miniconda3/bin/activate /home/husby036/installed/build/miniconda3/envs/pgc

## Bash settings
set -uo pipefail


vrt_to_tile_script="$ARG_JOBSCRIPT"
tile_row="$ARG_ROWNUM"  # Arctic has X rows, Antarctic has 65 rows

#for tile_col in {1..80}; do  # Arctic has 80 cols, Antarctic has 60 cols
for tile_col in {1..60}; do  # Arctic has 80 cols, Antarctic has 60 cols
    tilename=$(printf '%02d_%02d' "$tile_row" "$tile_col")

    echo "Running VRT to tile export for tile: ${tilename}"

#    export ARG_SRC_VRT="/mnt/pgc/data/thematic/landcover/esa_worldcover_2020/mosaics/esa_worldcover_2020_10m_adem-domain.vrt"
#    export ARG_OUT_DIR="/mnt/pgc/data/thematic/landcover/esa_worldcover_2020/mosaics/arctic_tiles"
#    export ARG_PREDICTOR=2
#    export ARG_ROUND_FLOAT=no
#    export ARG_RESAMPLING=near
#    export ARG_RESOLUTION=10
#    export ARG_BUFFER=100
#    export ARG_OUT_SUFFIX=_10m_cover.tif
#    export ARG_TILENAME="$tilename"
#    export ARG_DOMAIN=arctic

#    export ARG_SRC_VRT="/mnt/pgc/data/elev/dem/egm2008/mosaic/us_nga_egm2008_1_arctic-100m.vrt"
#    export ARG_OUT_DIR="/mnt/pgc/data/elev/dem/egm2008/mosaic/arctic_tiles"
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_BUFFER=100
#    export ARG_OUT_SUFFIX=_10m_egm08.tif
#    export ARG_TILENAME="$tilename"
#    export ARG_DOMAIN=arctic

#    export ARG_SRC_VRT="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/copernicus-dem-30m_adem-domain.vrt"
#    export ARG_OUT_DIR="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles"
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_BUFFER=100
#    export ARG_OUT_SUFFIX=_10m_cop30.tif
#    export ARG_TILENAME="$tilename"
#    export ARG_DOMAIN=arctic

    export ARG_SRC_VRT="/mnt/pgc/data/elev/dem/egm2008/mosaic/us_nga_egm2008_1_rema-100m.vrt"
    export ARG_OUT_DIR="/mnt/pgc/data/elev/dem/egm2008/mosaic/rema_tiles"
    export ARG_PREDICTOR=3
    export ARG_ROUND_FLOAT=yes
    export ARG_RESAMPLING=bilinear
    export ARG_RESOLUTION=10
    export ARG_BUFFER=100
    export ARG_OUT_SUFFIX=_10m_egm08.tif
    export ARG_TILENAME="$tilename"
    export ARG_DOMAIN=antarctic

#    export ARG_SRC_VRT="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/copernicus-dem-30m_rema-domain.vrt"
#    export ARG_OUT_DIR="/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles"
#    export ARG_PREDICTOR=3
#    export ARG_ROUND_FLOAT=yes
#    export ARG_RESAMPLING=bilinear
#    export ARG_RESOLUTION=10
#    export ARG_BUFFER=100
#    export ARG_OUT_SUFFIX=_10m_cop30.tif
#    export ARG_TILENAME="$tilename"
#    export ARG_DOMAIN=antarctic

    bash "$vrt_to_tile_script"
    echo

done
