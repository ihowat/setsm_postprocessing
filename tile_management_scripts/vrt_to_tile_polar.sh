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

#source ~/.bashrc; conda activate pgc

## Bash settings
set -uo pipefail


## Arguments

# Envvar arguments if run with job scheduler
set +u
src_vrt="$ARG_SRC_VRT"
set -u
if [ -n "$src_vrt" ]; then
    out_dir="$ARG_OUT_DIR"
    predictor="$ARG_PREDICTOR"
    round_float="$ARG_ROUND_FLOAT"
    resampling="$ARG_RESAMPLING"
    resolution="$ARG_RESOLUTION"
    buffer="$ARG_BUFFER"
    out_suffix="$ARG_OUT_SUFFIX"
    tilename="$ARG_TILENAME"
    domain="$ARG_DOMAIN"
fi

if [ -z "$src_vrt" ]; then
    if (( $# != 9 )); then
        echo "Script requires 9 arguments: SRC_VRT OUT_DIR PREDICTOR(3=float|2=int|1=off) ROUND_FLOAT(yes|no) RESAMPLING(near|bilinear|cubic) RESOLUTION(meters) BUFFER(meters) OUT_SUFFIX TILENAME(like '01_01')"
        exit 1
    fi
    # Command line arguments if run in bash
    src_vrt="$1"; shift
    out_dir="$1"; shift
    predictor="$1"; shift
    round_float="$1"; shift
    resampling="$1"; shift
    resolution="$1"; shift
    buffer="$1"; shift
    out_suffix="$1"; shift
    tilename="$1"; shift
    domain="$1"; shift
fi


## Validate arguments
if [ ! -f "$src_vrt" ]; then
    echo "SRC_VRT file does not exist: ${src_vrt}"
    exit 1
fi
if [ ! -d "$out_dir" ]; then
    echo "OUT_DIR directory does not exist: ${out_dir}"
    exit 1
fi
if ! { [ "$predictor" = 1 ] || [ "$predictor" = 2 ] || [ "$predictor" = 3 ]; }; then
    echo "PREDICTOR must be one of the following, but was '${predictor}': 1 2 3"
    exit 1
fi
if ! { [ "$round_float" = 'yes' ] || [ "$round_float" = 'no' ]; }; then
    echo "ROUND_FLOAT must be one of the following, but was '${round_float}': yes no"
    exit 1
fi
if ! { [ "$resampling" = 'near' ] || [ "$resampling" = 'bilinear' ] || [ "$resampling" = 'cubic' ]; }; then
    echo "RESAMPLING must be one of the following, but was '${resampling}': near bilinear cubic"
    exit 1
fi
if ! echo "$out_suffix" | grep -q -E '^_.*\.tif'; then
    echo "OUT_SUFFIX ('${out_suffix}') does not match expected pattern: '_*.tif'"
    exit 1
fi
if ! echo "$tilename" | grep -q -E '^(utm[0-9]{2}[ns]_)?[0-9]{2}_[0-9]{2}$'; then
    echo "TILENAME ('${tilename}') does not match expected pattern: '(utm[0-9]{2}[ns]_)?[0-9]{2}_[0-9]{2}'"
    exit 1
fi
if ! { [ "$domain" = 'arctic' ] || [ "$domain" = 'antarctic' ]; }; then
    echo "DOMAIN must be one of the following, but was '${domain}': arctic antarctic"
    exit 1
fi


# Mosaic tile domain settings
if [ "$domain" = 'arctic' ]; then
    mosaic_target_srs="EPSG:3413"  # gdalwarp argument '-t_srs' value
    lower_left_xmin=-4000000
    lower_left_ymin=-4000000
elif [ "$domain" = 'antarctic' ]; then
    mosaic_target_srs="EPSG:3031"  # gdalwarp argument '-t_srs' value
    lower_left_xmin=-3000000
    lower_left_ymin=-3000000
fi


# Parse tile row and col
IFS=_ read -r tile_row tile_col <<< "$tilename"

# First find proper tile bounds
tile_xmin="$(( lower_left_xmin + (10#$tile_col - 1) * 100000 ))"
tile_ymin="$(( lower_left_ymin + (10#$tile_row - 1) * 100000 ))"
tile_xmax="$(( tile_xmin + 100000 ))"
tile_ymax="$(( tile_ymin + 100000 ))"

# Then adjust tile bounds for buffer
tile_xmin="$(( tile_xmin - buffer ))"
tile_ymin="$(( tile_ymin - buffer ))"
tile_xmax="$(( tile_xmax + buffer ))"
tile_ymax="$(( tile_ymax + buffer ))"


out_tile="${out_dir}/${tilename}${out_suffix}"
if [ -f "$out_tile" ]; then
    echo "Removing existing output tile file: ${out_tile}"
    rm "$out_tile"
fi

if [ "$round_float" = 'no' ]; then
    set -x
    gdalwarp -overwrite -t_srs "$mosaic_target_srs" -tap -tr ${resolution} ${resolution} -te "$tile_xmin" "$tile_ymin" "$tile_xmax" "$tile_ymax" -of COG -co resampling=${resampling} -co overviews=ignore_existing -co compress=lzw -co predictor=${predictor} -co bigtiff=yes -co num_threads=all_cpus "$src_vrt" "$out_tile"
    cmd_status=$?
    set +x
    if (( cmd_status != 0 )); then
        echo "ERROR: gdalwarp returned non-zero exit status (${cmd_status})"
    fi

elif [ "$round_float" = 'yes' ]; then
    out_tile_tmp1="${out_dir}/${tilename}${out_suffix%.*}_tmp1.tif"
    out_tile_tmp2="${out_dir}/${tilename}${out_suffix%.*}_tmp2.tif"

    set -x
    gdalwarp -overwrite -t_srs "$mosaic_target_srs" -tap -tr ${resolution} ${resolution} -te "$tile_xmin" "$tile_ymin" "$tile_xmax" "$tile_ymax" -of GTiff -r ${resampling} -co tiled=yes -co compress=lzw -co bigtiff=yes -co num_threads=all_cpus "$src_vrt" "$out_tile_tmp1"
    cmd_status=$?
    set +x
    if (( cmd_status != 0 )); then
        echo "ERROR: gdalwarp returned non-zero exit status (${cmd_status})"
    fi

    nodata_val=$(gdalinfo "$out_tile_tmp1" | grep 'NoData Value' | cut -d= -f2)
    if [ -n "$nodata_val" ]; then
        nodata_arg="--NoDataValue=${nodata_val}"
    else
        nodata_arg=''
    fi

    set -x
    gdal_calc.py --overwrite -A "$out_tile_tmp1" --outfile="$out_tile_tmp2" --calc='round_(A*128.0)/128.0' ${nodata_arg} --co tiled=yes --co compress=lzw --co bigtiff=yes --co num_threads=all_cpus
    cmd_status=$?
    set +x
    if (( cmd_status != 0 )); then
        echo "ERROR: gdal_calc.py returned non-zero exit status (${cmd_status})"
    fi

    if [ -z "$nodata_arg" ]; then
        set -x
        gdal_edit.py -unsetnodata "$out_tile_tmp2"
        cmd_status=$?
        set +x
        if (( cmd_status != 0 )); then
            echo "ERROR: gdal_edit.py returned non-zero exit status (${cmd_status})"
        fi
    fi

    set -x
    gdal_translate -of COG -co resampling=${resampling} -co overviews=ignore_existing -co compress=lzw -co predictor=${predictor} -co bigtiff=yes -co num_threads=all_cpus "$out_tile_tmp2" "$out_tile"
    cmd_status=$?
    set +x
    if (( cmd_status != 0 )); then
        echo "ERROR: gdal_translate returned non-zero exit status (${cmd_status})"
    fi

    rm "$out_tile_tmp1" "$out_tile_tmp2"
fi

if [ -f "$out_tile" ]; then
    echo "Output tile file was created: ${out_tile}"
else
    echo "ERROR: output tile file was not created: ${out_tile}"
fi

echo "Done"
