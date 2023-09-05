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

#source /home/husby036/installed/build/miniconda3/bin/activate /home/husby036/installed/build/miniconda3/envs/pgc

## Bash settings
set -uo pipefail


## Arguments

# Envvar arguments if run with job scheduler
set +u
out_file="$ARG_OUT_FILE"
set -u
if [ -n "$out_file" ]; then
    raster_a="$ARG_RASTER_A"
    raster_b="$ARG_RASTER_B"
    add_or_subtract="$ARG_ADD_OR_SUBTRACT"
    predictor="$ARG_PREDICTOR"
    round_float="$ARG_ROUND_FLOAT"
    resampling="$ARG_RESAMPLING"
    resolution="$ARG_RESOLUTION"
    adjust_res="$ARG_ADJUST_RES"
    adjust_res_degrees="$ARG_ADJUST_RES_DEGREES"
    overwrite="$ARG_OVERWRITE"

fi

if [ -z "$out_file" ]; then
    if (( $# != 8 )); then
        echo "Script requires 8 arguments: OUT_FILE RASTER_A RASTER_B ADD_OR_SUBTRACT(add|subtract) PREDICTOR(3=float|2=int|1=off) ROUND_FLOAT(yes|no) RESAMPLING(near|bilinear|cubic) RESOLUTION(meters)"
        exit 1
    fi
    # Command line arguments if run in bash
    out_file="$1"; shift
    raster_a="$1"; shift
    raster_b="$1"; shift
    add_or_subtract="$1"; shift
    predictor="$1"; shift
    round_float="$1"; shift
    resampling="$1"; shift
    resolution="$1"; shift
    adjust_res="$1"; shift
    adjust_res_degrees="$1"; shift
    overwrite="$1"; shift
fi


## Validate arguments
if [ ! -f "$raster_a" ]; then
    echo "RASTER_A file does not exist: ${raster_a}"
    exit 1
fi
if [ ! -f "$raster_b" ]; then
    echo "RASTER_B file does not exist: ${raster_b}"
    exit 1
fi
if ! { [ "$add_or_subtract" = 'add' ] || [ "$add_or_subtract" = 'subtract' ]; }; then
    echo "ADD_OR_SUBTRACT must be one of the following, but was '${add_or_subtract}': add subtract"
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
if ! { [ "$adjust_res" = 'true' ] || [ "$adjust_res" = 'false' ]; }; then
    echo "ADJUST_RES must be one of the following, but was '${adjust_res}': true false"
    exit 1
fi
if ! { [ "$adjust_res_degrees" = 'true' ] || [ "$adjust_res_degrees" = 'false' ]; }; then
    echo "ADJUST_RES_DEGREES must be one of the following, but was '${adjust_res_degrees}': true false"
    exit 1
fi
if ! { [ "$overwrite" = 'true' ] || [ "$overwrite" = 'false' ]; }; then
    echo "OVERWRITE must be one of the following, but was '${overwrite}': true false"
    exit 1
fi


if [ -f "$out_file" ]; then
    if [ "$overwrite" = false ]; then
        echo "Output file already exists: ${out_file}"
        exit 0
    fi
    echo "Removing existing output file: ${out_file}"
    rm "$out_file"
fi

if [ "$add_or_subtract" = 'add' ]; then
    calc_expr="A+B"
elif [ "$add_or_subtract" = 'subtract' ]; then
    calc_expr="A-B"
fi
if [ "$round_float" = 'yes' ]; then
    calc_expr="round_((${calc_expr})*128.0)/128.0"
fi

out_file_tmp="${out_file%.*}_tmp.tif"
raster_b_tmp_vrt="${out_file%.*}_tmp_rasterB.vrt"


if [ "$adjust_res" = true ]; then
    pixel_size=$(gdalinfo "$raster_a" | grep '^Pixel Size =')
    pixel_size_x=$(echo "$pixel_size" | grep -Eo '\(-?[0-9]*\.?[0-9]*,' | sed -r 's|[(),-]||g')
    pixel_size_y=$(echo "$pixel_size" | grep -Eo ',-?[0-9]*\.?[0-9]*\)' | sed -r 's|[(),-]||g')

    if [ "$adjust_res_degrees" = true ]; then
        module load glibc
        pixels_per_degree_lat=$(printf %.0f "$(echo "1 / ${pixel_size_y}" | bc -l)")
        pixels_per_degree_lon=$(printf %.0f "$(echo "1 / ${pixel_size_x}" | bc -l)")
        pixels_global_lat=$(echo "scale=0 ; 180 * ${pixels_per_degree_lat}" | bc)
        pixels_global_lon=$(echo "scale=0 ; 360 * ${pixels_per_degree_lon}" | bc)
        module unload glibc

        set -x
        gdalwarp -overwrite -of VRT -r bilinear -te "-180" "-90" 180 90 -ts "$pixels_global_lon" "$pixels_global_lat" "$raster_b" "$raster_b_tmp_vrt"
        cmd_status=$?
        set +x
        if (( cmd_status != 0 )); then
            echo "ERROR: gdalwarp returned non-zero exit status (${cmd_status})"
        fi
    else
        set -x
        gdalwarp -overwrite -of VRT -r bilinear -tap -tr "$pixel_size_x" "$pixel_size_y" "$raster_b" "$raster_b_tmp_vrt"
        cmd_status=$?
        set +x
        if (( cmd_status != 0 )); then
            echo "ERROR: gdalwarp returned non-zero exit status (${cmd_status})"
        fi
    fi

    raster_b="$raster_b_tmp_vrt"
fi

nodata_val=$(gdalinfo "$raster_a" | grep 'NoData Value' | cut -d= -f2)
if [ -n "$nodata_val" ]; then
    nodata_arg="--NoDataValue=${nodata_val}"
else
    nodata_arg=''
fi

set -x
gdal_calc.py --overwrite -A "$raster_a" -B "$raster_b" --extent=intersect --outfile="$out_file_tmp" --calc="$calc_expr" ${nodata_arg} --co tiled=yes --co compress=lzw --co bigtiff=yes --co num_threads=all_cpus
cmd_status=$?
set +x
if (( cmd_status != 0 )); then
    echo "ERROR: gdal_calc.py returned non-zero exit status (${cmd_status})"
fi

if [ -z "$nodata_arg" ]; then
    set -x
    gdal_edit.py -unsetnodata "$out_file_tmp"
    cmd_status=$?
    set +x
    if (( cmd_status != 0 )); then
        echo "ERROR: gdal_edit.py returned non-zero exit status (${cmd_status})"
    fi
fi

if [ -n "$resolution" ]; then
    resolution_args="-tr ${resolution} ${resolution}"
else
    resolution_args=''
fi
set -x
gdal_translate -of COG ${resolution_args} -co resampling=${resampling} -co overviews=ignore_existing -co compress=lzw -co predictor=${predictor} -co bigtiff=yes -co num_threads=all_cpus "$out_file_tmp" "$out_file"
cmd_status=$?
set +x
if (( cmd_status != 0 )); then
    echo "ERROR: gdal_translate returned non-zero exit status (${cmd_status})"
fi

rm "$out_file_tmp"
if [ -f "$raster_b_tmp_vrt" ]; then
    rm "$raster_b_tmp_vrt"
fi

if [ -f "$out_file" ]; then
    echo "Output file was created: ${out_file}"
else
    echo "ERROR: output file was not created: ${out_file}"
fi

echo "Done"
