#!/bin/bash

## PGC settings
##PBS -l walltime=200:00:00,nodes=1:ppn=16,mem=64gb
##PBS -m n
##PBS -k oe
##PBS -j oe
##PBS -q old

## BW settings
##PBS -l nodes=1:ppn=16:xe,gres=shifter
##PBS -l walltime=96:00:00
##PBS -v CRAY_ROOTFS=SHIFTER,UDI="ubuntu:xenial"
##PBS -e $PBS_JOBID.err
##PBS -o $PBS_JOBID.out
##PBS -m n
##PBS -q high


set -uo pipefail

set +u; tileName="$ARG_TILENAME" ; set -u
system="$ARG_SYSTEM"
bst_jobscript="$ARG_BST_JOBSCRIPT"
mst_jobscript="$ARG_MST_JOBSCRIPT"
mst_pyscript="$ARG_MST_PYSCRIPT"
output_tiles_dir="$ARG_OUTPUT_TILES_DIR"
project="$ARG_PROJECT"
waterTileDir="$ARG_WATERTILEDIR"
make_10m_only="$ARG_MAKE_10M_ONLY"

if [ -z "$tileName" ]; then
    tileName="$1"
fi
if [ -z "$tileName" ]; then
    echo "argument 'tileName' not supplied, exiting"
    exit 0
fi

tile_subtiles_dir="${output_tiles_dir}/${tileName}/subtiles/"
bst_finfile_10m="${output_tiles_dir}/${tileName}/subtiles_10m"
bst_finfile_2m="${output_tiles_dir}/${tileName}/subtiles_2m"


#if [ "$system" = 'pgc' ]; then
#    module load gdal/2.1.3
#    python_exec='python'
#elif [ "$system" = 'bw' ]; then
#    module load bwpy/2.0.2
#    python_exec='bwpy-environ -- python'
#fi

bst_cmd="bash \"${bst_jobscript}\""
#mst_cmd_base="${python_exec} \"${mst_pyscript}\" \"${output_tiles_dir}\" ${tileName} --project ${project} --chain-mst-jobscript \"${mst_jobscript}\" --submit"
#mst_cmd_10m="${mst_cmd_base} 10"
#mst_cmd_2m="${mst_cmd_base} 2 --quads"
mst_cmd_base="bash \"${mst_jobscript}\""
mst_cmd_10m="${mst_cmd_base} 10"
mst_cmd_2m="${mst_cmd_base} 2"
quads_arr=( '1_1' '1_2' '2_1' '2_2' )


run_bst=true
if [ "$make_10m_only" = true ]; then
    if [ -f "$bst_finfile_10m" ] || [ -f "$bst_finfile_2m" ]; then
        run_bst=false
    fi
else
    if [ -f "$bst_finfile_2m" ]; then
        run_bst=false
    fi
fi
if [ "$run_bst" = true ]; then
    bst_txt="BST then "
else
    bst_txt=''
fi
if [ "$make_10m_only" = true ]; then
    echo "Will run ${bst_txt}MST 10m for tile: ${tileName}"
else
    echo "Will run ${bst_txt}MST 10m and 2m for tile: ${tileName}"
fi


if [ "$run_bst" = true ]; then
    cmd="$bst_cmd"
    echo "Running BST process with command: ${cmd}"
    export ARG_TILENAME="$tileName"
    export ARG_WATERTILEDIR="$waterTileDir"
    eval "$cmd"
    cmd_status=$?
    if (( cmd_status != 2 )); then
        echo "BST jobscript exited with non-success(2) exit status (${cmd_status})"
        exit
    fi
fi

cmd="${mst_cmd_10m}"
echo "Running 10m MST process with command: ${cmd}"
export ARG_TILENAME="$tileName"
eval "$cmd"
cmd_status=$?
if (( cmd_status != 2 )); then
    echo "MST jobscript exited with non-success(2) exit status (${cmd_status})"
    exit
fi

if [ "$make_10m_only" = false ]; then
    cmd="${mst_cmd_2m}"
    for quad in "${quads_arr[@]}"; do
        quadTileName="${tileName}_${quad}"
        echo "Running 2m MST process with command: ${cmd}"
        export ARG_TILENAME="$quadTileName"
        eval "$cmd"
        cmd_status=$?
        if (( cmd_status != 2 )); then
            echo "MST jobscript exited with non-success(2) exit status (${cmd_status})"
            exit
        fi
    done
fi


if [ ! -d "$tile_subtiles_dir" ]; then
    echo "Cannot remove tile subtiles directory, it does not exist: ${tile_subtiles_dir}"
else
    echo "Removing tile subtiles directory: ${tile_subtiles_dir}"
    rm -rf "$tile_subtiles_dir"
fi
