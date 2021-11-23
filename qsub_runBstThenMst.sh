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
keep_subtiles="$ARG_KEEP_SUBTILES"

if [ -z "$tileName" ]; then
    tileName="$1"
fi
if [ -z "$tileName" ]; then
    echo "argument 'tileName' not supplied, exiting"
    exit 0
fi

tile_subtiles_dir="${output_tiles_dir}/${tileName}/subtiles/"
bst_finfile_10m="${output_tiles_dir}/${tileName}/subtiles_10m.fin"
bst_finfile_2m="${output_tiles_dir}/${tileName}/subtiles_2m.fin"
mst_finfile_10m="${output_tiles_dir}/${tileName}/${tileName}_10m.fin"
mst_finfile_2m_template="${output_tiles_dir}/${tileName}/<quadTileName>_2m.fin"


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


echo
echo "Checking progress of MST step..."
run_mst=false
if [ -f "$mst_finfile_10m" ]; then
    echo "MST 10m finfile exists: ${mst_finfile_10m}"
else
    echo "MST 10m finfile does not exist: ${mst_finfile_10m}"
    run_mst=true
fi
if [ "$make_10m_only" = false ]; then
    all_2m_mst_finfiles_exist=true
    for quad in "${quads_arr[@]}"; do
        quadTileName="${tileName}_${quad}"
        mst_finfile_2m="${mst_finfile_2m_template//<quadTileName>/${quadTileName}}"
        if [ -f "$mst_finfile_2m" ]; then
            echo "MST 2m quad finfile exists: ${mst_finfile_2m}"
        else
            echo "MST 2m quad finfile does not exist: ${mst_finfile_2m}"
            all_2m_mst_finfiles_exist=false
            run_mst=true
        fi
    done
    if [ "$all_2m_mst_finfiles_exist" = true ]; then
        echo "All MST 2m quad finfiles exist, so --make-10m-only BST option will be applied"
        bst_cmd="${bst_cmd} --make-10m-only"
        make_10m_only=true
    fi
fi
if [ "$run_mst" = true ]; then
    echo "MST step is incomplete as expected"
else
    echo "All MST finfiles already exist, so BST and MST steps will not be run!"
    echo "Exiting"
    exit
fi


echo
echo "Checking progress of BST step..."
run_bst=false
if [ ! -d "$tile_subtiles_dir" ]; then
    echo "BST subtiles folder does not exist: ${tile_subtiles_dir}"
    run_bst=true
    if [ -f "$bst_finfile_2m" ]; then
        echo "Removing existing BST 2m finfile: ${bst_finfile_2m}"
        rm -vf "$bst_finfile_2m"
    fi
    if [ -f "$bst_finfile_10m" ]; then
        echo "Removing existing BST 10m finfile: ${bst_finfile_10m}"
        rm -vf "$bst_finfile_10m"
    fi
else
    echo "BST subtiles folder exists: ${tile_subtiles_dir}"
    if [ -f "$bst_finfile_2m" ]; then
        echo "BST 2m finfile exists, so BST step will be skipped: ${bst_finfile_2m}"
        run_bst=false
    elif [ "$make_10m_only" = true ]; then
        if [ -f "$bst_finfile_10m" ]; then
            echo "BST 10m finfile exists, so BST step will be skipped: ${bst_finfile_10m}"
            run_bst=false
        else
            echo "BST 10m finfile does not exist, so BST step will be run: ${bst_finfile_10m}"
            run_bst=true
        fi
    else
        echo "BST 2m finfile does not exist, so BST step will be run: ${bst_finfile_2m}"
        run_bst=true
    fi
fi

echo
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


# BST step
bst_failure=false
if [ "$run_bst" = true ]; then
    if [ ! -d "$tile_subtiles_dir" ]; then
        echo "Creating tile subtiles dir: ${tile_subtiles_dir}"
        mkdir -p "$tile_subtiles_dir"
    fi
    cmd="$bst_cmd"
    echo "Running BST process with command: ${cmd}"
    export ARG_TILENAME="$tileName"
    export ARG_WATERTILEDIR="$waterTileDir"
    eval "$cmd"
    cmd_status=$?
    if (( cmd_status != 2 )); then
        echo "BST jobscript exited with non-success(2) exit status (${cmd_status})"
        bst_failure=true
    fi
fi
if [ "$bst_failure" = true ]; then
    echo "Exiting before MST step due to BST failure"
    exit 1
fi


mst_failure=false

# MST 10m step
echo
if [ -f "$mst_finfile_10m" ]; then
    echo "MST 10m finfile already exists, so MST 10m step will be skipped: ${mst_finfile_10m}"
else
    cmd="${mst_cmd_10m}"
    echo "Running 10m MST process with command: ${cmd}"
    export ARG_TILENAME="$tileName"
    eval "$cmd"
    cmd_status=$?
    if (( cmd_status != 2 )); then
        echo "MST jobscript exited with non-success(2) exit status (${cmd_status})"
        mst_failure=true
    fi
fi


# MST 2m quads step
echo
if [ "$make_10m_only" = false ]; then
    cmd="${mst_cmd_2m}"
    for quad in "${quads_arr[@]}"; do
        quadTileName="${tileName}_${quad}"
        mst_finfile_2m="${mst_finfile_2m_template//<quadTileName>/${quadTileName}}"
        if [ -f "$mst_finfile_2m" ]; then
            echo "MST 2m finfile already exists for quadtile ${quadTileName}: ${mst_finfile_2m}"
            continue
        fi
        echo
        echo "Running 2m MST quadtile ${quadTileName} with command: ${cmd}"
        export ARG_TILENAME="$quadTileName"
        eval "$cmd"
        cmd_status=$?
        if (( cmd_status != 2 )); then
            echo "MST jobscript exited with non-success(2) exit status (${cmd_status})"
            mst_failure=true
        fi
    done
fi


if [ "$keep_subtiles" = false ]; then
    echo
    if [ "$mst_failure" = true ]; then
        echo "Will not remove tile subtiles directory due to MST failure"
        exit 1
    fi
    if [ -d "$tile_subtiles_dir" ]; then
        echo "Removing tile subtiles directory: ${tile_subtiles_dir}"
        rm -rf "$tile_subtiles_dir"
    else
        echo "Cannot remove tile subtiles directory, it does not exist: ${tile_subtiles_dir}"
    fi
fi
