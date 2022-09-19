#!/bin/bash

CORES_PER_NODE=$PBS_NUM_PPN
APRUN_PREFIX="aprun -b -N 1 -d ${CORES_PER_NODE} -cc none"
export ALREADY_IN_APRUN=true

tilelist_string="$ARG_TILENAME"
tilerun_jobscript="$TILERUN_JOBSCRIPT"
run_tiles_in_parallel="$IN_PARALLEL"

IFS='@' read -r -a tilelist_arr <<< "$tilelist_string"
total_ntiles=${#tilelist_arr[@]}

if [ "$run_tiles_in_parallel" = true ]; then
    echo "Will run ${total_ntiles} tiles in parallel"
else
    echo "Will run ${total_ntiles} tiles in serial"
fi

for tile_idx in "${!tilelist_arr[@]}"; do
    tile="${tilelist_arr[$tile_idx]}"
    echo -e "\n\nRunning tile $((tile_idx+1)) of ${total_ntiles}: ${tile}"
    export ARG_TILENAME="$tile"
    if [ "$run_tiles_in_parallel" = true ]; then
        ${APRUN_PREFIX} bash "$tilerun_jobscript" &
    else
        ${APRUN_PREFIX} bash "$tilerun_jobscript"
    fi
done

if [ "$run_tiles_in_parallel" = true ]; then
    wait
fi
