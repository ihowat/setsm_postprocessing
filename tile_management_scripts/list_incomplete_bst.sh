#!/bin/bash

if (( $# < 1 )); then
    echo 1>/dev/stderr "Script usage: TILENAME [OUTPUT_TILES_DIR]"
    exit 1
fi
tilename="$1"
output_tiles_dir="${2:-./results/output_tiles/}"
logdir="$(dirname "$output_tiles_dir")/logs/$(basename "$output_tiles_dir")/matlab/bst/2m/"

if [ ! -d "$output_tiles_dir" ]; then
    echo 1>/dev/stderr "Output tiles directory does not exist: ${output_tiles_dir}"
    exit 1
fi
if [ ! -d "$logdir" ]; then
    echo 1>/dev/stderr "Output tiles logging directory does not exist: ${logdir}"
    exit 1
fi


logfile="${logdir}/${tilename}.log"
nsubtiles=$(grep -c 'of 10000' "$logfile")


if grep -q -i 'error' "$logfile"; then
    error_msg=" -- ERROR"
else
    error_msg=''
fi

if grep -q 'no land in tile, returning' "$logfile"; then
    noland_msg=" (no land in tile)"
else
    noland_msg=''
fi
if grep -q 'No strips overlap tile' "$logfile"; then
    nostrips_msg=" (no strips overlap tile)"
else
    nostrips_msg=''
fi

if [ -z "$error_msg" ] && { [ -n "$noland_msg" ] || [ -n "$nostrips_msg" ]; }; then
    error_msg=" -- ERROR"
fi


echo "${logfile} -- subtiles processed: ${nsubtiles}${error_msg}${noland_msg}${nostrips_msg}"
