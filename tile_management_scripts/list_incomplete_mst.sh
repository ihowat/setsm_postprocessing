#!/bin/bash

if (( $# < 1 )); then
    echo 1>/dev/stderr "Script usage: TILENAME [OUTPUT_TILES_DIR]"
    exit 1
fi
tilename="$1"
output_tiles_dir="${2:-./results/output_tiles/}"
logdir="$(dirname "$output_tiles_dir")/logs/$(basename "$output_tiles_dir")/matlab/mst/2m/"

if [ ! -d "$output_tiles_dir" ]; then
    echo 1>/dev/stderr "Output tiles directory does not exist: ${output_tiles_dir}"
    exit 1
fi
if [ ! -d "$logdir" ]; then
    echo 1>/dev/stderr "Output tiles logging directory does not exist: ${logdir}"
    exit 1
fi


logfiles=$(find "${logdir}/" -mindepth 1 -maxdepth 1 -type f -name "${tilename}*.log")
nmatfiles=$(find "${output_tiles_dir}/${tilename}/" -mindepth 1 -maxdepth 1 -type f \( -name "${tilename}_10m.mat" -o -name "${tilename}_?_?_2m.mat" \) | wc -l)


error_msg=''
nosubtilesfound_msg=''
if [ -n "$logfiles" ]; then
    nologs_msg=''
    while IFS= read -r logfile; do
        if grep -q -i 'error' "$logfile"; then
            error_msg=" -- ERROR"
            if grep -q -E 'ERROR: (No files found matching|No subtile matfiles found matching)' "$logfile"; then
                if grep -q 'ERROR (more of a warning): Found all 10000 subtile-is-empty indicator files' "$logfile"; then
                    nosubtilesfound_msg=" (all subtiles were empty, check for little/poor strip coverage of whole tile)"
                elif grep -q -E 'ERROR: Found only [0-9]+ subtile-is-empty indicator files matching' "$logfile"; then
                    nosubtilesfound_msg=" (subtile-building was incomplete, check for BST/workflow issues)"
                else
                    nosubtilesfound_msg=" (no subtiles found, check for BST/workflow issues)"
                fi
            fi
        fi
    done <<< "$logfiles"
else
    nologs_msg=" (no MST logs)"
fi


echo "${logdir}/${tilename}*.log -- ${nmatfiles}/5 matfiles${nologs_msg}${error_msg}${nosubtilesfound_msg}"
