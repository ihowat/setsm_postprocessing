#!/bin/bash

tilelist_string="$ARG_TILENAME"
tilerun_jobscript="$TILERUN_JOBSCRIPT"

IFS='@' read -r -a tilelist_arr <<< "$tilelist_string"
total_ntiles=${#tilelist_arr[@]}

for tile_idx in "${!tilelist_arr[@]}"; do
    tile="${tilelist_arr[$tile_idx]}"
    echo -e "\n\nRunning tile $((tile_idx+1)) of ${total_ntiles}: ${tile}"
    export ARG_TILENAME="$tile"
    bash "$tilerun_jobscript"
done
