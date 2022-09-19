#!/bin/bash

if (( $# < 1 )); then
    echo 1>/dev/stderr "Script usage: TILELIST_FILE [OUTPUT_TILES_DIR]"
    exit 1
fi
tilelist_file="$1"
output_tiles_dir="${2:-./results/output_tiles/}"

if [ ! -d "$output_tiles_dir" ]; then
    echo 1>/dev/stderr "Output tiles directory does not exist: ${output_tiles_dir}"
    exit 1
fi

while IFS='' read -r tile || [ -n "$tile" ]; do
    if [ -z "$tile" ]; then continue; fi

    nfin=$(find "${output_tiles_dir}/${tile}/" -mindepth 1 -maxdepth 1 -type f -name "*.fin" | wc -l)

    if (( nfin != 6 )) && (( nfin != 7 )); then
        echo "$tile"
    fi
done < "$tilelist_file"
