#!/bin/bash

if (( $# != 1 )); then
    echo 1>/dev/stderr "Script requires one argument, the path to a tilelist text file"
    exit 1
fi
tilelist="$1"

tileset_name=$(basename "$tilelist" | cut -d"." -f1)
dstdir="results/tileset_deliveries/${tileset_name}/"

echo
echo "Will link out tile results here: ${dstdir}"
echo
python ~/scratch/repos/pyscript-utils/file_xfer.py -sl "$tilelist" -slp "results/output_tiles/" --mindepth 1 --maxdepth 1 -d "$dstdir" -tt --fmatch "*.tif" "*.txt"
echo
echo "Tile results have been linked out here: ${dstdir}"
