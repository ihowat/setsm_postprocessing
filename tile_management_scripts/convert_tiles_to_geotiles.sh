#!/bin/bash

if (( $# != 3 )); then
    echo 1>/dev/stderr "Script requires three arguments: project (arcticdem|earthdem|rema), path to mosaic (tileset) results, path to output geocell tilelist textfile"
    exit 1
fi
project="$1"
tile_results="$2"
geocell_tilelist="$3"

dstdir="results/gtp/"

echo
echo "Will output geotilesplus results here: ${dstdir}"
echo
python ~/scratch/repos/setsm-utils/tiles2geocell.py "$tile_results" "$dstdir" geotileplus "$geocell_tilelist" --egm /mnt/pgc/data/projects/nga/trex/egm08_gtx/egm08_25.gtx --pbs --project "$project"
