#!/bin/bash

if (( $# != 1 )); then
    echo 1>/dev/stderr "Script requires one argument, the path to a geocell list text file"
    exit 1
fi
geotilelist="$1"

tileset_name=$(basename "$geotilelist" | cut -d"." -f1)
dstdir="results/gtp_deliveries/${tileset_name}/"
density_table="${dstdir}/matchtag_density_table.txt"

echo
echo "Will link out geotile results here: ${dstdir}"
echo
python ~/scratch/repos/pyscript-utils/file_xfer.py -sl "$geotilelist" -slp "results/gtp/" -d "$dstdir" -tt --fexcl "*_count_*.tif" "*_countmt_*.tif" "*.log"

echo
echo "Creating matchtag density table here: ${density_table}"
while IFS='' read -r dfile; do
    dfname=$(basename "$dfile")
    echo "$(echo "$dfname" | sed -r 's#(_count|_countmt|_density\.txt)##g'),$(cat "$dfile")"
done <<< "$(find "$dstdir" -type f -name "*density.txt" | sort)" > "$density_table"
find "$dstdir" -type f -name "*density.txt" -delete

echo "Geotile results have been linked with density table here: ${dstdir}"
