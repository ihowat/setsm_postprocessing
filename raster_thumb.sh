#!/bin/bash

## Bash settings
set -uo pipefail

if (( $# < 2 )); then
    echo "Usage: SRCFILE OUTFILE [RES_METERS] [RESAMPLE_METHOD]"
    exit 1
fi

set +u
srcfile="$1"
outfile="$2"
res_meters="$3"
resampling="$4"
set -u

regexgrp_imgaux_minimum='<MDI key="STATISTICS_MINIMUM">(-?[0-9]+(\.[0-9]*)?)</MDI>'
regexgrp_imgaux_maximum='<MDI key="STATISTICS_MAXIMUM">(-?[0-9]+(\.[0-9]*)?)</MDI>'

imgaux_file="${srcfile}.aux.xml"
echo "Getting image min/max stats from source raster: ${srcfile} -> ${imgaux_file}"
cmd="gdalinfo -stats ${srcfile}"
echo "$cmd"
eval "$cmd" 1>/dev/null
if [ ! -f "$imgaux_file" ]; then
    echo "Expected image aux file does not exist: ${imgaux_file}"
    exit 1
fi

imgaux_min=$(grep -Eoi "$regexgrp_imgaux_minimum" "$imgaux_file" | sed -r "s|${regexgrp_imgaux_minimum}|\1|")
imgaux_max=$(grep -Eoi "$regexgrp_imgaux_maximum" "$imgaux_file" | sed -r "s|${regexgrp_imgaux_maximum}|\1|")
if [ -z "$imgaux_min" ] || [ -z "$imgaux_max" ]; then
    echo "Unable to parse image min/max stats from aux file: ${imgaux_file}"
    exit 1
fi
echo "Removing image aux file: ${imgaux_file}"
rm "$imgaux_file"

cmd="gdal_translate"
if [ -n "$res_meters" ]; then
    cmd="${cmd} -tr ${res_meters} ${res_meters}"
fi
if [ -n "$resampling" ]; then
    cmd="${cmd} -r ${resampling}"
fi
cmd="${cmd} -ot Byte -scale ${imgaux_min} ${imgaux_max} 0 255"
cmd="${cmd} -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW"
cmd="${cmd} \"${srcfile}\" \"${outfile}\""
echo "Running main GDAL command:"
echo "$cmd"
eval "$cmd"
