#!/bin/bash

#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=4gb
#PBS -m n
#PBS -k oe
#PBS -j oe

module load gdal/2.1.3

csvfile="$ARG_IS2_CSVFILE"
shpfile="${csvfile/.csv/.shp}"

if [ ! -f "$csvfile" ]; then
    echo "ERROR! is2 CSV file does not exist: ${csvfile}"
    exit 1
fi

echo "Converting is2 CSV to shapefile: ${csvfile}"
ogr2ogr -overwrite -s_srs EPSG:3031 -t_srs EPSG:3031 -oo X_POSSIBLE_NAMES=x -oo Y_POSSIBLE_NAMES=y  -f "ESRI Shapefile" "$shpfile" "$csvfile"
