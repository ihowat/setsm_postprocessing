#!/bin/bash

#demdir='/data2/ArcticDEM/region_04_greenland_central/strips/8m'
demdir=$1
f=(`find $demdir -name '*dem.tif'`)
for i in ${f[@]}; do
	
	hfout=${i/dem.tif/dem_browse.tif}

	if ! [ -e $hfout ]; then
	
			echo "writing "$hfout
			
			tempname=temp.$$.tif			
			gdal_translate -tr 20 20 -r bilinear -a_nodata -9999 $i $tempname
			gdaldem hillshade -z 4 -compute_edges  -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW $tempname $hfout
			
			rm $tempname

	fi

	ofin=${i/dem.tif/ortho.tif}
	ofout=${ofin/ortho.tif/ortho_browse.tif}
	
	if ! [ -e $ofout ]; then 
		echo "writing "$ofout
		gdal_translate -ot Byte -scale -co compress=lzw -a_nodata 0 -tr 20 20 -r bilinear $ofin $ofout	
	fi
	
done
