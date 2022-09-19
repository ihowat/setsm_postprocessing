
## Undo registration and get back to cleanest state of tile .mat files (as if they're fresh out of MST step)

# Remove all *.tif and *_meta.txt mosaic results files
# Remove all *_reg.mat mosaic tile files (we will be redoing registration)
find v2/results/output_tiles/ -mindepth 2 -maxdepth 2 -type f ! \( -name "*_10m.mat" -o -name "*_2m.mat" -o -name "*.fin" \) -delete

# Undo tile buffer merge on all existing 10m and 2m (non-reg) tile .mat files
# Crop 2m tiles, which had 1100m quad-internal buffer, to 100m buffer on all four edges
for tiledir in v2/results/output_tiles/*/ ; do qsub -v "ARG_TILEDIR=${tiledir}" ~/scratch/repos/setsm_postprocessing_pgc/qsub_cropTiles.sh ; done

# Rename cropped output .mat files, backing up uncropped
find v2/results/output_tiles/ -mindepth 2 -maxdepth 2 -type f -name "*_cropped.mat" -exec bash -c 'tilef_crop={}; tilef_orig="${tilef_crop%_cropped.mat}.mat"; mv "$tilef_orig" "${tilef_orig}.bak"; mv "$tilef_crop" "$tilef_orig";' \;



## Add land mask to tile .mat files (unregistered)
# Should only need to be done once for each new tile out of MST, could be added into MST step

# Add land mask to both 10m and 2m tile .mat files
for tiledir in v2/results/output_tiles/*/ ; do qsub -v "ARG_TILEDIR=${tiledir}" ~/scratch/repos/setsm_postprocessing_pgc/qsub_addLandMask2REMATiles.sh ; done



## Register tiles to ICESat-2

# Create *_reg.mat copies of the tile .mat files with registration to ICESat-2.
# - Registration info is calculated and stored in the "reg" variable of both registered and unregistrered .mat tile files.
# - *_reg.mat files are created by the "applyRegistration" function with the 3D vector translation from registration applied to all data arrays.
# - The "fit2is2" function adds quadratic surface "sf" (fit equation) and "dzfit" (mesh grid) variables of point GCP offsets from ICESat-2, and applies offset to the "z" array (DEM data).
for tiledir in v2/results/output_tiles/*/ ; do qsub -v "ARG_TILEDIR=${tiledir}" ~/scratch/repos/setsm_postprocessing_pgc/qsub_batchRegisterTiles.sh ; done
# Alternatively, use new python batch_registerTiles.py script:
python ~/scratch/repos/setsm_postprocessing_pgc/batch_registerTiles.py v2/results/output_tiles/ v2/tilelists/tiles_over_ice_shelves_to_skip_dzfit.txt 2 --pbs --skip-dzfit
python ~/scratch/repos/setsm_postprocessing_pgc/batch_registerTiles.py v2/results/output_tiles/ v2/tilelists/tiles_not_over_ice_shelves_to_apply_dzfit.txt 2 --pbs

# NOTE: The below 10-to-2 registration method did not turn out well.
#       We ended up using the 2m independent quad tile registration method above.
# Create *_reg.mat copies of the tile .mat files with registration to ICESat-2.
# - Registration info is calculated and stored in the "reg" variable of both registered and unregistrered .mat tile files.
# - *_reg.mat files are created by the "applyRegistration" function with the 3D vector translation from registration applied to all data arrays.
# - The "fit2is2" function adds quadratic surface "sf" (fit equation) and "dzfit" (mesh grid) variables of point GCP offsets from ICESat-2, and applies offset to the "z" array (DEM data).
python ~/scratch/repos/setsm_postprocessing_pgc/batch_registerTiles.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 10 --pbs
# Alternatively, use old qsub script method:
for tiledir in v2/results/output_tiles/*/ ; do qsub -v "ARG_TILEDIR=${tiledir},ARG_RESOLUTION=10m" ~/scratch/repos/setsm_postprocessing_pgc/qsub_batchRegisterTiles.sh ; done

# Apply registration from the 10m tiles to their 2m quad tiles:
for tiledir in v2/results/output_tiles/*/ ; do qsub -v "ARG_TILEDIR=${tiledir}" ~/scratch/repos/setsm_postprocessing_pgc/qsub_batchApplyRegTo2m.sh ; done



## Create tile neighbor index file, needed by subsequent steps

# Create 2m index file at: /mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/tile_index_files/output_tiles/tileNeighborIndex_2m.mat
qsub ~/scratch/repos/setsm_postprocessing_pgc/qsub_make_tileNeighborIndex.sh -v "ARG_TILE_ROOTDIR=/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles,ARG_RESOLUTION=2m"



## Boundary adjustment step

# The batch_boundaryAdjust.py script will spin off jobs that run these two matlab scripts in serial:
# batch_boundaryAdjustCalc: Calculate DEM offset between tile and all its neighbors, save as 'dz0' variable in tile .mat files
# batch_boundaryAdjustApply: Apply 'z = z - dz0' to tile .mat files ('dz0' calculated in the previous step)
python ~/scratch/repos/setsm_postprocessing_pgc/batch_boundaryAdjust.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 2 --pbs --process-group stripe-1
# wait until all jobs finish...
python ~/scratch/repos/setsm_postprocessing_pgc/batch_boundaryAdjust.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 2 --pbs --process-group stripe-2
# wait until all jobs finish...
python ~/scratch/repos/setsm_postprocessing_pgc/batch_boundaryAdjust.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 2 --pbs --process-group stripe-3



## Merge tile buffers

python ~/scratch/repos/setsm_postprocessing_pgc/batch_mergeTileBuffer.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 2 --pbs --process-group row
# wait until all jobs finish...
python ~/scratch/repos/setsm_postprocessing_pgc/batch_mergeTileBuffer.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt 2 --pbs --process-group column



## Export tile data array GeoTIFFs and meta.txt file

python ~/scratch/repos/setsm_postprocessing_pgc/batch_tiles2tif_v4.py v2/results/output_tiles/ v2/tilelists/rema_tiles_all_FINAL.txt rema 2 --pbs --process-by supertile-dir
