# EarthDEM Mosaic

An installable python package that contains postprocessing stages for creating EarthDEM v1.1 tiles.

## Set up the processing environment

Create the conda environment:
```shell
mamba env create --file environment.yml
mamba activate earthdem-mosaic
```

Create and configure the .env for the project
```shell
cp .env.template .env           # start from the template
vim .env                        # replace placeholders with real values
earthdem-mosaic show-settings   # verify the CLI reads the variables correctly
```
## Set the $UTM_ZONE environment variable

Many of the command in the `earthdem-mosaic` CLI take the UTM Zone name as an argument. To enable the copy-and-pasting
of the example commands in this document, the UTM Zone name is assumed to be set in the $UTM_ZONE environement variable.

```shell
# For UTM 10N, store the value utm10n
export UTM_ZONE=utm10n
# Check the stored value anytime by echoing the variable
echo $UTM_ZONE
```

## Initialize UTM Zone for processing

Create empty working directories for the UTM zone using the `create-working-dirs` command:

```shell
# Preform a dryrun
earthdem-mosaic create-working-dirs $UTM_ZONE --verbose --dryrun
# Create the working directories
earthdem-mosaic create-working-dirs $UTM_ZONE --verbose
```

Hardlink the `.mat` and `.fin` files from for the zone from `$EARTHDEM_MOSAIC_MATFILE_SOURCE_DIR` to the working directories:

```shell
# Preform a dryrun
earthdem-mosaic link-source-matfiles $UTM_ZONE --verbose --dryrun
# Preform the linking
earthdem-mosaic link-source-matfiles $UTM_ZONE --verbose
```

Create a line-delimited list of supertiles to process:
```shell
find ./00-matfiles/ -maxdepth 1 -type d -name "utm*" | sed "s|./00-matfiles/||" | sort > ./all_supertiles.txt
```

## Review coregistration offsets

Hardlink the matfiles and finfiles into the `10-coregistration-debug` stage directory:

```shell
# Preform a dryrun of linking the matfiles
earthdem-mosaic link-files-to-stage --src ./00-matfiles --dst ./10-coregistration-debug --src-suffix ".mat" --verbose --dryrun
# Link the matfiles
earthdem-mosaic link-files-to-stage --src ./00-matfiles --dst ./10-coregistration-debug --src-suffix ".mat" --verbose

# Preform a dryrun of linking the fin files
earthdem-mosaic link-files-to-stage --src ./00-matfiles --dst ./10-coregistration-debug --src-suffix ".fin" --verbose --dryrun
# Link the fin files
earthdem-mosaic link-files-to-stage --src ./00-matfiles --dst ./10-coregistration-debug --src-suffix ".fin" --verbose
```

Submit the coregistraion debug process to the cluster:

```shell
# Inspect the generated command
earthdem-mosaic coreg-debug $UTM_ZONE ./all_supertiles.txt --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic coreg-debug $UTM_ZONE ./all_supertiles.txt --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic coreg-debug $UTM_ZONE ./all_supertiles.txt --slurm
```

Once complete, build a VRT mosaic of the `_offset.tif` files for review:

```shell
# Run from within the 10-coregistration-debug directory
gdalbuildvrt ./10-coregistration-debug_offset.vrt $(find -type f -name "*_offset.tif" | paste -sd " ")
```

Create a polygon shapefile named `<UTM ZONE>_skipreg.shp` in the zone working directory with an integer field called `skipreg`. 
Add features to this shapefile that overlap blobs that should not have coregistration applied. Set the `skipreg` value
to 1 for these features.

## Coregister the matfiles

Apply the coregistration to the matfiles, excluding the blobs intersecting the `skipreg` features, to produce `_reg.mat` files:

```shell
# Inspect the generated command
earthdem-mosaic coreg-matfiles $UTM_ZONE ./all_supertiles.txt --skipreg-shp "${UTM_ZONE}_skipreg.shp" --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic coreg-matfiles $UTM_ZONE ./all_supertiles.txt --skipreg-shp "${UTM_ZONE}_skipreg.shp" --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic coreg-matfiles $UTM_ZONE ./all_supertiles.txt --skipreg-shp "${UTM_ZONE}_skipreg.shp" --slurm
```

## Water flatten the matfiles

Apply water flattening to the `_reg.mat` files to produce the `_rem_fill.mat` files:

```shell
# Inspect the generated command
earthdem-mosaic water-flatten-matfiles $UTM_ZONE ./all_supertiles.txt --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic water-flatten-matfiles $UTM_ZONE ./all_supertiles.txt --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic water-flatten-matfiles $UTM_ZONE ./all_supertiles.txt --slurm
```

## Merge the overlapping matfile buffers

This next stage mutates the matfiles it operates on, rather than producing a separate copy. So first we will make a copy
of the `_reg_fill.mat` files with the suffix `_reg_fill_merge.mat` and then operate on those.

```shell
# Preform a dryrun
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 00-matfiles/ --src-suffix _reg_fill.mat --dst-suffix _reg_fill_merge.mat --verbose --dryrun
    
# Copy the matfiles
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 00-matfiles/ --src-suffix _reg_fill.mat --dst-suffix _reg_fill_merge.mat --verbose
```

Build the tile neighbor index. This will create a file at <UTM ZONE>/tile_index_files/tileNeighborIndex_2m.mat that is
used as a lookup table during the merge steps.

```shell
# Inspect the generated command
earthdem-mosaic create-neighbor-index $UTM_ZONE --show-command
# Run on the head node
earthdem-mosaic create-neighbor-index $UTM_ZONE
```

Merge tiles to their east and west neighbors (process by row):

```shell
# Inspect the generated command
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt row --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt row --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt row --slurm
```
Merge tiles to their north and south neighbors (process by column):

```shell
# Inspect the generated command
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt column --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt column --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt column --slurm
```

## Create final products

Link the `_reg_fill_merge.mat` file to the stage directory (`20-no-slope-filter` or `30-yes-slope-filter`), renaming the 
destination files to exclude the `_reg_fill_merge` part of the file name. Link the `.fin` files from `00-matfiles` to 
the destination folders as well. Then submit the TIF creation process to the cluster.

### No slope filter

```shell
# Link files to the 20-no-slope-filter directory
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 20-no-slope-filter/ --src-suffix _reg_fill_merge.mat --dst-suffix .mat
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 20-no-slope-filter/ --src-suffix .fin

# Inspect the generated command
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --slurm
```

### Yes slope filter

```shell
# Link files to the 30-yes-slope-filter directory
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 30-yes-slope-filter/ --src-suffix _reg_fill_merge.mat --dst-suffix .mat
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 30-yes-slope-filter/ --src-suffix .fin

# Inspect the generated command
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --apply-slope-filter --slurm --dryrun --show-command
# Preform a dryrun
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --apply-slope-filter --slurm --dryrun
# Submit the jobs to the cluster
earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --apply-slope-filter --slurm
```

## Review slope filter

Create VRT mosaics of the `*_browse.tif` files:

```shell
# Run from within the 20-no-slope-filter directory
gdalbuildvrt ./20-no-slope-filter_browse.vrt $(find -type f -name "*_browse.tif" | paste -sd " ")

# Run from within the 30-yes-slope-filter directory
gdalbuildvrt ./30-yes-slope-filter_browse.vrt $(find -type f -name "*_browse.tif" | paste -sd " ")
```

Create the slope filter review geopackage:

```shell
earthdem-mosaic slope-filter-review $UTM_ZONE --verbose
```

## Link final products to destination directory

Link the selected final products to the `$EARTHDEM_MOSAIC_FINAL_PRODUCTS_DIR`:

```shell
# Preform a dryrun
earthdem-mosaic link-final-products $UTM_ZONE "${UTM_ZONE}_slope_filter_review.gpkg" "${UTM_ZONE}_skipreg.shp" --verbose --dryrun
# Link the files
earthdem-mosaic link-final-products $UTM_ZONE "${UTM_ZONE}_slope_filter_review.gpkg" "${UTM_ZONE}_skipreg.shp" --verbose
```
