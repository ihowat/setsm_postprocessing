## Set up the processing environment

### Clone repos

The processing requires two branches of the same repo to be cloned. 

```shell
cd /path/to/repos/dir
git clone -b 'earthdem-v1.1' git@github.com:ihowat/setsm_postprocessing.git setsm_postprocessing_pgc
git clone -b '4.0-pgc' git@github.com:ihowat/setsm_postprocessing.git setsm_postprocessing4
```

### Configure environment

The CLI uses environment variables to obtain the paths of the repos that were just cloned. These can be read from
an `.env` file. A `.env.template` file is provided in the `setsm_postprocessing_pgc` repo that specifies the 
expected variables.

```shell
cd setsm_postprocessing_pgc/rema-mosaic
cp .env.template .env
# Open .env in your favorite editor and specify the variables. Vim for example.
vim .env
```

If you would rather the `.env` file be located somewhere outside the repo, you can specify its location by setting
the `REMA_MOSAIC_ENV_FILE` to its path. If this environment variable is set, it will load from the pointed file.
Otherwise, it will look for `.env` in the package root.

### Create environment

Create the mamba environment `mosaic-production` using the environment file located in the `setsm_postprcessing_pgc`
repo. This environment will contain the `rema-mosaic` CLI that is used throughout the processing.
```shell
cd /path/to/repos/setsm_postprocessing_pgc

# If make is available:
make create_env

# Else:
mamba env create -n mosaic-production -f environment.mosaic-production-clis.yml
mamba run -n mosaic-production python -m pip install -e ./earthdem-mosaic
mamba run -n mosaic-production python -m pip install -e ./rema-mosaic
```

### Validate setup

Verify that the environment is created and configured correctly by running the `show-settings` command of the 
`rema-mosaic` CLI. This command will read the `.env` file and display the configured values if everything is
set up correctly.

```shell
mamba activate mosaic-production
rema-mosaic show-settings
# Expected output
# {
#   "SETSM_POSTPROCESSING_PYTHON_DIR": "/path/to/setsm_postprocessing_pgc",
#   "SETSM_POSTPROCESSING_MATLAB_DIR": "/path/to/setsm_postprocessing4"
# }
```

## Setup and run BST & MST

Start by navigating to the project directory and creating a subdirectory for the area of interest (AOI) that you will
be processing tiles for. Then move into the AOI directory.

```shell
cd /path/to/project
mkdir <AOI DIR>
cd <AOI DIR>
```
Create a list of tiles with one tile name per line. 

```shell
touch tilelist.txt
# Add tile names to the file with a text editor or whatever method you prefer
vim tilelist.txt
```

Set up the working directory and configuration for a given date range. 

For example, to run a single season mosaic using strips between July 1st, 2012 and June 30th, 2013 possible flag values 
would be `--dstdir 2012_2013`, `--datefilt-start 20120701`, and `--datefilt-end 20130630`. This will create a directory 
named after the value passed to`--dstdir`. Move into that directory.

```shell
# Preform a dryrun first. If the output is acceptable, remove the --dryrun flag and rerun the command
rema-mosaic init-bst --dstdir YYYY_YYYY --datefilt-start YYYYMMDD --datefilt-end YYYYMMDD --dryrun
cd YYYY_YYYY
```

Run the BST & MST stages using the `run-bst` command. The `--show-command` flag will display the generated script
invocation without running anything. If the output of `--show-command` matches your expectation, remove that flag to 
preform a dryrun. If no errors or unexpected warnings are displayed, remove the `--dryrun` flag and run the command 
again to submit the tasks to the cluster.

```shell
rema-mosaic run-bst --dstdir . --tiles ../tilelist.txt --dryrun --show-command
```

## Check BST & MST outputs

Once the slurm jobs created in the previous section have completed, use the `--dryrun` functionality of the `run-bst`
command to evaluate the completeness of the BST & MST stages.
```shell
rema-mosaic run-bst --dstdir . --tiles ../tilelist.txt --dryrun
```

Interpreting the output:
- Output: `Tile 24_15 seems complete (MST 10m finfile exists)`
  - Meaning: BST & MST stages are complete


- Output: `Tile seems complete with BST step (subtiles_10m.fin exists) ... but MST step is not complete, so tile will be run`
  - Meaning: BST complete. MST should be rerun


- Output: `Verifying tile 27_16 BST results before rerun`
  - Meaning: BST & MST incomplete. Likely no overlapping strips. Run the following `grep` to investigate further.

```shell
grep "overlapping strips match date filter" logs/slurm/bst/10m/*.o*
```

Rerun any tiles that have overlapping strips but are not complete by running the `run-bst` command again with the tiles
to rerun specified as comma separated values passed to the `--tiles` flag.

```shell
rema-mosaic run-bst --dstdir . --tiles <tile1>,<tile2>,... --dryrun
```

Once the rerun jobs are complete, repeat the steps above to access the completion of those tiles.

## Add land mask

Add the `land` field to the tiles. This field is used in the `export-products` command later to fill sea surface
heights while exporting to GeoTIFF.

```shell
rema-mosaic add-land-mask --dstdir . --tiles ../tilelist.txt --dryrun
```

## Merge tile buffers

### Make backup copies of the matfiles (OPTIONAL)

The BST & MST stages are the longest and this merge stage is not reversible. The following command will make backup 
copies of the `.mat` files found in the tile directories with the provided tag appended to the end.

```shell
rema-mosaic backup-matfiles --dstdir . --tiles ../tilelist.txt --tag bstmst --dryrun
```

### Create neighbor index

The merge buffers routine needs an index that identifies which tiles are adjacent to one another. Create one with the
`create-neighbor-index` command.

```shell
rema-mosaic create-neighbor-index --dstdir . --show-command
```

### Preform tile buffer merging

The following two commands merge tile buffers by row and then by column. Make sure that all slurm tasks for row merging
have completed before submitting the column merging tasks.

```shell
rema-mosaic merge-tile-buffers --dstdir . --tiles ../tilelist.txt --axis row --dryrun --show-command
rema-mosaic merge-tile-buffers --dstdir . --tiles ../tilelist.txt --axis column --dryrun --show-command
```

## Export products

Finally, export the matfiles to GeoTIFF. The `export-products` command will create the full 'output-set' for all
tiles.

```shell
rema-mosaic export-products --dstdir . --tiles ../tilelist.txt --dryrun --show-command
```

### Build VRTs (OPTIONAL)

Use the `build-vrts` command to mosaic all the GeoTIFF tiles in the working directory.

```shell
rema-mosaic build-vrts --dstdir . --dryrun
```