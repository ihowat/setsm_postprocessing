import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'

tileDefFile = 'PGC_Imagery_Mosaic_Tiles_Arctic.mat'
databaseFile = 'arcticdem_database_unf_pgcpaths.mat'
waterTileDir = '/mnt/pgc/data/scratch/claire/pgc/arcticdem/coastline/global_surface_water/tiles_v2/'
refDemFile = '/mnt/pgc/data/scratch/claire/pgc/arcticdem/mosaic/tandemx_alaska_mosaic_3413_tap90m.tif'

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")


    # parser.add_argument("region", choices=['arctic','antarctic','above'], help="region (arctic, antarctic, or above)")

    parser.add_argument("--tile-def", default=tileDefFile,
                        help="mosaic tile definition mat file(default={}".format(tileDefFile))
    parser.add_argument("--strip-db", default=databaseFile,
                        help="strip database mat file (default={}".format(databaseFile))
    parser.add_argument("--water-tile-dir", default=waterTileDir,
                        help="directory of water tifs (default={}".format(waterTileDir))
    parser.add_argument("--ref-dem", default=refDemFile,
                        help="reference DEM (default={}".format(refDemFile))
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))

    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--rerun", action='store_true', default=False,
            help="rerun tile, behavior determined by redoFlag in Matlab code")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_buildSubTiles.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_buildSubTiles.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    matlab_script = 'buildSubTiles'

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ## if output does not exist, add to task list
            tile_dstdir = os.path.join(dstdir,tile,'subtiles')
            if not os.path.isdir(tile_dstdir):
                if not args.dryrun:
                    os.makedirs(tile_dstdir)
            dstfps = glob.glob(os.path.join(tile_dstdir,'{}_*2m.mat'.format(tile)))

            if len(dstfps) > 0 and not args.rerun:
                print '{} subtiles exists, skipping'.format(tile)

            else:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'bst_{}'.format(tile)
                    cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8},p8={9},p9={10} {0}'.format(
                        qsubpath,
                        job_name,
                        scriptdir, #p1
                        args.lib_path,  #p2
                        matlab_script, #p3
                        tile, #p4
                        tile_dstdir, #p5
                        args.tile_def,  #p6
                        args.strip_db,  #p7
                        args.water_tile_dir,  #p8
                        args.ref_dem,  #p9
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}','{4}','{5}','{6}','{7}','{8}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        tile,
                        tile_dstdir,
                        args.tile_def,
                        args.strip_db,
                        args.water_tile_dir,
                        args.ref_dem,
                    )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()