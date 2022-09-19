import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')

waterTileDir = '/mnt/pgc/data/scratch/claire/pgc/arcticdem/coastline/global_surface_water/tiles_v2/'
gcpFile = '/mnt/pgc/data/scratch/claire/pgc/arcticdem/gcp/icesat/mat/GLA14_rel634.mat'
qsub_default = 'qsub_registerTileVert.sh'

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source directory (dir above tile subfolders)")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic tiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )


    # parser.add_argument("region", choices=['arctic','antarctic','above'], help="region (arctic, antarctic, or above)")

    parser.add_argument("--water-tile-dir", default=waterTileDir,
                        help="directory of water tifs (default={}".format(waterTileDir))
    parser.add_argument("--gcp-file", default=gcpFile,
                        help="GCP file in csv or mat format (default={}".format(gcpFile))
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))

    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(qsub_default))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            tiles = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        tiles = args.tiles.split(',')
    tiles = sorted(list(set(tiles)))

    srcdir = os.path.abspath(args.srcdir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,qsub_default)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(srcdir):
        parser.error("srcdir does not exist: {}".format(srcdir))

    matlab_script = 'registerTileVert'

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ## if output does not exist, add to task list
            tile_srcdir = os.path.join(srcdir,tile)

            ## if pbs, submit to scheduler
            i+=1
            if args.pbs:
                job_name = 'rtv_{}'.format(tile)
                cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7} {0}'.format(
                    qsubpath,
                    job_name,
                    scriptdir, #p1
                    args.lib_path,  #p2
                    matlab_script, #p3
                    tile_srcdir, #p4
                    args.water_tile_dir,  #p5
                    args.gcp_file,  #p6
                )
                print(cmd)
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)

            ## else run matlab
            else:
                cmd = """matlab -nojvm -nodisplay -nosplash -r "try; addpath('{0}'); addpath('{1}'); {2}('{3}','{4}','{5}'); catch e; disp(getReport(e)); exit(1); end; exit(0)" """.format(
                    scriptdir,
                    args.lib_path,
                    matlab_script,
                    tile_srcdir,
                    args.water_tile_dir,
                    args.gcp_file,
                )

                print("{}, {}".format(i, cmd))
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
