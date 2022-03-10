import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing3')
default_qsub = 'qsub_tiles2tif.sh'

#### TODO add projstring to passed args
def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic tiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )
    parser.add_argument("res", choices=['2','8','40'], help="resolution (2, 8, or 40)")
    parser.add_argument("region", choices=['arctic','antarctic','above'], help="region (arctic, antarctic, or above)")

    parser.add_argument("--rerun", action='store_true', default=False,
            help="run script even if target dem already exists")
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(default_qsub))
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

    dstdir = os.path.abspath(args.dstdir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,default_qsub)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    if args.region == 'arctic':
        projstr = 'polar stereo north'
    elif args.region == 'antarctic':
        projstr = 'polar stereo south'
    elif args.region == 'above':
        projstr = 'canada albers equal area conic'

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ## if output does not exist, add to task list
            dstfp = os.path.join(dstdir,tile,'{}_{}m_reg_dem.tif'.format(tile, args.res))
            dstfp2 = os.path.join(dstdir,tile,'{}_{}m_dem.tif'.format(tile, args.res))
            ## reg_dem matfile will be substituted in the matlab script if exists
            matfile = os.path.join(dstdir,tile,'{}_{}m_dem.mat'.format(tile, args.res))
            if not os.path.isfile(matfile):
                print('source matfile does not exist: {}'.format(matfile))

            else:
                if (os.path.isfile(dstfp) or os.path.isfile(dstfp2)) and not args.rerun:
                    print('{} or {} exists, skipping'.format(dstfp, dstfp2))

                else:
                    ## if pbs, submit to scheduler
                    i+=1
                    if args.pbs:
                        job_name = 't2t_{}'.format(tile)
                        cmd = r'qsub -N {} -v p1={},p2={},p3="{}",p4={},p5={} {}'.format(
                            job_name,
                            scriptdir,
                            matfile,
                            projstr,
                            args.lib_path,
                            args.res,
                            qsubpath
                        )
                        print(cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)

                    ## else run matlab
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); writeTileToTif('{}',{},'{}'); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matfile,
                            args.res,
                            projstr
                        )
                        print("{}, {}".format(i, cmd))
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
