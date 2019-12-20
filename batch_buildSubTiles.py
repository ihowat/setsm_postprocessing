import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
default_dstdir = '/mnt/pgc/data/scratch/claire/pgc/arcticdem/mosaic/2m_v4'


def main():

    ## args
    parser = argparse.ArgumentParser()
    # parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")
    # parser.add_argument("res", choices=['2','8','20','40'], help="resolution (2, 8, or 40)")
    # parser.add_argument("region", choices=['arctic','antarctic','above'], help="region (arctic, antarctic, or above)")

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
    #dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_buildSubTiles.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    # if not os.path.isdir(dstdir):
    #     parser.error("dstdir does not exist: {}".format(dstdir))

    matlab_script = 'buildSubTiles'

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ## if output does not exist, add to task list
            dstfps = glob.glob(os.path.join(default_dstdir,tile,'subtiles','{}_*2m.mat'.format(tile)))

            if len(dstfps) > 0 and not args.rerun:
                print '{} subtiles exists, skipping'.format(tile)

            else:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'bst_{}'.format(tile)
                    cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={} {}'.format(
                        job_name,
                        scriptdir,
                        matlab_script,
                        tile,
                        args.lib_path,
                        qsubpath
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    #cmd = """matlab -nodisplay -nosplash -r "addpath('{}'); parpool(4); selectTileByName('{}',{}); exit" """.format(scriptdir, tile, args.res)
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        tile,
                    )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()