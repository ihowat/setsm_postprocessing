import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
quadnames = ('1_1','1_2','2_1','2_2')
RESOLUTIONS = ['2','10']
REGIONS = ['arctic','antarctic']
default_qsub = 'qsub_tiles2tif_v4.sh'

#### TODO add projstring to passed args
def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")
    parser.add_argument("res", choices=RESOLUTIONS, help="resolution ({})".format(','.join(RESOLUTIONS)))
    parser.add_argument("region", choices=REGIONS, help="region ({})".format(','.join(REGIONS)))

    parser.add_argument("--meta-only", action='store_true', default=False,
                        help="build meta files only")
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

    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,default_qsub)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    # TODO: add earthdem projection handling
    if args.region == 'arctic':
        projstr = 'polar stereo north'
    elif args.region == 'antarctic':
        projstr = 'polar stereo south'

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            for q in quadnames:
                tq = "{}_{}".format(tile,q)
                dstfp = os.path.join(dstdir,tile,'{}_{}m_dem.tif'.format(tq, args.res))
                matfile = os.path.join(dstdir,tile,'{}_{}m.mat'.format(tq, args.res))
                if not os.path.isfile(matfile):
                    print "Tile {} {}m mat file does not exist: {}".format(tq,args.res,matfile)

                elif os.path.isfile(dstfp) and not args.rerun:
                    print '{} exists, skipping'.format(dstfp)

                else:
                    ## if pbs, submit to scheduler
                    i+=1
                    if args.pbs:
                        job_name = 't2t_{}'.format(tq)
                        cmd = r'qsub -N {} -v p1={},p2={},p3="{}",p4={},p5={} {}'.format(
                            job_name,
                            scriptdir,
                            matfile,
                            projstr,
                            args.lib_path,
                            'true' if args.meta_only else 'false',
                            qsubpath,
                        )
                        print cmd
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)

                    ## else run matlab
                    else:
                        if args.meta_only:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                            )
                        else:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); writeTileToTifv4('{2}','{3}'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                                projstr
                            )
                        print "{}, {}".format(i, cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()