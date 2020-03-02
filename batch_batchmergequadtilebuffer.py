import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
quadnames = ('1_1','1_2','2_1','2_2')
qsub_default = 'qsub_mergequadtilebuffer.sh'

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory")
    parser.add_argument("dimension", choices=['row','column'], help="dimension on which to group tiles for merging")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")

    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(qsub_default))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,qsub_default)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    # Test tiles exist
    existing_tiles = []
    for t in tiles:
        for q in quadnames:
            tq = "{}_{}".format(t,q)
            filename = "{}/{}/{}_2m.mat".format(dstdir,t,tq)
            if not os.path.isfile(filename):
                print "Tile {} 2m mat file does not exist: {}".format(tq,filename)
            else:
                existing_tiles.append(tq)

    #  group tiles by dimension
    groups = {}
    for quad in existing_tiles:
        row = quad[0:2]+'_'+quad[6:7]
        col = quad[3:5]+'_'+quad[8:9]

        if args.dimension == 'row':
            key = row
        else:
            key = col

        if key not in groups:
            groups[key] = [quad]
        else:
            groups[key].append(quad)

    i=0
    if len(groups) > 0:
        keys = groups.keys()
        keys.sort()

        for key in keys:

            print "Submitting tile group from {} {}".format(args.dimension,key)
            quads = groups[key]

            if len(quads) < 2:
                print "Tile group {} has only 1 member: {}. Skipping".format(key, quads)
            else:
                tile_str = ";".join(quads)

                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'merge_{}'.format(key)
                    cmd = r'qsub -N {} -v p1={},p2={},p3="{}",p4={} {}'.format(
                        job_name,
                        scriptdir,
                        dstdir,
                        tile_str,
                        args.lib_path,
                        qsubpath
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); batch_batchMergeQuadTileBuffer('{}',{{'{}'}}); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        dstdir,
                        tile_str.replace(";","','")
                    )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()