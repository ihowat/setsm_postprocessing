import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing3'


def main():
    
    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")
    
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_addUnreg2Reg.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')
    
    args = parser.parse_args()
    
    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_addUnreg2Reg.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
           
    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))
    
    matlab_script = 'addUnreg2Reg'
    
    i=0
    if len(tiles) > 0:
    
        for tile in tiles:
            tiledir=os.path.join(dstdir,tile)
                 
            ## if pbs, submit to scheduler
            i+=1
            if args.pbs:
                job_name = 'addunreg{}'.format(tile)
                cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={} {}'.format(
                    job_name,
                    scriptdir,
                    matlab_script,
                    tiledir,
                    args.lib_path,
                    qsubpath
                )
                print cmd
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)
            
            ## else run matlab
            else:
                cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}'); exit" """.format(
                    scriptdir,
                    args.lib_path,
                    matlab_script,
                    tiledir,
                )
                print "{}, {}".format(i, cmd)
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()