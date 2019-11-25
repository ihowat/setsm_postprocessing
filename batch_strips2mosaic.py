import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing'


def main():
    
    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")
    parser.add_argument("res", choices=['2','8','20','40'], help="resolution (2, 8, or 40)")
    parser.add_argument("region", choices=['arctic','antarctic','above'], help="region (arctic, antarctic, or above)")
    
    parser.add_argument("--rebuild", action='store_true', default=False, help="rebuild DEM from 40m template. 40m version must already exist)")
    parser.add_argument("--gcpfile", help="csv file of GCP points (must have x, y, and z, no headers)")
    parser.add_argument("--rerun", action='store_true', default=False,
                        help="run script even if target dem already exists")
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_strips2mosaic.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')
    
    args = parser.parse_args()
    
    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_strips2mosaic.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
        
    if args.gcpfile:
        if not os.path.isfile(args.gcpfile):
            parser.error("gcpfile does not exist: {}".format(args.gcpfile))
    
    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))
    
    if args.rebuild: 
        if args.region == 'arctic':
            matlab_script = 'selectTileByNameRebuild'    
        elif args.region == 'antarctic':
            matlab_script = 'selectTileByNameRebuild_Antarctic'
        elif args.region == 'above':
            matlab_script = 'selectTileByNameRebuild_Above'
            
    else:
        if args.region == 'arctic':
            matlab_script = 'selectTileByName'    
        elif args.region == 'antarctic':
            matlab_script = 'selectTileByName_Antarctic'
        elif args.region == 'above':
            matlab_script = 'selectTileByName_Above'
    
    i=0
    if len(tiles) > 0:
    
        for tile in tiles:
                 
            ## if output does not exist, add to task list
            if args.rebuild:
                dstfp = os.path.join(dstdir,tile,'{}_{}m_reg_dem.mat'.format(tile, args.res))
                dstfp2 = os.path.join(dstdir,tile,'{}_{}m_dem.mat'.format(tile, args.res))
            else:
                dstfp = os.path.join(dstdir,tile,'{}_{}m_reg_dem.mat'.format(tile, args.res))
                dstfp2 = os.path.join(dstdir,tile,'{}_{}m_dem.mat'.format(tile, args.res))
                
            if (os.path.isfile(dstfp) or os.path.isfile(dstfp2)) and not args.rerun:
                print '{} or {} exists, skipping'.format(dstfp, dstfp2)
            
            else:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'mos_{}'.format(tile)
                    if args.gcpfile:
                        cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={},p5={},p6={},p7={} {}'.format(
                            job_name,
                            scriptdir,
                            matlab_script,
                            dstdir,
                            tile,
                            args.res,
                            args.lib_path,
                            args.gcpfile,
                            qsubpath
                        )
                    else:
                        cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={},p5={},p6={} {}'.format(
                            job_name,
                            scriptdir,
                            matlab_script,
                            dstdir,
                            tile,
                            args.res,
                            args.lib_path,
                            qsubpath
                        )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)
                
                ## else run matlab
                else:
                    #cmd = """matlab -nodisplay -nosplash -r "addpath('{}'); parpool(4); selectTileByName('{}',{}); exit" """.format(scriptdir, tile, args.res)
                    if args.gcpfile:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); addpath('{1}/intersections'); {2}('{3}','{4}',{5},'{6}'); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            dstdir,
                            tile,
                            args.res,
                            args.gcpfile
                        )
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); addpath('{1}/intersections'); {2}('{3}','{4}',{5}); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            dstdir,
                            tile,
                            args.res
                        )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
