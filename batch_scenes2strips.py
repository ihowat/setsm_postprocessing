import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing3'

def main():
    
    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("src", help="source directory")
    parser.add_argument("dst", help="destination directory")
    parser.add_argument("res", choices=['2','8'], help="resolution of target dems (2 or 8)")
    
    parser.add_argument("--rema2a", action='store_true', default=False,
            help="use filter rema2a")
    parser.add_argument("--noentropy", action='store_true', default=False,
            help="use no entropy filter")
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_scenes2strips.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.src):
        parser.error('src must be a directory')

    src = os.path.abspath(args.src)
    dstdir = os.path.abspath(args.dst)
    scriptdir = os.path.dirname(sys.argv[0])

    if args.noentropy and args.rema2a:
        parser.error('no entropy and rema2a filters are incompatible')

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_scenes2strips.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
    
    scene_dems = glob.glob(os.path.join(src,'*dem.tif'))
    
    if src == dstdir:
        parser.error('src dir is the same as the dst dir: {}'.format(dstdir))
    
    if not os.path.isdir(dstdir):
        if not args.dryrun:
            os.makedirs(dstdir)
            
    ## find uniue strip IDs (<catid1_catid2>)    
    if len(scene_dems) > 0:
        stripids = list(set([os.path.basename(s)[14:47] for s in scene_dems]))
        stripids.sort()
        
        i=0
        for stripid in stripids:
            i+=1            
                      
            ## if output does not exist, add to task list
            dst_dems = glob.glob(os.path.join(dstdir,'*'+stripid+'_seg*_dem.tif'))
            if len(dst_dems) > 0:
                print '{} output files exist, skipping'.format(stripid)
            
            else:
                if args.rema2a:
                    script = 'scenes2strips_single_rema2a'
                elif args.noentropy:
                    script = 'scenes2strips_single_noentropy'
                else:
                    script = 'scenes2strips_single'
                
                ## if pbs, submit to scheduler
                if args.pbs:
                    job_name = 's2s{:04g}'.format(i)
                    cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={},p5={},p6={},p7={} {}'.format(
                        job_name,
                        scriptdir,
                        src,
                        stripid,
                        args.res,
                        script,
                        args.lib_path,
                        dstdir,
                        qsubpath
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)
                
                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); {}('{}','{}','{}','{}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        script,
                        src,
                        stripid,
                        args.res,
                        dstdir
                    )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)
                            


if __name__ == '__main__':
    main()