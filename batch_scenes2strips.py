import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing'

def main():
    
    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("src", help="source directory")
    parser.add_argument("res", choices=['2','8'], help="resolution of target dems (2 or 8)")
    
    parser.add_argument("--no-entropy", action='store_true', default=False,
            help="use filter without entropy protection")
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
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_scenes2strips.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
    
    scene_dems = glob.glob(os.path.join(src,'*dem.tif'))
    
    ## find uniue strip IDs (<catid1_catid2>)    
    if len(scene_dems) > 0:
        stripids = list(set([os.path.basename(s)[14:47] for s in scene_dems]))
        stripids.sort()
        
        i=0
        for stripid in stripids:
            i+=1
            ## if dem date is less than existing strip date, remove existing output, else skip
            dstdir = src.replace('tif_results','strips')
            dst_dems = glob.glob(os.path.join(dstdir,'*'+stripid+'_seg*_dem.tif'))
            if len(dst_dems) > 0:
                src_datamasks = glob.glob(os.path.join(src,'*'+stripid+'*_matchtag.tif'))
                if len(src_datamasks) > 0:
                    a = min([os.path.getmtime(f) for f in dst_dems])
                    b = max([os.path.getmtime(f) for f in src_datamasks])
                    if a < b:
                        print '{} old strip exists, deleting and reprocessing'.format(stripid)
                        for fp in glob.glob(os.path.join(dstdir,'*'+stripid+'_seg*')):
                            #print 'removed {}'.format(fp)
                            if not args.dryrun:
                                os.remove(fp)
            
            if not os.path.isdir(dstdir):
                if not args.dryrun:
                    os.makedirs(dstdir)            
            
            ## if output does not exist, add to task list
            dst_dems = glob.glob(os.path.join(dstdir,'*'+stripid+'_seg*_dem.tif'))
            if len(dst_dems) > 0:
                print '{} output files exist, skipping'.format(stripid)
            
            else:
                
                ## if pbs, submit to scheduler
                if args.pbs:
                    if args.no_entropy:
                        script = 'scenes2strips_single_noentropy'
                    else:
                        script = 'scenes2strips_single'
                    job_name = 's2s{:04g}'.format(i)
                    cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={},p5={},p6={} {}'.format(
                        job_name,
                        scriptdir,
                        src,
                        stripid,
                        args.res,
                        script,
                        args.lib_path,
                        qsubpath
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)
                
                ## else run matlab
                else:
                    if args.no_entropy:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); scenes2strips_single_noentropy('{}','{}','{}'); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            src,
                            stripid,
                            args.res
                        )
                        print "{}, {}".format(i, cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); scenes2strips_single('{}','{}','{}'); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            src,
                            stripid,
                            args.res
                        )
                        print "{}, {}".format(i, cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)
                            


if __name__ == '__main__':
    main()