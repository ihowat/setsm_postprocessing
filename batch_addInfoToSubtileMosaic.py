import os, string, sys, argparse, glob, subprocess
from collections import namedtuple

matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
quadnames = ('1_1','1_2','2_1','2_2')
qsub_default = 'qsub_addInfoToSubtileMosaic.sh'

Task = namedtuple('Task', 't st')

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source directory")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")
    parser.add_argument("res", choices=['2','10'], help="resolution (2 or 10)")

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
    srcdir = os.path.abspath(args.srcdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,qsub_default)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(srcdir):
        parser.error("srcdir does not exist: {}".format(srcdir))

    matlab_script = 'addInfoToSubtileMosaic'
    tasks = []

    for t in tiles:
        if args.res == '2':
            for q in quadnames:
                tasks.append(Task(t, q))

        else:
            tasks.append(Task(t, ''))

    i=0
    if len(tasks) > 0:
        for task in tasks:
            if task.st == '':
                dstfn = "{}_{}m.mat".format(task.t,args.res)
            else:
                dstfn = "{}_{}_{}m.mat".format(task.t,task.st,args.res)
            dstfp = os.path.join(srcdir, task.t, dstfn)
            subtile_dir = os.path.join(srcdir,task.t,'subtiles')

            ## if pbs, submit to scheduler
            i+=1
            if args.pbs:
                job_name = 'addinfo_{}'.format(task.t)
                cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8} {0}'.format(
                    qsubpath,
                    job_name,
                    scriptdir,
                    args.lib_path,
                    matlab_script,
                    subtile_dir,
                    args.res,
                    dstfp,
                    task.st
                )
                print(cmd)
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)

            ## else run matlab
            else:
                if task.st == '':
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}',{4},'{5}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        dstfp,
                    )
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}',{4},'{5}','{6}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        dstfp,
                        task.st,
                    )
                print("{}, {}".format(i, cmd))
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)



if __name__ == '__main__':
    main()
