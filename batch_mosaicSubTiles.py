import os, string, sys, argparse, glob, subprocess
from collections import namedtuple
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
quads = ['1_1','1_2','2_1','2_2']

Task = namedtuple('Task', 't st')

tileDefFile = 'PGC_Imagery_Mosaic_Tiles_Arctic.mat'

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source root dir (level above tile name dir)")
    parser.add_argument("tiles", help="list of tiles, comma delimited")
    parser.add_argument("res", choices=['2','10'], help="resolution (2 or 10)")

    
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--tile-def", default=tileDefFile,
            help="mosaic tile definition mat file(default={}".format(tileDefFile))
    parser.add_argument('--quads', action='store_true', default=False,
            help="build into quad subtiles")
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_mosaicSubTiles.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    tiles = args.tiles.split(',')
    srcdir = os.path.abspath(args.srcdir)
    scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_mosaicSubTiles.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    # if not os.path.isdir(dstdir):
    #     parser.error("dstdir does not exist: {}".format(dstdir))

    matlab_script = 'mosaicSubTiles'

    tasks = []

    i=0
    if len(tiles) > 0:

        for tile in tiles:
            if args.quads:
                for quad in quads:
                    tasks.append(Task(tile, quad))
            else:
                tasks.append(Task(tile, ''))


    print("{} tasks found".format(len(tasks)))

    if len(tasks) > 0:
        for task in tasks:
            if task.st == '':
                dstfn = "{}_{}m.mat".format(task.t,args.res)
            else:
                dstfn = "{}_{}_{}m.mat".format(task.t,task.st,args.res)
            dstfp = os.path.join(srcdir, task.t, dstfn)
            sem = os.path.join(srcdir, task.t, dstfn.replace('.mat','_empty.txt'))
            subtile_dir = os.path.join(srcdir,task.t,'subtiles')

            if os.path.isfile(dstfp):
                print 'Output exists, skipping {}'.format(dstfn)
            elif os.path.isfile(sem):
                print 'N array was empty on last run, skipping {}'.format(dstfn)

            else:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'mst_{}'.format(task.t)
                    cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8},p8={9},p9={10} {0}'.format(
                        qsubpath,
                        job_name,
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        dstfp,
                        task.t,
                        args.tile_def,
                        task.st,
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    if task.st == '':
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{6}','{7}'); {2}('{3}',{4},'{5}','extent',[x0,x1,y0,y1]); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.t,
                            args.tile_def
                        )
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{7}','{8}','quadrant','{6}'); {2}('{3}',{4},'{5}','quadrant','{6}','extent',[x0,x1,y0,y1]); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.st,
                            task.t,
                            args.tile_def
                        )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()