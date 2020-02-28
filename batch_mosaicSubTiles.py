import os, string, sys, argparse, glob, subprocess
from collections import namedtuple
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
default_buffer = 100

MosaicScheme = namedtuple('MosaicOrigin', 'xorigin yorigin xsize ysize')
mosaic_schemes = {
    'arctic' :      MosaicScheme(-4000000, -4000000, 100000, 100000),
    'artarctic' :   MosaicScheme(-4000000, -4000000, 100000, 100000),
}

Task = namedtuple('Task', 't st x0 x1 y0 y1')

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source root dir (level above tile name dir)")
    parser.add_argument("tiles", help="list of tiles, comma delimited")
    parser.add_argument("res", choices=['2','10'], help="resolution (2 or 10)")
    parser.add_argument("region", choices=['arctic','antarctic'], help="region (arctic, antarctic)")


    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))

    parser.add_argument('--num-rows', type=int, default=1,
            help="number of subtile rows")
    parser.add_argument('--num-cols', type=int, default=1,
            help="number of subtile columns")
    parser.add_argument('--buffer', type=int, default=default_buffer,
            help="tile overlap buffer distance in meters (default={})".format(default_buffer))
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

    ## Get origin for region
    scheme = mosaic_schemes[args.region]

    tasks = []

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            trow,tcol = tile.split('_')
            txorigin = scheme.xorigin + scheme.xsize * (int(tcol) - 1)
            tyorigin = scheme.yorigin + scheme.ysize * (int(trow) - 1)
            dx = scheme.xsize / args.num_cols
            dy = scheme.ysize / args.num_rows

            ## If single row/col, proceed with no subtile ID in name
            if args.num_rows == 1 and args.num_cols == 1:
                subtile_name = ''
                x0 = txorigin - args.buffer
                y0 = tyorigin - args.buffer
                x1 = txorigin + dx + args.buffer
                y1 = tyorigin + dy + args.buffer
                tasks.append(Task(tile, subtile_name, x0, x1, y0, y1))

            else:
                for i in range(args.num_rows):
                    for j in range(args.num_cols):
                        subtile_name = '_{}_{}'.format(i+1,j+1)
                        y0 = tyorigin + i * dy - args.buffer
                        x0 = txorigin + j * dx - args.buffer
                        y1 = tyorigin + i * dy + dy + args.buffer
                        x1 = txorigin + j * dx + dx + args.buffer
                        tasks.append(Task(tile, subtile_name, x0, x1, y0, y1))

    print("{} tasks found".format(len(tasks)))

    if len(tasks) > 0:
        for task in tasks:
            dstfn = "{}{}_{}m.mat".format(task.t,task.st,args.res)
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
                    cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={},p5={},p6={},p7={},p8={},p9={},p10={} {}'.format(
                        job_name,
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        task.x0,task.x1,
                        task.y0,task.y1,
                        dstfp,
                        qsubpath
                    )
                    print cmd
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}',{4},{5},{6},{7},{8},'{9}'); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        task.x0,task.x1,
                        task.y0,task.y1,
                        dstfp,
                    )
                    print "{}, {}".format(i, cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()