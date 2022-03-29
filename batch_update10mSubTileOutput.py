import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
qsubscript = 'qsub_update10mSubTileOutput.sh'


def main():
    
    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory where tiles exist")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic tiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )
    
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(qsubscript))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')
    
    args = parser.parse_args()

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            tiles = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        tiles = args.tiles.split(',')
    tiles = sorted(list(set(tiles)))

    dstdir = os.path.abspath(args.dstdir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,qsubscript)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
           
    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))
    
    matlab_script = 'call_update10mSubTileOutput'
    
    i=0
    if len(tiles) > 0:
    
        for tile in tiles:
            tiledir=os.path.join(dstdir,tile,'subtiles')
                 
            ## if pbs, submit to scheduler
            i+=1
            if args.pbs:
                job_name = 'upd10m{}'.format(tile)
                cmd = r'qsub -N {} -v p1={},p2={},p3={},p4={} {}'.format(
                    job_name,
                    scriptdir,
                    matlab_script,
                    tiledir,
                    args.lib_path,
                    qsubpath
                )
                print(cmd)
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)
            
            ## else run matlab
            else:
                cmd = """matlab -nojvm -nodisplay -nosplash -r "try; addpath('{0}'); addpath('{1}'); {2}('{3}'); catch e; disp(getReport(e)); exit(1); end; exit(0)" """.format(
                    scriptdir,
                    args.lib_path,
                    matlab_script,
                    tiledir,
                )
                print("{}, {}".format(i, cmd))
                if not args.dryrun:
                    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
