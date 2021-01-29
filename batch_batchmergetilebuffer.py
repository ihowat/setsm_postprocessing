import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing3'

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
            help="qsub script to use in PBS submission (default is qsub_mergetilebuffer.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')
    
    args = parser.parse_args()
    
    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.dirname(sys.argv[0])

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_mergetilebuffer.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)
        
    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    # Test tiles exist
    existing_tiles = []
    for t in tiles:
        filename = "{0}/{1}/{1}_2m_reg_dem.mat".format(dstdir,t)
        if not os.path.isfile(filename):
            print("Tile {} 2m dem does not exist: {}".format(t,filename))
        else:
            existing_tiles.append(t)

        dstfps_old_pattern = [
            "{0}/{1}/{1}*2m*.tif".format(dstdir,t),
            "{0}/{1}/{1}*2m*meta.txt".format(dstdir,t)
        ]
        dstfps_old = [fp for pat in dstfps_old_pattern for fp in glob.glob(pat)]
        if dstfps_old:
            if not os.path.isfile(filename):
                print("ERROR! Tile mat file does not exist, but other MST results exist matching {}".format(dstfps_old_pattern))
                continue
            print("{}Removing old MST results matching {}".format('(dryrun) ' if args.dryrun else '', dstfps_old_pattern))
            if not args.dryrun:
                for dstfp_old in dstfps_old:
                    os.remove(dstfp_old)

    #  group tiles by dimension
    groups = {}
    for tile in existing_tiles:
        row = tile[:2]
        col = tile[3:]
        
        if args.dimension == 'row':
            key = row
        else:
            key = col
            
        if key not in groups:
            groups[key] = [tile]
        else:
            groups[key].append(tile)
        
    i=0
    if len(groups) > 0:
        keys = list(groups.keys())
        keys.sort()
        
        for key in keys:
            
            print("Submitting tile group from {} {}".format(args.dimension,key))
            tiles = groups[key]
            
            if len(tiles) < 2:
                print("Tile group {} has only 1 member: {}. Skipping".format(key, tiles))
            else:
                tile_str = ";".join(tiles)
                    
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
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)
                
                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{}'); addpath('{}'); batch_batchMergeTileBuffer('{}',{{'{}'}}); exit" """.format(
                        scriptdir,
                        args.lib_path,
                        dstdir,
                        tile_str.replace(";","','")
                    )
                    print("{}, {}".format(i, cmd))
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
