import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
quadname_list = ['1_1', '1_2', '2_1', '2_2']
RESOLUTIONS = ['2', '10']
REGIONS = ['arcticdem', 'earthdem', 'rema']
TIF_OUTPUT_CHOICES = ['browse-LZW', 'browse-COG', 'full-LZW', 'full-COG']
default_qsub = 'qsub_tiles2tif_v4.sh'


class RawTextArgumentDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass

def main():
    global quadname_list

    ## args
    parser = argparse.ArgumentParser(
        formatter_class=RawTextArgumentDefaultsHelpFormatter,
        description=""
    )
    parser.add_argument("tiledir", help="tile root directory containing supertile subfolders")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic supertiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )
    parser.add_argument("region", choices=REGIONS,
            help="processing domain of tiles")
    parser.add_argument("res", choices=RESOLUTIONS,
            help="resolution of source tiles to be exported")

    parser.add_argument("--tif-output", choices=TIF_OUTPUT_CHOICES, default='full-COG',
            help="type of tifs to create for raster output")
    parser.add_argument("--meta-only", action='store_true', default=False,
            help="build meta files only")
    parser.add_argument("--rerun", action='store_true', default=False,
            help="run script even if target dem already exists")
    parser.add_argument("--keep-old-results", action='store_true', default=False,
            help="do not remove existing results before submitting jobs")
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--hold", action='store_true', default=False,
            help="when submitting PBS jobs, submit them with held (H) status")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(default_qsub))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'supertile_list' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            supertile_list = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        supertile_list = args.tiles.split(',')
    supertile_list = sorted(list(set(supertile_list)))

    root_tiledir = os.path.abspath(args.tiledir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,default_qsub)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(root_tiledir):
        parser.error("tiledir does not exist: {}".format(root_tiledir))

    if args.region == 'arcticdem':
        projstr = 'polar stereo north'
    elif args.region == 'earthdem':
        projstr = None
    elif args.region == 'rema':
        projstr = 'polar stereo south'
    else:
        parser.error("unexpected region name")

    if args.res == '10':
        quadname_list = ['']

    tif_output_is_browse = args.tif_output in ('browse-LZW', 'browse-COG')

    i=0
    if len(supertile_list) > 0:

        for supertile in supertile_list:

            tile_projstr = projstr

            if tile_projstr is None:
                assert args.region == 'earthdem'

                utm_tilename_parts = supertile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(supertile))

                tile_projstr = utm_tilename_prefix

            for quadname in quadname_list:
                tile_name = '{}_{}'.format(supertile, quadname) if quadname != '' else supertile
                tile_rootpath = os.path.join(root_tiledir, supertile, '{}_{}m'.format(tile_name, args.res))

                unregmatfile    = '{}.mat'.format(tile_rootpath)
                regmatfile      = '{}_reg.mat'.format(tile_rootpath)
                finfp           = '{}.fin'.format(tile_rootpath)
                demfp           = '{}_dem.tif'.format(tile_rootpath)
                browsefp        = '{}_browse.tif'.format(tile_rootpath)
                metafp          = '{}_meta.txt'.format(tile_rootpath)
                matfile         = regmatfile if os.path.isfile(regmatfile) else unregmatfile

                run_tile = True

                if not os.path.isfile(unregmatfile) and not os.path.isfile(regmatfile):
                    print(
                        "Tile {} {}m mat and reg.mat files do not exist{}: {}".format(
                            tile_name, args.res,
                            " (AND .fin file also does not exist!!)" if not os.path.isfile(finfp) else '',
                            matfile
                        )
                    )
                    run_tile = False

                elif args.meta_only:
                    if os.path.isfile(metafp):
                        if args.rerun:
                            print("Removing existing meta file: {}".format(metafp))
                            if not args.dryrun:
                                os.remove(metafp)
                        else:
                            print("{} exists, skipping".format(metafp))
                            run_tile = False

                else:
                    if args.rerun:
                        assume_complete = False
                    if tif_output_is_browse:
                        assume_complete = os.path.isfile(browsefp) and os.path.isfile(metafp)
                    else:
                        assume_complete = os.path.isfile(demfp) and os.path.isfile(browsefp) and os.path.isfile(metafp)

                    if args.rerun or not assume_complete:
                        if not args.keep_old_results:
                            dstfps_old_pattern = [
                                demfp.replace('_dem.tif', '*.tif'),
                                metafp
                            ]
                            dstfps_old = [fp for pat in dstfps_old_pattern for fp in glob.glob(pat)]
                            if dstfps_old:
                                print("{}Removing existing tif tile results matching {}".format('(dryrun) ' if args.dryrun else '', dstfps_old_pattern))
                                if not args.dryrun:
                                    for dstfp_old in dstfps_old:
                                        os.remove(dstfp_old)

                    elif assume_complete:
                        "{} {} exist, skipping".format(
                            metafp,
                            "browse and meta" if tif_output_is_browse else "dem, browse, and meta"
                        )
                        run_tile = False

                if run_tile:
                    ## if pbs, submit to scheduler
                    i+=1
                    if args.pbs:
                        job_name = 't2t_{}'.format(tile_name)
                        cmd = r'qsub {} -N {} -v p1={},p2={},p3="{}",p4={},p5={},p6={} {}'.format(
                            '-h' if args.hold else '',
                            job_name,
                            scriptdir,
                            args.lib_path,
                            'true' if args.meta_only else 'false',
                            matfile,
                            tile_projstr,
                            args.tif_output,
                            qsubpath,
                        )
                        print(cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)

                    ## else run matlab
                    else:
                        if args.meta_only:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                            )
                        else:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); writeTileToTifv4('{2}','{3}','outRasterType','{4}'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                                tile_projstr,
                                args.tif_output,
                            )
                        print("{}, {}".format(i, cmd))
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
