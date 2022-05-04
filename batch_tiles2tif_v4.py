
import argparse
import glob
import os
import subprocess
import sys

from lib import jobscript_utils
from lib.jobscript_utils import wrap_multiline_str
from lib.jobscript_utils import UnsupportedMethodError


# This script paths
SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

# Paths relative to this script
MATLAB_LIBDIR = os.path.join(SCRIPT_DIR, "../setsm_postprocessing4")


class RawTextArgumentDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass

def get_arg_parser():
    
    parser = argparse.ArgumentParser(
        formatter_class=RawTextArgumentDefaultsHelpFormatter,
        description=wrap_multiline_str("""
            Export raster data from Matlab '<tilename>_<resolution>.mat'
            DEM mosaic tile files as GeoTIFF raster images, along with
            associated processing metadata in text format.
        """)
    )
    
    ## Positional args
    
    parser.add_argument(
        'tiledir',
        type=str,
        help="Root directory of tiles to be processed."
    )
    parser.add_argument(
        'tiles',
        help=wrap_multiline_str("""
            List of mosaic supertiles (10m tile names) to process,
            either specified on command line (comma delimited)
            or a text file list (each tile on separate line).
        """)
    )
    parser.add_argument(
        'domain',
        type=str,
        choices=['arcticdem', 'earthdem', 'rema'],
        help="DEM production domain of source tiles"
    )
    parser.add_argument(
        'resolution',
        type=int,
        choices=[10, 2],
        help=wrap_multiline_str("""
            Resolution class of existing tiles in source tiledir
            to be processed (meters).
        """)
    )

    ## Optional args

    parser.add_argument(
        '--tif-output',
        type=str,
        choices=['browse-LZW', 'browse-COG', 'full-LZW', 'full-COG'],
        default='full-COG',
        help="Package of raster data to output."
    )
    parser.add_argument(
        '--meta-only',
        action='store_true',
        help="Export tile metadata only."
    )
    
    parser.add_argument(
        '--rerun',
        action='store_true',
        help=wrap_multiline_str("""
            Submit processing jobs even if exported results files
            already exist.
        """)
    )
    parser.add_argument(
        '--keep-old-results',
        action='store_true',
        help=wrap_multiline_str("""
            Do not remove existing exported results files before submitting
            processing jobs.
        """)
    )

    parser.add_argument(
        '--tile-org',
        type=str,
        choices=['pgc', 'osu'],
        default='pgc',
        help="Tile file organization scheme."
    )
    parser.add_argument(
        '--process-by',
        type=str,
        choices=['supertile-dir', 'tile-file'],
        default='tile-file',
        help=wrap_multiline_str(r"""
            \nIf 'supertile-dir', if one tile within the supertile folder
            needs to be processed, hand whole supertile folder to Matlab script
            to be processed as a single task.
            \nIf 'tile-file', send each tile to be processed individually
            to Matlab script as a single task.
        """)
    )
    
    parser.add_argument(
        '--matlib',
        type=str,
        default=MATLAB_LIBDIR,
        help=wrap_multiline_str("""
            Path to directory containing necessary Matlab scripts
            and functions.
        """)
    )
    
    parser.add_argument(
        '--dryrun',
        action='store_true',
        help="Print actions without executing."
    )
    
    ## Argument groups
    
    jobscript_utils.argparse_add_job_scheduler_group(parser, config_group=SCRIPT_FNAME)

    return parser


def main():

    ## Parse and adjust arguments

    arg_parser = get_arg_parser()
    script_args = arg_parser.parse_args()

    root_tiledir = os.path.abspath(script_args.tiledir)

    if script_args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(script_args.tiles):
        tilelist_file = script_args.tiles
        if not os.path.isfile(script_args.tiles):
            arg_parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            supertile_list = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        supertile_list = script_args.tiles.split(',')
    supertile_list = sorted(list(set(supertile_list)))

    if script_args.domain == 'arcticdem':
        global_projstr = 'polar stereo north'
    elif script_args.domain == 'earthdem':
        global_projstr = None
    elif script_args.domain == 'rema':
        global_projstr = 'polar stereo south'
    else:
        raise UnsupportedMethodError(
            "No projstr mapping for argument 'domain': {}".format(script_args.domain)
        )

    res_name = '{}m'.format(script_args.resolution)

    if script_args.resolution == 2:
        quadname_list = ['1_1', '1_2', '2_1', '2_2']
    else:
        quadname_list = ['']

    tif_output_is_browse = script_args.tif_output in ('browse-LZW', 'browse-COG')

    batch_t2t_metaonly_args = ", 'metaOnly'" if script_args.meta_only else ''

    jobscript_utils.adjust_args(script_args, arg_parser)
    jobscript_utils.create_dirs(script_args, arg_parser)


    ## Verify arguments

    if not os.path.isdir(root_tiledir):
        arg_parser.error("Argument 'tiledir' is not an existing directory: {}".format(root_tiledir))

    if script_args.tile_org == 'osu' and script_args.process_by == 'supertile-dir':
        arg_parser.error("--process-by must be set to to 'tile-file' when --tile-org='osu'")

    jobscript_utils.verify_args(script_args, arg_parser)


    ## Main processing loop

    num_supertiles = len(supertile_list)
    if script_args.resolution == 2 and script_args.process_by == 'tile-file':
        est_num_tasks = num_supertiles * 4
    else:
        est_num_tasks = num_supertiles

    job_handler = jobscript_utils.JobHandler(
        script_args,
        est_num_tasks,
        num_tasks_is_estimate=True,
        init_env_requests='matlab gdal'
    )

    run_tile_list_all = []

    for supertile in supertile_list:

        tile_projstr = global_projstr
        if tile_projstr is None:
            assert script_args.domain == 'earthdem'

            utm_tilename_parts = supertile.split('_')
            utm_tilename_prefix = utm_tilename_parts[0]
            if not utm_tilename_prefix.startswith('utm'):
                arg_parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(supertile))

            tile_projstr = utm_tilename_prefix

        run_tile_matlist = []

        for quadname in quadname_list:
            tile_name = '{}_{}'.format(supertile, quadname) if quadname != '' else supertile
            tile_rootpath = os.path.join(
                root_tiledir,
                supertile if script_args.tile_org == 'pgc' else '',
                '{}_{}'.format(tile_name, res_name)
            )

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
                        tile_name, script_args.resolution,
                        " (AND .fin file also does not exist!!)" if not os.path.isfile(finfp) else '',
                        matfile
                    )
                )
                run_tile = False

            elif script_args.meta_only:
                if os.path.isfile(metafp):
                    if script_args.rerun:
                        print("Removing existing meta file: {}".format(metafp))
                        if not script_args.dryrun:
                            os.remove(metafp)
                    else:
                        print("{} exists, skipping".format(metafp))
                        run_tile = False

            else:
                if script_args.rerun:
                    assume_complete = False
                if tif_output_is_browse:
                    assume_complete = os.path.isfile(browsefp) and os.path.isfile(metafp)
                else:
                    assume_complete = os.path.isfile(demfp) and os.path.isfile(browsefp) and os.path.isfile(metafp)

                if script_args.rerun or not assume_complete:
                    if not script_args.keep_old_results:
                        dstfps_old_pattern = [
                            demfp.replace('_dem.tif', '*.tif'),
                            metafp
                        ]
                        dstfps_old = [fp for pat in dstfps_old_pattern for fp in glob.glob(pat)]
                        if dstfps_old:
                            print("{}Removing existing tif tile results matching {}".format('(dryrun) ' if script_args.dryrun else '', dstfps_old_pattern))
                            if not script_args.dryrun:
                                for dstfp_old in dstfps_old:
                                    os.remove(dstfp_old)

                elif assume_complete:
                    print("{} exist; skipping tile: {}".format(
                        "browse and meta" if tif_output_is_browse else "dem, browse, and meta",
                        matfile
                    ))
                    run_tile = False

            if run_tile:
                run_tile_matlist.append(matfile)

        run_tile_list = []

        if script_args.process_by == 'tile-file':
            run_tile_list = run_tile_matlist
        elif script_args.process_by == 'supertile-dir' and len(run_tile_matlist) > 0:
            supertile_dir = os.path.join(root_tiledir, supertile)
            run_tile_list = [supertile_dir]

        run_tile_list_all.extend(run_tile_list)
        # Un-comment the following line and un-indent the proceeding
        # `for` loop if you want to check & gather the full list of
        # processing tasks before looping through list for job submission.
    # for tile_path in run_tile_list_all:
        for tile_path in run_tile_list:
            job_id = None

            matlab_addpath = wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
            """)

            if not tile_path.endswith('.mat'):
                supertile_dir = tile_path
                job_id = os.path.basename(supertile_dir)

                task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                    {matlab_addpath}
                    batch_tiles2tif_v4('{supertile_dir}', '{tile_projstr}', 'resolution','{res_name}', 'outRasterType','{script_args.tif_output}' {batch_t2t_metaonly_args});
                """))

            else:
                tile_matfile = tile_path
                job_id, _ = os.path.splitext(os.path.basename(tile_matfile))

                if script_args.meta_only:
                    task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                        {matlab_addpath}
                        tileMetav4('{tile_matfile}');
                    """))
                else:
                    task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                        {matlab_addpath}
                        writeTileToTifv4('{tile_matfile}', '{tile_projstr}', 'outRasterType','{script_args.tif_output}');
                        tileMetav4('{tile_matfile}');
                    """))

            submit_cmd = job_handler.add_task_cmd(task_cmd, job_id)
            if submit_cmd is not None:
                print("Submitting job [{}/{}(est)]:\n  {}".format(
                    job_handler.job_num, job_handler.num_jobs, submit_cmd
                ))
                if not script_args.dryrun:
                    subprocess.call(submit_cmd, shell=True)

    if job_handler.bundle_tasks:
        submit_cmd = job_handler.get_last_bundle_submit_cmd()
        if submit_cmd is not None:
            print("Submitting job [{}/{}]:\n  {}".format(
                job_handler.job_num, job_handler.job_num, submit_cmd
            ))
            if not script_args.dryrun:
                subprocess.call(submit_cmd, shell=True)

    sys.exit(0)


if __name__ == '__main__':
    main()
