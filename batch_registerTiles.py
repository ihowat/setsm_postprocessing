
import argparse
import glob
import os
import subprocess
import sys
import time

from lib import jobscript_utils
from lib import tilelist_utils
from lib.jobscript_utils import wrap_multiline_str


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
            Generate "registered" '<tilename>_<resolution>_reg.mat'
            DEM mosaic tile files next to original "unregistered" files,
            where the registered files have their raster data arrays
            corrected by adding a quadratic surface fit to input ICESat-2
            altimetry data.
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
        '--is2dir',
        type=str,
        default="/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile/",
        help="Path ICESat-2 altimetry-by-supertile files."
    )

    parser.add_argument(
        '--skip-dzfit',
        action='store_true',
        help=wrap_multiline_str("""
            Skip the dzfit step of ICESat-2 registration application.
            You should skip applying dzfit to tiles over the ice shelves of Antarctica.
        """)
    )
    parser.add_argument(
        '--dzfit-minpoints',
        type=int,
        default=1000,
        help=wrap_multiline_str("""
            Minimum number of valid ICESat-2 points overalapping tile
            for dzfit to be applied.
        """)
    )

    parser.add_argument(
        '--tile-org',
        type=str,
        # choices=['pgc', 'osu'],
        choices=['pgc'],
        default='pgc',
        help="Tile file organization scheme."
    )
    parser.add_argument(
        '--process-by',
        type=str,
        # choices=['supertile-dir', 'tile-file'],
        choices=['supertile-dir'],
        default='supertile-dir',
        help=wrap_multiline_str(r"""
            \nIf 'supertile-dir', if one tile within the supertile folder
            needs to be processed, hand whole supertile folder to Matlab script
            to be processed as a single task.
            \nIf 'tile-file', send each tile to be processed individually
            to Matlab script as a single task.
        """)
    )
    parser.add_argument(
        '--process-group',
        type=str,
        # choices=['separate', 'all'],
        choices=['separate'],
        default='separate',
        help="Per-job grouping of tiles for processing."
    )
    
    parser.add_argument(
        '--matlib',
        type=str,
        default=MATLAB_LIBDIR,
        help="Path to directory containing necessary Matlab functions."
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

    res_name = '{}m'.format(script_args.resolution)

    jobscript_utils.adjust_args(script_args, arg_parser)
    jobscript_utils.create_dirs(script_args, arg_parser)


    ## Verify arguments

    if not os.path.isdir(root_tiledir):
        arg_parser.error("Argument 'tiledir' is not an existing directory: {}".format(root_tiledir))

    if not os.path.isdir(script_args.is2dir):
        arg_parser.error("--is2dir is not an existing directory: {}".format(script_args.is2dir))

    jobscript_utils.verify_args(script_args, arg_parser)


    ## Test tiles exist and group into mosaic groups

    run_supertiles = set()
    for supertile in supertile_list:
        run_tile = True

        supertile_folder = root_tiledir if script_args.tile_org == 'osu' else os.path.join(root_tiledir, supertile)
        tile_basename_pat = '{}{}{}'.format(supertile, '_*_' if script_args.resolution == 2 else '_', res_name)

        matfile_unreg_pattern = os.path.join(supertile_folder, '{}.mat'.format(tile_basename_pat))
        matfile_unreg_list = glob.glob(matfile_unreg_pattern)
        if len(matfile_unreg_list) == 0:
            print("Tile {} {} unreg .mat file does not exist: {}".format(supertile, res_name, matfile_unreg_pattern))
            run_tile = False
        else:
            matfile_reg_list = [f.replace('.mat', '_reg.mat') for f in matfile_unreg_list]
            matfile_reg_exist_list = [f for f in matfile_reg_list if os.path.isfile(f)]
            if len(matfile_reg_exist_list) == len(matfile_unreg_list):
                run_tile = False

        finfile_pattern = os.path.join(supertile_folder, '{}.fin'.format(tile_basename_pat))
        finfile_list = glob.glob(finfile_pattern)
        num_finfiles_exist = len(finfile_list)
        num_finfiles_expected = 4 if script_args.resolution == 2 else 1
        if num_finfiles_exist != num_finfiles_expected:
            print("WARNING: {}/{} of expected {} MST finfiles exist for tile {}".format(
                num_finfiles_exist, num_finfiles_expected, res_name, finfile_pattern))
            
        if not run_tile:
            continue
            
        altimetry_file = os.path.join(script_args.is2dir, '{}_is2.mat'.format(supertile))
        if not os.path.isfile(altimetry_file):
            print("Tile {} altimetry file does not exist: {}".format(supertile, altimetry_file))
            run_tile = False
            
        if not run_tile:
            continue

        tif_meta_patterns = [
            os.path.join(supertile_folder, '{}_*.tif'.format(tile_basename_pat)),
            os.path.join(supertile_folder, '{}_meta.txt'.format(tile_basename_pat))
        ]
        tif_meta_files = [f for pat in tif_meta_patterns for f in glob.glob(pat)]
        if len(tif_meta_files) > 0:
            if len(matfile_unreg_list) == 0:
                print("ERROR! Tile .mat file does not exist, but tif/meta results exist matching {}".format(tif_meta_patterns))
                continue
            print("{}Removing old tif/meta results matching {}".format('(dryrun) ' if script_args.dryrun else '', tif_meta_patterns))
            if not script_args.dryrun:
                for f in tif_meta_files:
                    os.remove(f)
                    
        if run_tile:
            run_supertiles.add(supertile)

    run_supertiles = sorted(list(run_supertiles))


    ## Group tiles for processing

    group_tiles_dict, group_args_dict = tilelist_utils.group_supertiles_for_processing(script_args, run_supertiles)
    group_key_list = sorted(list(group_tiles_dict.keys()))
    num_groups = len(group_key_list)
    num_tiles = sum([len(tile_list) for tile_list in list(group_tiles_dict.values())])

    if num_groups == 0:
        print("No tiles to process")
        sys.exit(0)
    else:
        print("Processing {} tiles in {} groups".format(num_tiles, num_groups))

    wait_seconds = 4 if script_args.dryrun else 6
    if wait_seconds > 0:
        print("Pausing for {} seconds before job submission ({})".format(
            wait_seconds,
            'dryrun' if script_args.dryrun else 'REAL run'
        ))
        time.sleep(wait_seconds)


    ## Main processing loop

    num_tasks = num_groups
    job_handler = jobscript_utils.JobHandler(
        script_args, num_tasks,
        init_env_requests='matlab'
    )

    for group_key in group_key_list:
        job_id = group_key

        supertile = None
        tile_list = sorted(list(group_tiles_dict[group_key]))
        if len(tile_list) == 1 and script_args.process_by == 'supertile-dir' and script_args.process_group == 'separate':
            supertile = tile_list[0]
        else:
            arg_parser.error("Cannot process more than one supertile per processing group (see --process-by and --process-group args)")

        supertile_dir = os.path.join(root_tiledir, supertile)

        matscript_args = ''
        if group_args_dict is not None:
            matscript_args += ", {}".format(group_args_dict[group_key])
        if script_args.process_group in ('row', 'column'):
            if script_args.process_group == 'row':
                matscript_args += ", 'rows'"
            elif script_args.process_group == 'column':
                matscript_args += ", 'cols'"
        if script_args.skip_dzfit:
            matscript_args += ", 'skipDzfit'"
        matscript_args += ", 'dzfitMinPoints',{}".format(script_args.dzfit_minpoints)

        task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
            addpath('{SCRIPT_DIR}');
            addpath('{script_args.matlib}');
            batchRegisterTiles('{supertile_dir}', '{script_args.is2dir}', 'resolution','{res_name}' {matscript_args});
        """))

        submit_cmd = job_handler.add_task_cmd(task_cmd, job_id)
        if submit_cmd is not None:
            print("Submitting job ({}/{}):\n  {}".format(
                job_handler.job_num, job_handler.num_jobs, submit_cmd
            ))
            if not script_args.dryrun:
                subprocess.call(submit_cmd, shell=True)

    if job_handler.bundle_tasks:
        submit_cmd = job_handler.get_last_bundle_submit_cmd()
        if submit_cmd is not None:
            print("Submitting job ({}/{}):\n  {}".format(
                job_handler.job_num, job_handler.num_jobs, submit_cmd
            ))
            if not script_args.dryrun:
                subprocess.call(submit_cmd, shell=True)

    sys.exit(0)


if __name__ == '__main__':
    main()
