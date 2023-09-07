
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


## Argument defaults by 'project'

domain_choices = [
    'arcticdem',
    'rema',
    'earthdem',
]

epsg_projstr_dict = {
    None: '',
    3413: 'polar stereo north',
    3031: 'polar stereo south',
}
domain_epsg_dict = {
    'arcticdem': 3413,
    'earthdem':  None,
    'rema':      3031,
}

domain_regMethod_dict = {
    'arcticdem': 'cop30',
    'earthdem':  'cop30',
    'rema':      'is2'
}
supertile_key = '<supertile>'
domain_is2path_dict = {
    'arcticdem': None,
    'earthdem': None,
    'rema': "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile/<supertile>_is2.mat",
}
domain_refDemPath_dict = {
    'arcticdem': "/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles_wgs84/<supertile>_10m_cop30_wgs84.tif",
    'earthdem': None,
    'rema': "/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles_wgs84/<supertile>_10m_cop30_wgs84.tif",
}
domain_coverTifPath_dict = {
    'arcticdem': "/mnt/pgc/data/thematic/landcover/esa_worldcover_2020/mosaics/arctic_tiles/<supertile>_10m_cover.tif",
    'earthdem': None,
    'rema': None,
}
domain_cop30SkipregShp_dict = {
    'arcticdem': "/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/arcticdem_v4.1_mosaic_reg-cop30_skip-reg.shp",
    'earthdem': None,
    'rema': None,
}


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
        '--epsg',
        type=int,
        choices=list(epsg_projstr_dict.keys()),
        default=None,
        help=wrap_multiline_str(r"""
            Mosaic tile projection EPSG code.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_epsg_dict.items()]
        )))
    )

    parser.add_argument(
        '--reg-method',
        type=str,
        choices=['is2', '2to10', 'cop30', 'fillWater', 'seaSurf'],
        default=None,
        help=wrap_multiline_str(r"""
            Registration method to employ.
            \n'is2': Register tile raster data to ICESat-2 altimetry data.
            \n'2to10': Register 2m 50x50km quad tile raster data to corresponding 100x100km 10m supertile raster data.
            \n'cop30': Vertically register contiguous chunks ("blobs") of tile raster data to --ref-dem-path DEM data
            (Copernicus GLO-30 30m dataset tiles cut to match mosaic tile schema).
            \n'fillWater': Run 'fillWater.m' to generate output '_fill.mat' files, where the 'z' DEM array has
            areas marked as water in --cover-tif-path (ESA WorldCover dataset tiles) filled with the corresponding
            surface from --ref-dem-path DEM data.
            \n'seaSurf': Run 'addSeaSurfaceHeight.m' to to generated output '_fill.mat' files, where the 'z' DEM array
            has pixels over ocean (indicated by the tile matfile internal 'land' array) filled with the EGM2008
            geoid height. That tile matfile 'land' array can be set through a script like 'addLandMask2REMATiles.m'
            or through the '2to10' 2m to 10m tile registration process, which copies+subsets an existing 'land' array
            from the 10m tile matfile to the output 2m '_reg.mat' file.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_regMethod_dict.items()]
        )))
    )

    parser.add_argument(
        '--ref-dem-path',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Path to reference DEM used for 'cop30' registration method, or 'fillWater'.
            The path should include the '{}' substring to be replaced with supertile name.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_refDemPath_dict.items()]
        )))
    )
    parser.add_argument(
        '--cover-tif-path',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Path to landcover raster used for 'cop30' registration method or 'fillWater'.
            The path should include the '{}' substring to be replaced with supertile name.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_coverTifPath_dict.items()]
        )))
    )

    parser.add_argument(
        '--is2-path',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Path to '*_is2.mat' ICESat-2 altimetry-by-supertile files.
            The path should include the '{}' substring to be replaced with supertile name.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_is2path_dict.items()]
        )))
    )

    parser.add_argument(
        '--is2-skip-dzfit',
        action='store_true',
        help=wrap_multiline_str(r"""
            (Only applicable when --reg-method='is2')
            \nSkip the dzfit step of ICESat-2 registration application.
            \nYou should skip applying dzfit to tiles over the ice shelves of Antarctica.
        """)
    )
    parser.add_argument(
        '--is2-dzfit-minpoints',
        type=int,
        default=1000,
        help=wrap_multiline_str(r"""
            (Only applicable when --reg-method='is2')
            \nMinimum number of valid ICESat-2 points overalapping tile
            for dzfit to be applied.
        """)
    )
    
    parser.add_argument(
        '--cop30-skipreg-shp',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            (Only applicable when --reg-method='cop30')
            \nPath to shapefile of polygon features with 'skipreg' field value set to 1
            where, in the COP30 registration process, DEM data blobs should not be
            vertically adjusted.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_cop30SkipregShp_dict.items()]
        )))
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

    if script_args.epsg is None:
        script_args.epsg = domain_epsg_dict[script_args.domain]
    projection_string = epsg_projstr_dict[script_args.epsg]
    if script_args.reg_method is None:
        script_args.reg_method = domain_regMethod_dict[script_args.domain]
    if script_args.is2_path is None:
        script_args.is2_path = domain_is2path_dict[script_args.domain]
    if script_args.ref_dem_path is None:
        script_args.ref_dem_path = domain_refDemPath_dict[script_args.domain]
    if script_args.cover_tif_path is None:
        script_args.cover_tif_path = domain_coverTifPath_dict[script_args.domain]

    jobscript_utils.adjust_args(script_args, arg_parser)
    jobscript_utils.create_dirs(script_args, arg_parser)

    script_args.job_name_prefix = f"{script_args.job_name_prefix}_{script_args.reg_method}"


    ## Verify arguments

    if not os.path.isdir(root_tiledir):
        arg_parser.error("Argument 'tiledir' is not an existing directory: {}".format(root_tiledir))

    if script_args.reg_method == 'is2':
        if not script_args.is2_path:
            arg_parser.error("--is2-path cannot be empty --reg-method='is2'")
        else:
            check_dir, _, _ = script_args.is2_path.partition(supertile_key)
            check_dir = os.path.dirname(check_dir)
            if not os.path.isdir(check_dir):
                arg_parser.error("--is2-path directory does not exist: {}".format(check_dir))

    if script_args.reg_method in ('cop30', 'fillWater'):
        if (not script_args.ref_dem_path) or (not script_args.cover_tif_path):
            arg_parser.error("--ref-dem-path and --cover-tif-path cannot be empty"
                             " when --reg-method is one of ['cop30', 'fillWater']")
        else:
            check_dir, _, _ = script_args.ref_dem_path.partition(supertile_key)
            check_dir = os.path.dirname(check_dir)
            if not os.path.isdir(check_dir):
                arg_parser.error("--ref-dem-path directory does not exist: {}".format(check_dir))
            check_dir, _, _ = script_args.cover_tif_path.partition(supertile_key)
            check_dir = os.path.dirname(check_dir)
            if not os.path.isdir(check_dir):
                arg_parser.error("--cover-tif-path directory does not exist: {}".format(check_dir))

    if script_args.reg_method == 'cop30':
        check_file = script_args.cop30_skipreg_shp
        if check_file and not os.path.isfile(check_file):
            arg_parser.error("--cop30-skipreg-shp file does not exist: {}".format(check_file))

    if script_args.reg_method == 'seaSurf':
        if not projection_string:
            arg_parser.error("Projection reference string is not defined for --epsg='{}'".format(script_args.epsg))

    jobscript_utils.verify_args(script_args, arg_parser)


    ## Test tiles exist and group into mosaic groups

    run_supertiles = set()
    for supertile in supertile_list:
        run_tile = True

        if supertile.startswith('utm'):
            if script_args.domain != 'earthdem':
                arg_parser.error("domain should be 'earthdem' when 'utm*' prefix tilenames are provided")
        else:
            if script_args.domain == 'earthdem':
                arg_parser.error("domain should NOT be 'earthdem' when tilenames do not have 'utm*' prefix")


        supertile_folder = root_tiledir if script_args.tile_org == 'osu' else os.path.join(root_tiledir, supertile)
        tile_basename_pat = '{}{}{}'.format(supertile, '_*_' if script_args.resolution == 2 else '_', res_name)

        if script_args.reg_method in ('fillWater', 'seaSurf'):
            matfile_unreg_pattern = os.path.join(supertile_folder, '{}.mat'.format(tile_basename_pat))
            matfile_unreg_list = glob.glob(matfile_unreg_pattern)

            matfile_reg_pattern = os.path.join(supertile_folder, '{}_reg.mat'.format(tile_basename_pat))
            matfile_reg_files = glob.glob(matfile_reg_pattern)
            matfile_unreg_files = glob.glob(matfile_unreg_pattern)
            matfile_unreg_files = list(set(matfile_unreg_files).difference(set([f.replace('_reg.mat', '.mat') for f in matfile_reg_files])))
            matfile_list = matfile_reg_files + matfile_unreg_files
            matfile_fill_list = [f.replace('.mat', '_fill.mat') for f in matfile_list]
            matfile_fill_exist_list = [f for f in matfile_fill_list if os.path.isfile(f)]
            if len(matfile_fill_exist_list) == len(matfile_list):
                run_tile = False

        else:
            matfile_unreg_pattern = os.path.join(supertile_folder, '{}.mat'.format(tile_basename_pat))
            matfile_unreg_list = glob.glob(matfile_unreg_pattern)
            if len(matfile_unreg_list) == 0:
                print("Tile {} {} unreg .mat file does not exist: {}".format(supertile, res_name, matfile_unreg_pattern))
                run_tile = False
            else:
                matfile_regcomplete_list = []
                for f_unreg in matfile_unreg_list:
                    f_reg = f_unreg.replace('.mat', '_reg.mat')
                    f_regfail = f_unreg.replace('.mat', '_reg.mat.regfail')
                    if os.path.isfile(f_reg):
                        matfile_regcomplete_list.append(f_reg)
                    elif os.path.isfile(f_regfail):
                        matfile_regcomplete_list.append(f_regfail)
                if len(matfile_unreg_list) == len(matfile_regcomplete_list):
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

        if script_args.reg_method == 'is2':
            altimetry_file = script_args.is2_path.replace(supertile_key, supertile)
            if not os.path.isfile(altimetry_file):
                print("ERROR: Tile {} altimetry file does not exist: {}".format(supertile, altimetry_file))
                print("Running this tile anyways in order to produce tile *_unreg.mat.bak backup tile matfile")
                # run_tile = False

        elif script_args.reg_method == '2to10':
            matfile_10m_unreg = os.path.join(supertile_folder, '{}_10m.mat'.format(supertile))
            matfile_10m_reg = os.path.join(supertile_folder, '{}_10m_reg.mat'.format(supertile))
            if not (os.path.isfile(matfile_10m_unreg) or os.path.isfile(matfile_10m_reg)):
                print("ERROR: Tile {} 10m DEM matfile ('_reg.mat' or unreg) does not exist: {}".format(supertile, matfile_10m_unreg))
                run_tile = False

        elif script_args.reg_method in ('cop30', 'fillWater'):
            auxfile_dne = False
            ref_dem_file = script_args.ref_dem_path.replace(supertile_key, supertile)
            cover_tif_file = script_args.cover_tif_path.replace(supertile_key, supertile)
            if not os.path.isfile(ref_dem_file):
                print("ERROR: Tile {} reference DEM file does not exist: {}".format(supertile, ref_dem_file))
                auxfile_dne = True
            if not os.path.isfile(cover_tif_file):
                print("ERROR: Tile {} water mask file does not exist: {}".format(supertile, cover_tif_file))
                auxfile_dne = True
            if auxfile_dne:
                if script_args.reg_method == 'cop30':
                    print("Running this tile anyways in order to produce tile *_unreg.mat.bak backup tile matfile")
                    # run_tile = False
                else:
                    run_tile = False
            
        if not run_tile:
            continue

        tif_meta_patterns = [
            os.path.join(supertile_folder, '{}_*.tif'.format(tile_basename_pat)),
            os.path.join(supertile_folder, '{}_meta.txt'.format(tile_basename_pat))
        ]
        tif_meta_files = [f for pat in tif_meta_patterns for f in glob.glob(pat)]
        tif_meta_files = [f for f in tif_meta_files if '_debug-reg_' not in f]
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
        init_env_requests='matlab gdal'
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

        if script_args.reg_method == 'is2':
            altimetry_file = script_args.is2_path.replace(supertile_key, supertile)

            if script_args.is2_skip_dzfit:
                matscript_args += ", 'skipDzfit'"
            matscript_args += ", 'dzfitMinPoints',{}".format(script_args.is2_dzfit_minpoints)

            task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
                batchRegisterTiles('{supertile_dir}', '{altimetry_file}', 'resolution','{res_name}' {matscript_args});
            """))

        elif script_args.reg_method == '2to10':

            task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
                batchRegister2mTileTo10mTile('{supertile_dir}' {matscript_args});
            """))

        elif script_args.reg_method in ('cop30', 'fillWater'):
            ref_dem_file = script_args.ref_dem_path.replace(supertile_key, supertile)
            cover_tif_file = script_args.cover_tif_path.replace(supertile_key, supertile)

            if script_args.reg_method == 'cop30':
                matlab_script = 'batchRegisterTilesToCOP30'
                matscript_args += ", 'registerBlobs'"
                if script_args.cop30_skipreg_shp:
                    matscript_args += ", 'registerBlobsSkipregShp',{}".format(script_args.cop30_skipreg_shp)

            elif script_args.reg_method == 'fillWater':
                matlab_script = 'batchWaterFillTiles'

            task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
                {matlab_script}('{supertile_dir}', '{ref_dem_file}', '{cover_tif_file}', 'resolution','{res_name}' {matscript_args});
            """))

        elif script_args.reg_method == 'seaSurf':

            task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
                batchAddSeaSurf('{supertile_dir}', '{projection_string}', 'resolution','{res_name}' {matscript_args});
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
