import shlex
import subprocess
from pathlib import Path

import click

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


def _base_cmd(
    settings: Settings, tiles: Path, utm_zone: UtmZone, reg_method: str
) -> str:
    tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "00-matfiles"
    ref_dem_path = (
        settings.REFERENCE_DEM_DIR / str(utm_zone) / "<supertile>_10m_cop30_wgs84.tif"
    )
    cover_tif_path = (
        settings.LANDCOVER_DIR
        / str(utm_zone)
        / "<supertile>_10m_esa_worldcover_2021.tif"
    )

    return "\n\t".join(
        [
            f"python {settings.SETSM_POSTPROCESSING_PYTHON_DIR}/batch_registerTiles.py",
            f"{tiledir} {tiles} earthdem 2",
            f"--reg-method {reg_method}",
            f"--ref-dem-path '{ref_dem_path}'",
            f"--cover-tif-path '{cover_tif_path}'",
            "--tile-org pgc",
            "--process-by supertile-dir",
            "--process-group separate",
            f"--matlib {settings.SETSM_POSTPROCESSING_MATLAB_DIR}",
        ],
    )


def _skipreg_option(skipreg_shp: Path) -> str:
    return f"--cop30-skipreg-shp {skipreg_shp}"


def _dryrun_option() -> str:
    return "--dryrun"


def _slurm_options(settings: Settings, job_name_prefix: str, utm_zone: UtmZone) -> str:
    job_config_file = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_config.ini"
    )
    job_script = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_<scheduler>.sh"
    )
    job_log_dir = settings.WORKING_ZONES_DIR / f"{utm_zone}" / "processing_logs"

    return "\n\t".join(
        [
            "--slurm",
            f"--job-config-file {job_config_file}",
            "--job-config-on",
            f"--job-script '{job_script}'",
            f"--job-name-prefix '{job_name_prefix}'",
            f"--job-log-dir {job_log_dir}",
            "--job-ncpus 4",
            "--job-mem 200",
            "--job-walltime 48",
            "--job-queue big_mem",
            "--tasks-per-job 1",
            "--tasks-per-job-mode serial",
            "--task-bundle-dir $HOME/scratch/task_bundles",
        ],
    )


@click.command(short_help="Apply coregistration to matfiles")
@click.option(
    "--skipreg-shp",
    type=EXISTING_FILE,
    help="Path to shapefile of polygon features with 'skipreg' field value set to 1 where, in the COP30 registration process, DEM data blobs should not be vertically adjusted.",
)
@click.option("--slurm", is_flag=True, help="Submit jobs to slurm cluster")
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.option(
    "--show-command",
    is_flag=True,
    help="Print the generated command without executing",
)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.argument(
    "tiles",
    nargs=1,
    type=EXISTING_FILE,
)
@click.pass_obj
def coreg_matfiles(
    settings: Settings,
    skipreg_shp: Path,
    slurm: bool,
    dryrun: bool,
    show_command: bool,
    utm_zone: UtmZone,
    tiles: Path,
) -> None:
    """Apply coregistration to the matfiles directly to produce '_reg.mat' versions of
    the tiles.
    """
    reg_method = "cop30"
    job_name_prefix = "reg<resolution>m"
    cmd = _base_cmd(
        settings=settings, tiles=tiles, utm_zone=utm_zone, reg_method=reg_method
    )
    if skipreg_shp:
        cmd = "\n\t".join([cmd, _skipreg_option(skipreg_shp)])
    if slurm:
        slurm_options = _slurm_options(
            settings=settings, job_name_prefix=job_name_prefix, utm_zone=utm_zone
        )
        cmd = "\n\t".join([cmd, slurm_options])
    if dryrun:
        cmd = "\n\t".join([cmd, _dryrun_option()])
    if show_command:
        print(cmd)
        return

    subprocess.run(shlex.split(cmd), check=True)


@click.command(short_help="Apply water flattening to matfiles")
@click.option("--slurm", is_flag=True, help="Submit jobs to slurm cluster")
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.option(
    "--show-command",
    is_flag=True,
    help="Print the generated command without executing",
)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.argument(
    "tiles",
    nargs=1,
    type=EXISTING_FILE,
)
@click.pass_obj
def water_flatten_matfiles(
    settings: Settings,
    slurm: bool,
    dryrun: bool,
    show_command: bool,
    utm_zone: UtmZone,
    tiles: Path,
) -> None:
    """Apply water flattening to the matfiles directly to produce '_reg_fill.mat'
    versions of the tiles.
    """
    reg_method = "fillWater"
    job_name_prefix = "reg<resolution>m"
    cmd = _base_cmd(
        settings=settings, tiles=tiles, utm_zone=utm_zone, reg_method=reg_method
    )
    if slurm:
        slurm_options = _slurm_options(
            settings=settings, job_name_prefix=job_name_prefix, utm_zone=utm_zone
        )
        cmd = "\n\t".join([cmd, slurm_options])
    if dryrun:
        cmd = "\n\t".join([cmd, _dryrun_option()])
    if show_command:
        print(cmd)
        return

    subprocess.run(shlex.split(cmd), check=True)
