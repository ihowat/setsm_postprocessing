from pathlib import Path
import subprocess

import click

from rema_mosaic.commands.utils import make_dirs_if_not_exist, Command
from rema_mosaic.config import Settings


@click.command(
    no_args_is_help=True, short_help="Create GeoTIFFs and metadata from matfiles"
)
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Target directory (where tile subfolders will be created)",
)
@click.option(
    "--tiles",
    required=True,
    type=click.STRING,
    help="List of mosaic tiles; either specified on command line (comma delimited), or a text file list (each tile on separate line)",
)
@click.option(
    "--dryrun", is_flag=True, default=False, help="Print actions without executing"
)
@click.option(
    "--show-command",
    is_flag=True,
    default=False,
    help="Print the generated command and exit",
)
@click.pass_obj
def export_products(
    settings: Settings, dstdir: Path, tiles: str, dryrun: bool, show_command: bool
) -> None:
    """Create final GeoTIFF and metadata products from the tile matfiles"""
    script = settings.SETSM_POSTPROCESSING_PYTHON_DIR / "batch_tiles2tif_v4.py"
    matlib = settings.SETSM_POSTPROCESSING_MATLAB_DIR
    slurm_log_dir = dstdir / "logs" / "slurm" / "export_products"
    task_bundle_dir = Path.home() / "scratch" / "task_bundles"
    job_name_prefix = f"t2t_{dstdir.name}"

    must_exist = [slurm_log_dir]
    if not dryrun:
        make_dirs_if_not_exist(must_exist)

    args = [
        f"{dstdir}",
        tiles,
        "rema",
        "10",
    ]

    options = {
        "--tif-format": "cog",
        "--tile-buffer-meters": "100",
        "--ref-dem-path": "/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles/<supertile>_10m_cop30_wgs84.tif",
        # "--cover-tif-path": None,
        # "--fill-water-interp-method": "2",
        # "--register-cop30-skipreg-shp": None,
        "--add-sea-surface-height": "true",
        "--final-qc-mask": "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/final_qc_mask/rema_v2.1/rema_final_mask_v2_1.mat",
        "--use-final-qc-mask": "true",
        "--keep-old-results": "output-set",
        "--tile-org": "pgc",
        "--process-by": "tile-file",
        "--matlib": f"{matlib}",
        "--job-config-file": "/project/vida/scratch/devin/projects/earthdem-v1-1/setsm_postprocessing_pgc/lib/jobscript_config.ini",
        "--job-script": "/project/vida/scratch/devin/projects/earthdem-v1-1/setsm_postprocessing_pgc/lib/jobscript_<scheduler>.sh",
        "--job-name-prefix": job_name_prefix,
        "--job-log-dir": f"{slurm_log_dir}",
        "--job-ncpus": "2",
        "--job-mem": "16",
        "--job-walltime": "4",
        "--job-queue": "batch",
        "--tasks-per-job": "1",
        "--tasks-per-job-mode": "serial",
        "--task-bundle-dir": f"{task_bundle_dir}",
    }

    flags = [
        # "--tile-nocrop",
        # "--apply-ref-filter",
        # "--apply-topo-filter",
        # "--apply-water-fill",
        # "--register-to-ref",
        # "--register-to-ref-debug",
        # "--rerun",
        # "--dryrun",
        "--slurm",
    ]

    if dryrun:
        flags.append("--dryrun")

    cmd = Command(
        program="python", script=f"{script}", args=args, options=options, flags=flags
    )
    if show_command:
        print(cmd.to_pretty_str())
        return

    subprocess.run(cmd.to_list())


@click.command(no_args_is_help=True, short_help="Create VRTs of GeoTIFF tiles")
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Target directory (where tile subfolders will be created)",
)
@click.option(
    "--dryrun", is_flag=True, default=False, help="Print actions without executing"
)
def build_vrts(dstdir: Path, dryrun: bool) -> None:
    """Build VRTs from all GeoTIFF products in the tile directories."""
    suffixes = [
        "_browse.tif",
        "_countmt.tif",
        "_count.tif",
        "_datamask.tif",
        "_dem.tif",
        "_mad.tif",
        "_maxdate.tif",
        "_mindate.tif",
    ]

    for suffix in suffixes:
        files = [str(f) for f in dstdir.glob(f"**/*{suffix}")]
        vrt = dstdir / f"{dstdir.name}{suffix.replace('.tif', '.vrt')}"

        click.echo(f"Found {len(files)} files with suffix '{suffix}'")
        if not dryrun and len(files) > 0:
            cmd = Command(program="gdalbuildvrt", args=[f"{vrt}", *files])
            click.echo(f"Creating {vrt.name}")
            subprocess.run(cmd.to_list())
