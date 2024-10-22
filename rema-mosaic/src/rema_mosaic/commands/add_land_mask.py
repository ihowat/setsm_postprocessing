from pathlib import Path
import subprocess

import click

from rema_mosaic.commands.utils import (
    parse_tiles_input,
    Command,
    make_dirs_if_not_exist,
)
from rema_mosaic.config import Settings


@click.command(no_args_is_help=True, short_help="Add land mask to matfiles")
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
@click.pass_obj
def add_land_mask(settings: Settings, dstdir: Path, tiles: str, dryrun: bool) -> None:
    """Adds the 'land' field to the tile matfiles.

    Note: The project-specific source tile definitions matfile is hardcoded in qsub_addLandMask2REMATiles.sh."""

    script = settings.SETSM_POSTPROCESSING_PYTHON_DIR / "qsub_addLandMask2REMATiles.sh"
    matlib = settings.SETSM_POSTPROCESSING_MATLAB_DIR
    slurm_log = dstdir / "logs" / "slurm" / "land_mask" / r"%x.o%j"

    must_exist = [slurm_log.parent]
    if not dryrun:
        make_dirs_if_not_exist(must_exist)

    for tile in parse_tiles_input(tiles):
        options = {
            "--job-name": f"add_land_mask_{dstdir.name}_{tile}",
            "--export": f"ARG_TILEDIR={dstdir / tile},ARG_MATLIB={matlib}",
            "--output": f"{slurm_log}",
            "--error": f"{slurm_log}",
        }
        cmd = Command(
            program="sbatch",
            options=options,
            # Use postfix to force the script to be at the end of the generated command
            postfix=f"{script}",
        )
        click.echo(cmd.to_pretty_str())
        if not dryrun:
            subprocess.run(cmd.to_list(), check=True)
