from __future__ import annotations

from dataclasses import dataclass, asdict
import datetime as dt
import json
from pathlib import Path
import subprocess

import click
from rema_mosaic.commands.utils import make_dirs_if_not_exist, Command
from rema_mosaic.config import Settings

BST_PARAMS_FILENAME = "bst_params.json"


@dataclass
class BSTParams:
    datefilt_start: str
    datefilt_end: str
    tempdir: str
    logdir: str

    def write_json(self, dst: Path) -> None:
        with open(dst, "w") as f:
            content = json.dumps(asdict(self), indent=2)
            f.write(content)

    @classmethod
    def read_json(cls, src: Path) -> BSTParams:
        with open(src, "r") as f:
            return cls(**json.loads(f.read()))


def _format_date(yyyymmdd: str) -> str:
    # Parse the input string to a datetime object
    date_object = dt.datetime.strptime(yyyymmdd, "%Y%m%d")
    # Format the datetime object to the desired output format
    formatted_date = date_object.strftime("%B %d, %Y")
    return formatted_date


@click.command(
    no_args_is_help=True, short_help="Create and configure working directory"
)
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Target directory (where tile subfolders will be created)",
)
@click.option(
    "--datefilt-start",
    required=True,
    type=click.STRING,
    help="Filter strip database to (first image) acquisition date on or after this date in 'YYYYMMDD' format",
)
@click.option(
    "--datefilt-end",
    required=True,
    type=click.STRING,
    help="Filter strip database to (first image) acquisition date on or before this date in 'YYYYMMDD' format",
)
@click.option(
    "--dryrun", is_flag=True, default=False, help="Print actions without executing"
)
def init_bst(
    dstdir: Path,
    datefilt_start: str,
    datefilt_end: str,
    dryrun: bool,
) -> None:
    """Create working directory and write BST configuration variables to JSON."""
    params_file = dstdir / BST_PARAMS_FILENAME
    tempdir = Path(dstdir) / "jobfiles"
    logdir = Path(dstdir) / "logs"

    must_exist = [dstdir, tempdir, logdir]

    # Print some of the inputs in a human-readable manner to help the user verify that
    # they entered what they intended
    click.echo(f"Directory Name: {dstdir.name}")
    click.echo(f"Start Date:     {_format_date(datefilt_start)}")
    click.echo(f"End Date:       {_format_date(datefilt_end)}")

    params = BSTParams(
        datefilt_start=datefilt_start,
        datefilt_end=datefilt_end,
        tempdir=f"{tempdir}",
        logdir=f"{logdir}",
    )

    click.echo(f"Saving config to: {params_file}")
    if not dryrun:
        if params_file.exists():
            click.echo(f"ERROR: {params_file} already exists", err=True)
            exit(1)

        make_dirs_if_not_exist(must_exist)
        params.write_json(params_file)


@click.command(no_args_is_help=True, short_help="Run the BST and MST stages")
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, resolve_path=True, path_type=Path
    ),
    help="Target directory (where tile subfolders will be created)",
)
@click.option(
    "--tiles",
    required=True,
    type=click.STRING,
    help="List of mosaic tiles; either specified on command line (comma delimited), or a text file list (each tile on separate line)",
)
@click.option(
    "--chain-mst/--no-chain-mst",
    default=True,
    show_default=True,
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
def run_bst(
    settings: Settings,
    dstdir: Path,
    tiles: str,
    chain_mst: bool,
    dryrun: bool,
    show_command: bool,
) -> None:
    """Run Build Subtiles (BST) followed by Mosaic Subtiles (MST)"""
    script = settings.SETSM_POSTPROCESSING_PYTHON_DIR / "batch_buildSubTiles.py"
    libdir = settings.SETSM_POSTPROCESSING_MATLAB_DIR
    jobscript = settings.SETSM_POSTPROCESSING_PYTHON_DIR / "qsub_buildSubTiles.sh"

    params_file = dstdir / BST_PARAMS_FILENAME
    if not params_file.exists():
        click.echo(f"ERROR: {params_file} does not exist", err=True)
        exit(1)
    params = BSTParams.read_json(params_file)

    # If any of the path options need to be changed in the future, migrate them to
    # the Settings object and configure via the .env file
    options = {
        "--project": "rema",
        "--epsg": "3031",
        "--tile-def": "/mnt/pgc/data/projects/earthdem/tiledef_files/rema_tile_definitions_plus_sgssi2.mat",
        "--strip-db": "/mnt/pgc/data/projects/earthdem/strip_databases/REMAdatabase4_2m_v4.1_20230120.mat",
        "--strips-dir": "/mnt/pgc/data/elev/dem/setsm/REMA/region",
        "--ref-dem": "/mnt/pgc/data/elev/dem/tandem-x/90m/mosaic/TanDEM-X_Antarctica_90m/old_mosaic/TanDEM_Antarctica_Mosaic.tif",
        # "--water-tile-dir": None,
        "--tileqc-dir": "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_qc/rema_v2.1_rema_strip_automosaic_qc_tile_v4_1",
        "--tileparam-list": "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_params/tileParamList_v13e.txt",
        "--datefilt-start": params.datefilt_start,
        "--datefilt-end": params.datefilt_end,
        "--libdir": f"{libdir}",
        "--jobscript": f"{jobscript}",
        "--tempdir": params.tempdir,
        "--logdir": params.logdir,
        "--queue": "batch",
        # "--swift-program": "/projects/sciteam/bazu/tools/swift-2/bin/swift",
        # "--swift-logdir": None,
        "--tasks-per-job": "1",
    }

    flags = [
        "--make-10m-only",
        "--rerun",
        # "--rerun-without-cleanup",
        # "--chain-mst",
        # "--chain-mst-keep-subtiles",
        # "--chain-mst-no-local",
        # "--bypass-finfile-req",
        # "--skip-missing-watermasks",
        # "--pbs",
        "--slurm",
        # "--swift",
        # "--dryrun",
    ]

    if chain_mst:
        flags.append("--chain-mst")
    if dryrun:
        flags.append("--dryrun")

    cmd = Command(
        program="python",
        script=f"{script}",
        args=[f"{dstdir}", tiles],
        options=options,
        flags=flags,
    )

    if show_command:
        click.echo(cmd.to_pretty_str())
        exit(0)

    subprocess.run(cmd.to_list())
