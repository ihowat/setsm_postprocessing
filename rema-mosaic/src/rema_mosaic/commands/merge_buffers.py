import subprocess
import click
import shlex
import shutil
from pathlib import Path

from rema_mosaic.commands.utils import (
    parse_tiles_input,
    Command,
    make_dirs_if_not_exist,
)
from rema_mosaic.config import Settings


@click.command(no_args_is_help=True, short_help="Make backup copies of matfiles")
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, resolve_path=True, path_type=Path
    ),
    help="Target directory (where tile subfolders have been created)",
)
@click.option(
    "--tiles",
    required=True,
    type=click.STRING,
    help="List of mosaic tiles; either specified on command line (comma delimited), or a text file list (each tile on separate line)",
)
@click.option(
    "--tag",
    required=True,
    type=click.STRING,
    help="Text to append to the end of the copied filename. (e.g. '--tag bstmst' will create <filename>.mat.bstmst)",
)
@click.option(
    "--dryrun", is_flag=True, default=False, help="Print actions without executing"
)
def backup_matfiles(dstdir: Path, tiles: str, tag: str, dryrun: bool) -> None:
    tile_list = parse_tiles_input(tiles)
    for tile in tile_list:
        for src in (dstdir / tile).glob("*_10m.mat"):
            dst = src.with_suffix(f".mat.{tag}")

            if dst.exists():
                click.echo(f"Skip: {src.name} -> {dst.name} already exists")
                continue

            click.echo(f"Copy: {src.name} -> {dst.name}")
            if not dryrun:
                shutil.copy2(src, dst)


@click.command(no_args_is_help=True, short_help="Build index for merge-tile-buffers")
@click.option(
    "--dstdir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, resolve_path=True, path_type=Path
    ),
    help="Target directory (where tile subfolders have been created)",
)
@click.option(
    "--show-command",
    is_flag=True,
    default=False,
    help="Print the generated command and exit",
)
@click.pass_obj
def create_neighbor_index(settings: Settings, dstdir: Path, show_command: bool) -> None:
    """Create the tile neighbor index matfile used to merge tile buffers"""
    libdir = settings.SETSM_POSTPROCESSING_MATLAB_DIR
    outfile = dstdir / "tile_index_files" / "tileNeighborIndex_10m.mat"
    priority_suffix = "_10m.mat"
    secondary_suffix = "_match_nothing.mat"

    must_exist = [outfile.parent]

    cmd = Command(
        program="matlab",
        flags=[
            "-nojvm",
            "-nodisplay",
            "-nosplash",
        ],
        options={
            "-r": f""" "try; addpath('{libdir}'); tileNeighborIndex( '{dstdir}', 'org','pgc', 'resolution','10m', 'outfile','{outfile}', 'priority_suffix','{priority_suffix}', 'secondary_suffix','{secondary_suffix}' );; catch e; disp(getReport(e)); exit(1); end; exit(0);" """,
        },
    )

    if show_command:
        click.echo(cmd.to_pretty_str())
        exit(0)

    if not show_command:
        make_dirs_if_not_exist(must_exist)
        click.echo(cmd.to_pretty_str())
        # This pattern of creating a string just to re-split it is to resolve quoting
        # peculiarities of the matlab script string.
        cmd_str = " ".join(cmd.to_list())
        subprocess.run(shlex.split(cmd_str), check=True)


@click.command(no_args_is_help=True, short_help="Merge overlapping tile areas")
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
    "--axis",
    type=click.Choice(["row", "column"]),
    nargs=1,
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
def merge_tile_buffers(
    settings: Settings,
    dstdir: Path,
    tiles: str,
    axis: str,
    dryrun: bool,
    show_command: bool,
) -> None:
    """Adjust overlapping tile buffers to merge adjacent tiles along a given axis"""
    python_lib = settings.SETSM_POSTPROCESSING_PYTHON_DIR
    script = python_lib / "batch_mergeTileBuffer.py"
    job_config_file = python_lib / "lib" / "jobscript_config.ini"
    job_script = python_lib / "lib" / "jobscript_<scheduler>.sh"
    matlib = settings.SETSM_POSTPROCESSING_MATLAB_DIR
    resolution = "10"
    tile_index = dstdir / "tile_index_files" / "tileNeighborIndex_10m.mat"
    job_log_dir = dstdir / "logs" / "slurm" / "merge" / "10m"
    task_bundle_dir = Path.home() / "scratch" / "task_bundles"

    must_exist = [job_log_dir]
    if not dryrun:
        make_dirs_if_not_exist(must_exist)

    args = [f"{dstdir}", tiles, resolution]

    options = {
        "--tile-index": f"{tile_index}",
        "--tile-org": "pgc",
        "--process-group": axis,
        "--matlib": f"{matlib}",
        "--job-config-file": f"{job_config_file}",
        "--job-script": f"{job_script}",
        "--job-name-prefix": f"mtb_{axis}",
        "--job-log-dir": f"{job_log_dir}",
        "--job-ncpus": "4",
        "--job-mem": "16",
        "--job-walltime": "48",
        "--job-queue": "batch",
        "--tasks-per-job": "1",
        "--tasks-per-job-mode": "serial",
        "--task-bundle-dir": f"{task_bundle_dir}",
    }

    flags = [
        "--slurm",
    ]

    if dryrun:
        flags.append("--dryrun")

    cmd = Command(
        program="python",
        script=f"{script}",
        args=args,
        flags=flags,
        options=options,
    )

    if show_command:
        click.echo(cmd.to_pretty_str())
        exit(0)

    subprocess.run(cmd.to_list())
