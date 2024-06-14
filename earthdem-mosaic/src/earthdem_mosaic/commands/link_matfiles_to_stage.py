import os
from pathlib import Path

import click

from earthdem_mosaic.commands._utils import EXISTING_DIR


@click.command(
    short_help="Hardlink files from one stage to another",
)
@click.option(
    "--src",
    help="Stage directory to search for matfiles to link.",
    type=EXISTING_DIR,
    required=True,
)
@click.option(
    "--dst",
    help="Stage directory to copy matfiles to.",
    type=EXISTING_DIR,
    required=True,
)
@click.option(
    "--src-suffix",
    help="Filename pattern used to find source matfiles. (e.g. '_reg_fill.mat')",
    type=str,
    required=True,
)
@click.option(
    "--dst-suffix",
    help="Replace the provided --src-suffix with this value when naming the linked file.",
    type=str,
    required=False,
)
@click.option("-v", "--verbose", is_flag=True)
def link_files_to_stage(
    src: Path,
    dst: Path,
    src_suffix: str,
    dst_suffix: str,
    verbose: bool,
) -> None:
    """Create hardlinks for files in one stage to the coresponding location in another
    stage. Optionally replace the suffix of the destination file during linking.
    """
    click.echo(
        f"Searching {src} for files with suffix: {src_suffix}",
    )
    matfiles = src.glob(pattern=f"**/*{src_suffix}")

    for matfile in matfiles:
        dst_str_path = str(matfile).replace(str(src), str(dst))
        if dst_suffix:
            dst_str_path = dst_str_path.replace(src_suffix, dst_suffix)

        destination = Path(dst_str_path)
        destination.parent.mkdir(parents=True, exist_ok=True)

        click.echo(f"cp --link {matfile} {destination}") if verbose else None
        os.link(src=matfile, dst=destination)

    click.echo(f"Matched files linked to: {dst}")
    if dst_suffix:
        click.echo(f"Suffix of linked files changed from {src_suffix} to {dst_suffix}")
