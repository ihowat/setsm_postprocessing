import os
import shutil
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
@click.option(
    "--link/--copy",
    default=True,
    help="Choice between hardlink or full copy",
    show_default=True,
)
@click.option("-v", "--verbose", is_flag=True)
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
def link_files_to_stage(
    src: Path,
    dst: Path,
    src_suffix: str,
    dst_suffix: str,
    link: bool,
    verbose: bool,
    dryrun: bool,
) -> None:
    """Create hardlinks (or full copies with --copy) for files in one stage to the
    coresponding location in another stage. Optionally replace the suffix of the
    destination file.
    """
    click.echo(
        f"Searching {src} for files with suffix: {src_suffix}",
    )
    src_files = src.glob(pattern=f"**/*{src_suffix}")

    for src_file in src_files:
        dst_str_path = str(src_file).replace(str(src), str(dst))
        if dst_suffix:
            dst_str_path = dst_str_path.replace(src_suffix, dst_suffix)

        dst_file = Path(dst_str_path)
        dst_file.parent.mkdir(parents=True, exist_ok=True)

        if dryrun or verbose:
            mock_cmd = [
                "cp",
                "--link" if link else "",
                f"{src_file}",
                f"{dst_file}",
            ]
            click.echo(" ".join(mock_cmd))
        if not dryrun:
            if link:
                os.link(src=src_file, dst=dst_file)
            else:
                shutil.copy2(src=src_file, dst=dst_file)

    click.echo(f"Matched files linked to: {dst}")
    if dst_suffix:
        click.echo(f"Suffix of linked files changed from {src_suffix} to {dst_suffix}")
