from pathlib import Path

import click

EXISTING_FILE = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
    path_type=Path,
)

EXISTING_DIR = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    resolve_path=True,
    path_type=Path,
)
