from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from matlib.funcs import utils

PROJECT_CHOICES = Literal["ArcticDEM", "EarthDEM", "REMA"]


@dataclass
class TileMetaV4Options:
    """
    Configuration options for TileMetaV4 operations.

    Attributes:
        project (PROJECT_CHOICES | None): The DEM project to use (ArcticDEM, EarthDEM, or REMA).
                                         Defaults to None.
        tileVersion (str | None): The version of the tile to use. Defaults to None.
        overwrite (bool): Whether to overwrite existing files. Defaults to False.
    """

    project: PROJECT_CHOICES | None = None
    tileVersion: str | None = None
    overwrite: bool = False

    @staticmethod
    def flag_keys() -> set[str]:
        return {"overwrite"}


@dataclass
class TileMetaV4:
    """
    Handles the creation of tileMetaV4 commands for digital elevation model operations.

    Attributes:
        tilef (str): The tile file path or identifier.
        options (TileMetaV4Options): Configuration options for the tile metadata operation.
    """

    tilef: str
    options: TileMetaV4Options

    def __post_init__(self):
        if not Path(self.tilef).is_absolute():
            raise ValueError(
                f"Expected absolute path for tilef, received relative path: {self.tilef}"
            )

    def build(self) -> str:
        func_name = "tileMetav4"
        args = [self.tilef]
        opts = utils.get_activated_options(self.options)
        flags = utils.get_activated_flags(self.options)
        return utils.to_str(func_name=func_name, args=args, options=opts, flags=flags)

    @property
    def input_files(self) -> list[Path]:
        return [Path(self.tilef)]

    @property
    def output_files(self) -> list[Path]:
        return [Path(self.tilef.replace(".mat", "_meta.txt"))]
