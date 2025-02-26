from dataclasses import dataclass
from typing import Literal

from matlib.funcs import utils

PROJECT_CHOICES = Literal["ArcticDEM", "EarthDEM", "REMA"]


@dataclass
class TileMetaV4Options:
    project: PROJECT_CHOICES | None = None
    tileVersion: str | None = None
    overwrite: bool = False

    @staticmethod
    def flag_keys() -> set[str]:
        return {"overwrite"}


@dataclass
class TileMetaV4:
    tilef: str
    options: TileMetaV4Options

    def build(self) -> str:
        func_name = "tileMetav4"
        args = [self.tilef]
        opts = utils.get_activated_options(self.options)
        flags = utils.get_activated_flags(self.options)
        return utils.to_str(func_name=func_name, args=args, options=opts, flags=flags)
