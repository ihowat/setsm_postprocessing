from dataclasses import dataclass
from typing import Literal

from matlib.funcs import utils

OUT_FORMAT_CHOICES = Literal["LZW", "COG"]
OUT_SET_CHOICES = Literal["full", "browseOnly", "demOnly", "demAndBrowse"]
REGISTER_TO_REF_CHOICES = Literal["none", "tiles", "blobs", "reportOffsetOnly"]


@dataclass
class WriteTileToTifV4Options:
    noCrop: bool = False
    bufferMeters: int | float = 100
    overwrite: bool = False
    outFormat: OUT_FORMAT_CHOICES = "LZW"
    outSet: OUT_SET_CHOICES = "full"
    refDemFile: str | None = None
    waterMaskFile: str | None = None
    registerToRef: REGISTER_TO_REF_CHOICES = "none"
    registerToRefDebug: bool = False
    registerBlobsSkipregShp: str | None = None
    filterFillDebug: bool = False
    applySlopeDiffFilt: bool = False
    applyResidualTopographyFractionalDifferenceFilter: bool = False
    applyWaterFill: bool = False
    addSeaSurface: bool = False
    qcMaskFile: str | None = None
    fillWaterInterpMethod: int | None = 2

    @staticmethod
    def flag_keys() -> set[str]:
        return {
            "noCrop",
            "overwrite",
            "registerToRefDebug",
            "filterFillDebug",
            "applySlopeDiffFilt",
            "applyResidualTopographyFractionalDifferenceFilter",
            "applyWaterFill",
            "addSeaSurface",
        }


@dataclass
class WriteTileToTifV4:
    tilef: str
    projstr: str
    options: WriteTileToTifV4Options

    def build(self) -> str:
        func_name = "writeTileToTifv4"
        args = [self.tilef, self.projstr]
        opts = utils.get_activated_options(self.options)
        flags = utils.get_activated_flags(self.options)
        return utils.to_str(func_name=func_name, args=args, options=opts, flags=flags)
