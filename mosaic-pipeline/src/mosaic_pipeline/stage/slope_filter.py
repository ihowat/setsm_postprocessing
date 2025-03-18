from pathlib import Path

from matlib.build import build_bash_command
from matlib.funcs.tile_meta_v4 import TileMetaV4, TileMetaV4Options
from matlib.funcs.add_path import AddPath
from matlib.funcs.write_tile_to_tif_v4 import WriteTileToTifV4, WriteTileToTifV4Options


def yes_slope_filter(matfile: Path) -> list[str]:
    funcs = [
        AddPath(path="/path/to/setsm_postprocessing_pgc"),
        AddPath(path="/path/to/setsm_postprocessing4"),
        WriteTileToTifV4(
            tilef=f"{matfile}",
            projstr="utm..n",
            options=WriteTileToTifV4Options(
                outFormat="COG",
                outSet="full",
                bufferMeters=100,
                registerToRef="none",
                fillWaterInterpMethod=2,
                refDemFile="/path/to/ref/dem",
                waterMaskFile="/path/to/land/cover",
                applySlopeDiffFilt=True,
            ),
        ),
        TileMetaV4(tilef=f"{matfile}", options=TileMetaV4Options()),
    ]

    return build_bash_command(funcs)
