from matlib.funcs.write_tile_to_tif_v4 import (
    WriteTileToTifV4Options,
    WriteTileToTifV4,
)
from matlib.funcs.add_path import AddPath
from matlib.funcs.tile_meta_v4 import TileMetaV4, TileMetaV4Options
from matlib.build import build_matlab_script, build_bash_command


def test_build_matlab_script():
    tilef = "/path/to/utm02s_53_02_1_1_2m.mat"
    projstr = "utm02s"
    tile_to_tif_opts = WriteTileToTifV4Options(
        outFormat="COG",
        outSet="full",
        bufferMeters=100,
        registerToRef="none",
        fillWaterInterpMethod=2,
        refDemFile="/path/to/cop30.tif",
        waterMaskFile="/path/to/worldcover.tif",
    )

    result = build_matlab_script(
        [
            AddPath("/path/to/setsm_postprocessing_pgc"),
            AddPath("/path/to/setsm_postprocessing4"),
            WriteTileToTifV4(tilef, projstr, tile_to_tif_opts),
            TileMetaV4(tilef, TileMetaV4Options()),
        ]
    )
    expected = "try; addpath('/path/to/setsm_postprocessing_pgc'); addpath('/path/to/setsm_postprocessing4'); writeTileToTifv4('/path/to/utm02s_53_02_1_1_2m.mat', 'utm02s', 'bufferMeters',100, 'outFormat','COG', 'outSet','full', 'refDemFile','/path/to/cop30.tif', 'waterMaskFile','/path/to/worldcover.tif', 'registerToRef','none', 'fillWaterInterpMethod',2); tileMetav4('/path/to/utm02s_53_02_1_1_2m.mat');; catch e; disp(getReport(e)); exit(1); end; exit(0);"
    assert expected == result


def test_build_bash_command():
    result = build_bash_command(
        [
            AddPath("/path/to/setsm_postprocessing_pgc"),
            AddPath("/path/to/setsm_postprocessing4"),
        ]
    )
    expected = [
        "matlab",
        "-nojvm",
        "-nodisplay",
        "-r",
        "\"try; addpath('/path/to/setsm_postprocessing_pgc'); addpath('/path/to/setsm_postprocessing4');; catch e; disp(getReport(e)); exit(1); end; exit(0);\"",
    ]
    assert expected == result
