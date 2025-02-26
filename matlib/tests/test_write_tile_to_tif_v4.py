from matlib.funcs.write_tile_to_tif_v4 import (
    WriteTileToTifV4, WriteTileToTifV4Options,
)
from matlib.funcs.utils import get_activated_flags, get_activated_options


def test_activated_flags():
    options = WriteTileToTifV4Options(noCrop=True, overwrite=True)
    assert {"overwrite", "noCrop"} == get_activated_flags(options)


def test_default_activated_options():
    options = WriteTileToTifV4Options()
    expected = {
        "bufferMeters": 100,
        "outFormat": "LZW",
        "outSet": "full",
        "registerToRef": "none",
        "fillWaterInterpMethod": 2,
    }
    assert expected == get_activated_options(options)


def test_activated_options():
    options = WriteTileToTifV4Options(refDemFile="/path/to/cop30")
    expected = {
        "bufferMeters": 100,
        "outFormat": "LZW",
        "outSet": "full",
        "registerToRef": "none",
        "fillWaterInterpMethod": 2,
        "refDemFile": "/path/to/cop30",
    }
    assert expected == get_activated_options(options)


def test_write_tile_to_tif_v4():
    options = WriteTileToTifV4Options(overwrite=True, refDemFile="/path/to/cop30")
    func = WriteTileToTifV4(
        tilef="/path/to/18_20_1_1.mat", projstr="polar stereo north", options=options
    )
    expected = "writeTileToTifv4('/path/to/18_20_1_1.mat', 'polar stereo north', 'bufferMeters',100, 'outFormat','LZW', 'outSet','full', 'refDemFile','/path/to/cop30', 'registerToRef','none', 'fillWaterInterpMethod',2, 'overwrite');"
    assert expected == func.build()
