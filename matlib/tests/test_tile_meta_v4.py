from matlib.funcs.tile_meta_v4 import TileMetaV4Options, TileMetaV4
from matlib.funcs.utils import get_activated_flags, get_activated_options


def test_activated_flags():
    options = TileMetaV4Options(overwrite=True)
    assert {"overwrite"} == get_activated_flags(options)


def test_default_activated_options():
    options = TileMetaV4Options()
    assert {} == get_activated_options(options)


def test_activated_options():
    options = TileMetaV4Options(project="EarthDEM", tileVersion="1.1")
    expected = {"project": "EarthDEM", "tileVersion": "1.1"}
    assert expected == get_activated_options(options)


def test_write_tile_to_tif_v4():
    options = TileMetaV4Options(project="EarthDEM", tileVersion="1.1")
    func = TileMetaV4(tilef="/path/to/18_20_1_1.mat", options=options)
    expected = "tileMetav4('/path/to/18_20_1_1.mat', 'project','EarthDEM', 'tileVersion','1.1');"
    assert expected == func.build()
