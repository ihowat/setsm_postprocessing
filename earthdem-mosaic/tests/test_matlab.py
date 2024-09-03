from pathlib import Path

from earthdem_mosaic.pipeline.matlab import (
    AddPath,
    BatchRegisterTilesToCOP30,
    CopyFile,
    Varargin,
    build_matlab_command,
)


def test_build_matlab_command():
    funcs = [
        AddPath(Path("/path/to/matlib").as_posix()),
        CopyFile(
            src=Path("/path/to/src/file").as_posix(),
            dst=Path("/path/to/dst/file").as_posix(),
        ),
    ]
    assert build_matlab_command(funcs) == [
        "matlab",
        "-nojvm",
        "-nodisplay",
        "-nosplash",
        "-r",
        '''"try; addpath('/path/to/matlib'); copyfile('/path/to/src/file', '/path/to/dst/file'); catch e; disp(getReport(e)); exit(1); end; exit(0);"''',
    ]


def test_varargin_as_str_with_flags_and_kwargs():
    varargin = Varargin(
        flags=["flag_one", "flag_two"],
        kwargs={
            "kwarg_one": "value_one",
            "kwarg_two": "value_two",
        },
    )

    # Default prefix
    assert (
        varargin.as_str()
        == ", 'flag_one', 'flag_two', 'kwarg_one','value_one', 'kwarg_two','value_two'"
    )
    # Custom prefix
    assert (
        varargin.as_str(prefix="tada")
        == "tada'flag_one', 'flag_two', 'kwarg_one','value_one', 'kwarg_two','value_two'"
    )


def test_varargin_as_str_with_flags_only():
    varargin = Varargin(
        flags=["flag_one", "flag_two"],
    )

    # Default prefix
    assert varargin.as_str() == ", 'flag_one', 'flag_two'"
    # Custom prefix
    assert varargin.as_str(prefix="tada") == "tada'flag_one', 'flag_two'"


def test_varargin_as_str_with_kwargs_only():
    varargin = Varargin(
        kwargs={
            "kwarg_one": "value_one",
            "kwarg_two": "value_two",
        },
    )

    # Default prefix
    assert varargin.as_str() == ", 'kwarg_one','value_one', 'kwarg_two','value_two'"
    # Custom prefix
    assert (
        varargin.as_str(prefix="tada")
        == "tada'kwarg_one','value_one', 'kwarg_two','value_two'"
    )


def test_batch_register_tiles_to_cop30_with_skipreg_shp():
    cmd = BatchRegisterTilesToCOP30(
        supertile_dir=Path("/path/to/supertile/dir").as_posix(),
        ref_dem=Path("/path/to/ref_dem").as_posix(),
        ref_landcover=Path("/path/to/ref_landcover").as_posix(),
        skipreg_shp=Path("/path/to/skipreg_shp").as_posix(),
    )
    assert (
        cmd.as_str()
        == "batchRegisterTilesToCOP30('/path/to/supertile/dir', '/path/to/ref_dem', '/path/to/ref_landcover', 'registerBlobs', 'resolution','2m', 'registerBlobsSkipregShp','/path/to/skipreg_shp');"
    )

def test_batch_register_tiles_to_cop30_without_skipreg_shp():
    cmd = BatchRegisterTilesToCOP30(
        supertile_dir=Path("/path/to/supertile/dir").as_posix(),
        ref_dem=Path("/path/to/ref_dem").as_posix(),
        ref_landcover=Path("/path/to/ref_landcover").as_posix(),
    )
    assert (
        cmd.as_str()
        == "batchRegisterTilesToCOP30('/path/to/supertile/dir', '/path/to/ref_dem', '/path/to/ref_landcover', 'registerBlobs', 'resolution','2m');"
    )
