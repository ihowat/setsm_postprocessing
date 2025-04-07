import datetime as dt
from pathlib import Path

import click
import hdf5storage
import numpy as np
from numpy.typing import NDArray
import rasterio
from rasterio.transform import from_origin

Y2K_DATENUM = float(730486)
Y2K_DATETIME = dt.datetime(2000, 1, 1, 0, 0, 0, 0)


@click.command()
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(
        dir_okay=True, file_okay=False, resolve_path=True, path_type=Path, exists=True
    ),
    help="Directory where final GeoTiffs will be written",
)
@click.option(
    "--matfile",
    required=False,
    type=click.Path(
        dir_okay=False, file_okay=True, resolve_path=True, path_type=Path, exists=True
    ),
    default=None,
    help="Process one matfile",
)
@click.option(
    "--textfile",
    required=False,
    type=click.Path(
        dir_okay=False, file_okay=True, resolve_path=True, path_type=Path, exists=True
    ),
    default=None,
    help="Process a list of newline-delimited matfile paths",
)
@click.option(
    "--input-dir",
    required=False,
    type=click.Path(
        dir_okay=True, file_okay=False, resolve_path=True, path_type=Path, exists=True
    ),
    default=None,
    help="Process all .mat files in the given directory (NOT recursive).",
)
def write_gpm_to_tif(
    output_dir: Path,
    matfile: Path | None,
    textfile: Path | None,
    input_dir: Path | None,
) -> None:
    input_sources = [matfile, textfile, input_dir]
    if input_sources.count(None) == 3:
        click.echo("One of --matfile, --textfile, or --input-dir is required")
        exit(1)
    if not input_sources.count(None) == len(input_sources) - 1:
        click.echo("Only one of --matfile, --textfile, or --input-dir may be specified")
        exit(1)

    if matfile:
        inputs = [matfile]
    if textfile:
        inputs = [Path(line).resolve() for line in textfile.read_text().splitlines()]
    if input_dir:
        inputs = [path.resolve() for path in input_dir.glob("*.mat")]

    for matfile in inputs:
        click.echo(f"Processing: {matfile}")

        # Write GeoTiffs to a subdirectory named after the matfile
        subdir = output_dir / matfile.stem
        subdir.mkdir(exist_ok=True, parents=True)

        spatial_props = calc_spatial_properties(matfile)

        click.echo("Writing timeseries z GeoTiff")
        out_file = subdir / f"{matfile.stem}_z.tif"
        write_timeseries_z(matfile, out_file, spatial_props)

        click.echo("Writing timeseries zerr GeoTiff")
        out_file = subdir / f"{matfile.stem}_zerr.tif"
        write_timeseries_zerr(matfile, out_file, spatial_props)

        click.echo("Writing timeseries N GeoTiff")
        out_file = subdir / f"{matfile.stem}_N.tif"
        write_timeseries_n(matfile, out_file, spatial_props)

        click.echo("Writing date range GeoTiff")
        out_file = subdir / f"{matfile.stem}_date_range.tif"
        write_date_range(matfile, out_file, spatial_props)


def load_mat_variable(matfile: Path, variable: str) -> NDArray:
    var_dict = hdf5storage.loadmat(
        file_name=path_as_str(matfile), variable_names=[variable]
    )
    return var_dict[variable]


def load_t(matfile: Path) -> NDArray[np.datetime64]:
    arr = load_mat_variable(matfile, "t")
    arr = arr.squeeze()
    days_since_y2k = arr - Y2K_DATENUM
    datetime_array = np.array(
        [Y2K_DATETIME + dt.timedelta(days=d) for d in days_since_y2k],
        dtype=np.datetime64,
    )
    return datetime_array


def load_t0(matfile: Path) -> NDArray[np.uint16]:
    arr = load_mat_variable(matfile, "t0")
    arr = arr - Y2K_DATENUM  # Convert to days since January 1, 2000
    arr[arr < 0] = 0  # Set any negative values (days before January 1, 2000) to zero
    arr[~np.isfinite(arr)] = 0  # Set NaN, +Inf & -Inf to zero
    arr = arr.astype(np.uint16)
    return arr


def load_t1(matfile: Path) -> NDArray[np.uint16]:
    arr = load_mat_variable(matfile, "t1")
    arr = arr - Y2K_DATENUM  # Convert to days since January 1, 2000
    arr[arr < 0] = 0  # Set any negative values (days before January 1, 2000) to zero
    arr[~np.isfinite(arr)] = 0  # Set NaN, +Inf & -Inf to NoData value
    arr = arr.astype(np.uint16)
    return arr


def load_z(matfile: Path) -> NDArray[np.float32]:
    arr = load_mat_variable(matfile, "z")
    arr[~np.isfinite(arr)] = -9999  # Set NaN, +Inf & -Inf to NoData value
    return arr


def load_zerr(matfile: Path) -> NDArray[np.float32]:
    arr = load_mat_variable(matfile, "zerr")
    arr[~np.isfinite(arr)] = -9999  # Set NaN, +Inf & -Inf to NoData value
    return arr


def load_n(matfile: Path) -> NDArray[np.uint16]:
    arr = load_mat_variable(matfile, "N")
    arr = arr.astype(np.uint16)
    return arr


def load_x(matfile: Path) -> NDArray[np.float32]:
    arr = load_mat_variable(matfile, "x")
    return arr.squeeze()


def load_y(matfile: Path) -> NDArray[np.float32]:
    arr = load_mat_variable(matfile, "y")
    return arr.squeeze()


def calc_spatial_properties(
    matfile: Path, x_res: int = 100, y_res: int = 100, epsg: int = 3031
):
    x = load_x(matfile)
    y = load_y(matfile)
    height = x.shape[0]
    width = y.shape[0]
    transform = from_origin(west=x[0], north=y[0], xsize=x_res, ysize=y_res)

    return dict(
        height=height,
        width=width,
        crs=f"EPSG:{epsg}",
        transform=transform,
    )


def write_timeseries_z(matfile: Path, out_file: Path, spatial_props: dict) -> None:
    t = load_t(matfile)
    z = load_z(matfile)
    metadata = dict(
        driver="COG",
        count=t.shape[0],
        dtype=z.dtype,
        nodata=-9999,
        **spatial_props,
    )

    date_str = np.datetime_as_string(t, unit="D")
    with rasterio.open(out_file, mode="w", **metadata) as dst:
        for i, date in enumerate(date_str):
            dst.write(z[:, :, i], i + 1)
            dst.set_band_description(i + 1, date)


def write_timeseries_zerr(matfile: Path, out_file: Path, spatial_props: dict) -> None:
    t = load_t(matfile)
    zerr = load_zerr(matfile)
    metadata = dict(
        driver="COG",
        count=t.shape[0],
        dtype=zerr.dtype,
        nodata=-9999,
        **spatial_props,
    )

    date_str = np.datetime_as_string(t, unit="D")
    with rasterio.open(out_file, mode="w", **metadata) as dst:
        for i, date in enumerate(date_str):
            dst.write(zerr[:, :, i], i + 1)
            dst.set_band_description(i + 1, date)


def write_timeseries_n(matfile: Path, out_file: Path, spatial_props: dict) -> None:
    t = load_t(matfile)
    n = load_n(matfile)
    metadata = dict(
        driver="COG",
        count=t.shape[0],
        dtype=n.dtype,
        nodata=None,
        **spatial_props,
    )

    date_str = np.datetime_as_string(t, unit="D")
    with rasterio.open(out_file, mode="w", **metadata) as dst:
        for i, date in enumerate(date_str):
            dst.write(n[:, :, i], i + 1)
            dst.set_band_description(i + 1, date)


def write_date_range(matfile: Path, out_file: Path, spatial_props: dict) -> None:
    t0 = load_t0(matfile)
    t1 = load_t1(matfile)
    metadata = dict(
        driver="COG",
        count=2,
        dtype=t0.dtype,
        nodata=0,
        **spatial_props,
    )

    with rasterio.open(out_file, mode="w", **metadata) as dst:
        dst.write(t0, 1)
        dst.set_band_description(1, "mindate")
        dst.write(t1, 2)
        dst.set_band_description(2, "maxdate")


def path_as_str(p: Path | str) -> str:
    if isinstance(p, Path):
        return str(p)
    if isinstance(p, str):
        return p
    raise ValueError(f"p must be of type str or Path. Got: {type(p)}")


if __name__ == "__main__":
    write_gpm_to_tif()
