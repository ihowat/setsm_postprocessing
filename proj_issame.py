
import os
import sys
import traceback
import warnings

from osgeo import gdal, osr


## GDAL error handler setup
def GdalErrorHandler(err_level, err_no, err_msg):
    error_message = (
        "Caught GDAL error (err_level={}, err_no={}) "
        "where level >= gdal.CE_Warning({}); error message below:\n{}".format(
            err_level, err_no, gdal.CE_Warning, err_msg
        )
    )
    if err_level == gdal.CE_Warning:
        warnings.warn(error_message)
    if err_level > gdal.CE_Warning:
        raise RuntimeError(error_message)
gdal.PushErrorHandler(GdalErrorHandler)
gdal.UseExceptions()


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


items = sys.argv[1:]

srs_base = None
try:
    for item in items:

        srs_compare = osr.SpatialReference()

        if os.path.isfile(item):
            rasterFile = item
            raster_ds = gdal.Open(rasterFile, gdal.GA_ReadOnly)
            srs_compare.ImportFromWkt(raster_ds.GetProjectionRef())

        elif item.startswith('EPSG:') or item.isdigit():
            epsg_code = item
            srs_compare.ImportFromEPSG(int(epsg_code.replace('EPSG:', '')))

        else:
            proj4_str = item
            srs_compare.ImportFromProj4(proj4_str)

        if srs_base is None:
            srs_base = srs_compare
            continue

        compare_result = srs_compare.IsSame(srs_base)
        if compare_result == 0:
            # Not all projections are the same, return 1
            sys.exit(1)
        elif compare_result != 1:
            raise RuntimeError("ERROR: Result of IsSame in proj_issame.py ({}) is not 0 or 1".format(compare_result))

except Exception as e:
    eprint("CAUGHT EXCEPTION in proj_issame.py, returning with exit status 2")
    eprint("\nScript arguments were as follows:\n  \"{}\"".format('" "'.join(sys.argv)))
    eprint("\nError trace is as follows:\n")
    traceback.print_exc()
    sys.exit(2)


# All projections are the same, return 0
sys.exit(0)
