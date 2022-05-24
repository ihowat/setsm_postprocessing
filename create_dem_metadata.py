# -*- coding: utf-8 -*-
"""
Extract polygon boundary from the metadata textfile of REMA, ArcticDEM, and EarthDEM
May 06, 2020: svalrd ArcticDEMs inside the EarthDEM folder; using argparser
Last Run: March 03, 2021 for arcticDEM
Sep 17, 2020 [Modified for EarthDEM]
Mar 31, 2022 [Updating for v4.1 strips ArcticDEM]
May 12, 2022 [Dry run for EarthDEM v4 and seems working]
@author: yadav.111@osu.edu

Manually, change/update this script to choose folder location, dem type, resolution before running

TODO:
=====
    Create a log file
    What to do with missing data? ie, there is folder for strip, but not (dem) data inside the strip folder
        Ignore them as they are just a placeholder when (good?) DEM was not created
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
# from rema_dem import get_rema_strip_polygon
from shapely.geometry import Polygon

import argparse
import logging


def extract_polygon(meta_file):
    ''' Extract coordinates, creation date, rmse etc from txt file for each DEM by reading the metadata text file
        Inputs
        ======
        meta_file: Full path to metadata text file
    '''
    with open(meta_file, 'r') as infile:
        data_dict = {}  # use dictionary to collect and organize attributes to go into geopandas dataframe
        line = ''  # just a starting point for variable to use in while loop
        line = infile.readline()
        dt1 = []
        dt2 = []
        # time_flag1 = ''
        # time_flag2 = ''
        # while 'Scene Metadata' not in line:
        while line:
            # if "Creation Date" in line: #perhaps not required, this is just day algorithm was run - has no meaning for metadata
            #     creationDate = line
            if 'Strip projection (proj4)' in line:
                proj4 = line.split(':')[-1]
                # Strip new line character, ', then extra space at end
                proj4 = proj4.strip().strip("""'""").strip()

            if 'X:' in line:
                X = line.split(':')[1].strip()  # second item is a string of X's
                X = np.fromstring(X, sep=' ')   # convert string to numpy array
            if 'Y:' in line:
                Y = line.split(':')[1].strip()
                Y = np.fromstring(Y, sep=' ')
            # New : Get the RMSE: can have one or more scenes with end denoted by blank line followed by "Filtering Applied:"
            if 'scene, rmse, dz, dx, dy, dz_err, dx_err, dy_err' in line:
                rmse_list = []
                line = infile.readline()
                while line:
                    rmse_val = line.split()[1]
                    rmse_list.append(rmse_val)
                    line = infile.readline()
                    if len(line.strip()) == 0:
                        break  # break the while loop
                # These list comprehension works for empty list as well
                rmse_list = [float(v) for v in rmse_list if v!='nan']  # remove nans and convert to float
                rmse_list = [v for v in rmse_list if v>0]  # remove all zeros
                if len(rmse_list) > 0:
                    rmse = np.array(rmse_list)
                    rmse = np.mean(rmse)
                else:
                    rmse = -1  # meaning it was not calculated
            # Append new version of mean time1 and time2 (Ian Howat Aug 08 2020)
            # for line in lines:
            if 'Image_1_Acquisition_time' in line:
                dt = line.split('=')[-1]
                dt = pd.to_datetime(dt)
                dt1.append(dt)
                # time_flag1 = 'new'
            if 'Image_2_Acquisition_time' in line:
                dt = line.split('=')[-1]
                dt = pd.to_datetime(dt)
                dt2.append(dt)
                # time_flag2 = 'new'
            line = infile.readline()  # read the next line
        # Cacluate average acquisition times for images
        if len(dt1) > 0:
            ser1 = pd.Series(dt1)
            time1 = ser1.mean()  # perhaps min is better?
            time1 = time1.strftime("%Y-%m-%d %H:%M:%S")
        else:
            # time1 = ''
            # TODO: Can call the function for reading old file directly here
            time1 = fill_missing_time(meta_file, img_time_string='Image 1=')
            # time_flag1 = 'old'
        if len(dt2) > 0:
            ser2 = pd.Series(dt2)
            time2 = ser2.mean()
            time2 = time2.strftime("%Y-%m-%d %H:%M:%S")
        else:
            # time2 = ''
            time2 = fill_missing_time(meta_file, img_time_string='Image 2=')
            # time_flag2 = 'old'
        # Collect data in a dictionary
        # To use numpy array as a cell inside a dataframe, Use a wrapper around the numpy array i.e. pass the numpy array as list
        # https://stackoverflow.com/questions/45548426/store-numpy-array-in-cells-of-a-pandas-dataframe
        data_dict['x'] = list(X)
        data_dict['y'] = list(Y)  # for pandas to export to csv
        data_dict['time1'] = time1
        data_dict['time2'] = time2
        data_dict['rmse'] = rmse
        # data_dict['time_flag1'] = time_flag1
        # data_dict['time_flag2'] = time_flag2

        # Create dataframe and finally a geodataframe
        df = pd.DataFrame([data_dict])
        geom = []
        for i in range(len(df)):
            # This can be made easier or more legible
            geom.append(Polygon(zip(df.x.values[i], df.y.values[i])))
        df['geometry'] = geom
        # gdf = gpd.GeoDataFrame(df[['rmse', 'time1', 'time2', 'time_flag1', 'time_flag2', 'geometry']], geometry='geometry', crs = proj4)
        gdf = gpd.GeoDataFrame(df[['rmse', 'time1', 'time2', 'geometry']], geometry='geometry', crs = proj4)
        # gdf.crs = proj4  # Same crs for all tiles
        # gdf = gdf.to_crs('EPSG:4326') # to change to lat/lon geographic crs
        return gdf


def fill_missing_time(meta_file, img_time_string):
    ''' Fill in the time1 and time2 for older version tiles
        ALERT: This is hard-coded to REMA (see meta_file path below)

        Inputs
        ======
        meta_file: Full path to metadata text file
        img_time_string: 'Image 1=' or 'Image 2='
            image acquisition time will be extracted from this substring
            only for older version of metadata
    '''
    # meta_file = f'{dem_folder}/{region}/strips_v4  /2m/{strip}_{version}/{strip}_seg{seg_id}_meta.txt'
    # meta_file = f'{dem_folder}/{region}/strips_v4.1/2m/{strip}_{version}/{strip}_seg{seg_id}_meta.txt'
    dt1 = []
    with open(meta_file, 'r') as infile:
        for line in infile.readlines():
            # if 'Image 1=' in line:
            if img_time_string in line:
                file_parts = line.split('/')
                img1 = file_parts[-1]
                img1 = img1.split('.tif')[0]  # to remove extension (.tif) from name
                # From this img1 (ie part of filename, extract the datetime of image acquisition)
                dt = img1.split('_')[1]  # sample '20091120125025' [YYYYMMDDhhmmss]
                dt = pd.to_datetime(dt)
                dt1.append(dt)
    # Cacluate average acquisition times for images
    if len(dt1) > 0:
        ser1 = pd.Series(dt1)
        time1 = ser1.mean()
        time1 = time1.strftime("%Y-%m-%d %H:%M:%S")
        return time1


def get_strip_polygon(strip_folder):
    ''' Mar 22, 2022: Renamed from get_rema_strip_polygon to incorporate processing the v4.1 DEMs
        Parses a REMA folder to reveal if one or more tiles are present in  folder
        strip_folder: Full path to 2m new REMA strips
        Assumptions:
            strip folder name always end with _v0xxxx
            Eg: W1W1_20081215_1020010005DB7A00_1020010005DE1800_2m_lsf_v040000
    '''
    # strip the trailing version number from folder names because actual DEM,
    # matchtag etc inside the folder does not have verison info
    # /fs/byo/howat-data5/ArcticDEM/arcticdem_03_greenland_southwest/strips_v4.1/2m/WV01_20180827_10200100771C9A00_10200100783E1000_2m_lsf_v030400
    strip = strip_folder.split('/')[-1]  # extract folder name; it has info in strip version
    strip = strip.split('_')
    strip_version = strip[-1]
    strip = '_'.join(strip[:-1])

    files = os.listdir(strip_folder)
    metas = [f for f in files if f.endswith('_meta.txt')]
    # sort in order of segmennts: 1, 2, 3, ..., 10, 11 etc.
    metas = sorted(metas, key=lambda x: int(x.split('_seg')[-1].split('_')[0]))  # no error even if list is empty
    if len(metas) == 0:
        # To Guard aginst empty/missing data inside the strip folder; but revise again; can cause other error downstream
        logging.error(f'Missing metadata: {strip_folder}')
        return
    meta_count = 0
    for meta in metas:
        meta_file = f'{strip_folder}/{meta}'  # fullpath to metadata text file
        seg_id = meta.split('_seg')[-1].split('_')[0]
        if meta_count == 0:
            gdf = extract_polygon(meta_file)
            gdf['seg_id'] = seg_id
        else:
            gdf_temp = extract_polygon(meta_file)
            gdf_temp['seg_id'] = seg_id
            gdf = pd.concat([gdf, gdf_temp])
        meta_count += 1
    # Perhaps dissolve to simplify the geometry if more than one tiles present within a strip
    # if len(gdf)>1:
    #    gdf = gdf.dissolve(by='name')
    gdf['strip'] = strip
    gdf['version'] = strip_version # Ian wanted in newer version for metadata [Aug 05, 2020]
    # gdf = gdf[['strip', 'version', 'seg_id', 'time1', 'time2', 'time_flag1', 'time_flag2', 'rmse', 'geometry']]  # Reorganize for presentation
    gdf = gdf[['strip', 'version', 'seg_id', 'time1', 'time2', 'rmse', 'geometry']]  # Reorganize for presentation
    return gdf


def main():
    # dem_type = sys.argv[1]  # pass REMA or ArcticDEM as an argument from the calling function
    parser = argparse.ArgumentParser(description='Create Shapefile for DEM.')
    parser.add_argument('dem_type', help='REMA, ArcticDEM, EarthDEM', type=str)
    parser.add_argument('--input', help='Input folder for DEMs', type=str)
    parser.add_argument('--output', help='Output folder for metadata shapefile', type=str)
    parser.add_argument('--log_name', help='Name of Log file', type=str, default='dem_metadata.log')

    args = parser.parse_args()
    dem_type = args.dem_type
    dem_folder = args.input
    metadata_folder = args.output
    log_name = args.log_name

    logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
    logging.info('-------------------------START LOGGING--------------------------------------')
    logging.info(f'dem_type: {dem_type}   dem_folder = {dem_folder}  metadata_folder = {metadata_folder}')

    if dem_type == 'REMA':
        dir_prefix = '/fs/byo/howat-data4'  # Main location where dems are staged
        prefix = 'rema'  # prefix for a given dem
        dem_folder = f'{dir_prefix}/{dem_type}'  # Path to DEM parent folder containing strips
        metadata_folder = f'{dem_folder}/metadata_{prefix}'  # To save metadata shapefile
        # List regions: subfolders in each region contain different segments of a DEM
        regions = os.listdir(dem_folder)
        regions = [f for f in regions if f.startswith(f'{prefix}_')]
        regions = [f for f in regions if not f.endswith('.mat')]  # TODO: better to check folder rather than individual files
        regions = [f for f in regions if not f.endswith('.txt')]
        regions = [f for f in regions if not f.endswith('.tif')]
        regions = [f for f in regions if not '_strip_' in f]
        strip_version = "v4"
        logging.info(f'Regions = {regions}')

    elif dem_type == 'ArcticDEM':
        dir_prefix = '/fs/byo/howat-data5'  # Example for ArcticDEM: /fs/byo/howat-data5/ArcticDEM/arcticdem_02_greenland_southeast/strips_v4/2m/WV01_20200416_10200100992F7600_10200100961F1200_2m_lsf_v040203
        # May 06, 2021: The 10m strips are now ready on Unity here:  /fs/project/howat.4/EarthDEM/arcticdem_14_svalbard
        prefix = 'arcticdem'
        dem_folder = f'{dir_prefix}/{dem_type}'  # Each strip inside this folder
        # Note: only for ArciticDEM create metadata in the base_folder because I do not have write priviledge inside the ArcticDEM folder
        # metadata_folder = f'{dir_prefix}/testing/metadata_{prefix}_v41'
        metadata_folder = f'{dir_prefix}/metadata_{prefix}_v41'
        # List regions: subfolders in each region contain different segments of a DEM
        regions = os.listdir(dem_folder)
        regions = [f for f in regions if f.startswith(f'{prefix}_')]
        regions = [f for f in regions if not f.endswith('.mat')]  # right now not required for arcticDEM
        regions = [f for f in regions if not f.endswith('_qc')]
        regions = [f for f in regions if not f.endswith('_shades')]
        # regions = ['arcticdem_02_greenland_southeast']  # Aside: for testing new version only
        strip_version = "v4.1"

    elif dem_type == 'EarthDEM':
        print('Creating metadata for EarthDEM')
        dir_prefix = '/fs/project/howat.4'  # Choose the main location where dems are staged by region
        prefix = 'EarthDEM'  #svalbard EarthDEM  # prefix for a given dem
        dem_folder = f'{dir_prefix}/{dem_type}'  # Path to DEM parent folder containing strips
        # dem_folder = f'{dir_prefix}/{dem_type}/{region}/strips_unf/2m' # EarthDEM
        metadata_folder = f'{dem_folder}/metadata_{prefix}'  # Shapefile metadata will be saved inside this folder
        # .../howat.4/EarthDEM/region_31_alaska_south/strips_v4/2m
        regions = ['region_31_alaska_south', 'region_34_alaska_north'] # 'region_01_iceland',  directly
        # regions = ['arcticdem_14_svalbard']
        strip_version = "v4"
    else:
        # Raise Error
        assert(False, 'Must pass either REMA, ArcticDEM, or EarthDEM')
        sys.exit(1)

    if not os.path.exists(metadata_folder):
        # Create Output folder to save metadata
        os.makedirs(metadata_folder)
    # -----------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------
    for region in regions:
        # First list all the dem strips
        # Mar 22, 2022: we have strips_v4.1 now
        # strip_subfolder = 'strips_v4.1'  # strips_v4.1 strips_v4 : name of the subfolder
        strip_subfolder = f"strips_{strip_version}"  # May 12, 2022
        strips_folder = f'{dem_folder}/{region}/{strip_subfolder}/2m'  # can be dangerous to make this assumption
        strips = os.listdir(f'{strips_folder}')
        strips = [f for f in strips if not f.endswith('.json')]  # Mar 30, 2022 (for new version of DEMs)
        logging.info(f'Number of Strips for {region} = {len(strips)}')
        count = 0
        for strip in strips:
            # strip_folder = f'{dir_prefix}/{dem_type}/{region}/strips_unf/2m/{strip}'
            strip_folder = f'{dem_folder}/{region}/{strip_subfolder}/2m/{strip}'
            if count == 0:
                # gdf = get_rema_strip_polygon(strip_folder)  #, option='hma'
                gdf = get_strip_polygon(strip_folder)  # New for strip version v4.1
            else:
                # gdf_temp = get_rema_strip_polygon(strip_folder)  #, option='hma'
                gdf_temp = get_strip_polygon(strip_folder)  # New for strip version v4.1
                # Concatenate
                gdf = pd.concat([gdf, gdf_temp])
            count += 1
        gdf.index = range(len(gdf))  # Required: assingn index, else subsetting for missing won't work

        # # Commented on March 30, 2022 since missing datetime assumed to be non-issue in new version v4.1
        # # Check if any row is missing time1 time2 due to strips being in old format
        # # subset = gdf[gdf.time1.isnull()] #when saved shp file is open, empty string to read as null
        # subset = gdf[gdf.time1=='']  # Since the shp is not yet saved, it still contains empty string
        # if len(subset) > 0:
        #     '''This ad-hoc code is only for old files. Remove once only new files are processed PGC'''
        #     logging.info(f'Missing time {len(subset)}')
        #     t = subset.apply(fill_missing_time, args=(region,), axis=1)
        #     gdf.loc[t.index, 'time1'] = t.apply(lambda x: x[0])
        #     gdf.loc[t.index, 'time2'] = t.apply(lambda x: x[1])
        #     # gdf.loc[t.index, 'flag'] = 'old'
        gdf.to_file(f'{metadata_folder}/{region}.shp')
        # gdf.to_file(f'{metadata_folder}/{region}.gpkg', driver='GPKG')
    # ------------------------------------------------------------------------------------------------------------

    shp_files = os.listdir(metadata_folder)
    # shp_files = [f for f in shp_files if (f.endswith('.shp') and f.startswith('region_'))]
    # shp_files = [f for f in shp_files if (f.endswith('.shp') and f.startswith('rema_'))] #REMA
    # Need to make the listing of shapefiles more complex to ignore rema_regions.shp and rema_strips.shp in subsequent processing
    # shp_files = [f for f in shp_files if (f.endswith('.shp') and f.startswith(f'{prefix}_'))] #ArcticDEM
    shp_files = [f for f in shp_files if (f.endswith('.shp'))]
    shp_files = [f for f in shp_files if not f.endswith('_regions.shp')]
    shp_files = [f for f in shp_files if not f.endswith('_strips.shp')]

    # ------------------------------------------------------------------------------------------------------------
    # Part2: Get Outlines of each region and Merge into Pan-Antarctic Map
    # Derive Outline of shapefile (by dissolving) and concatenating region outlines
    count = 0
    for shp_file in shp_files:
        gdf = gpd.read_file(f'{metadata_folder}/{shp_file}')
        # Dissolve Polygons within a shapefile to get only the Outline
        # Need to add some attribute based on which to dissolve
        gdf['region'] = shp_file.split('.')[0]#shp_file[:9]
        if count == 0:
            # No concatenation
            outline_gdf = gdf.dissolve(by='region')  ## This is a geodataframe with single row
            # Attach more metadata from the original shapefile
            # outline_gdf[['xmax', 'xmin', 'ymax', 'ymin']] = gdf.xmax.max(), gdf.xmin.min(), gdf.ymax.max(), gdf.ymin.min()
            # outline_gdf['Nscenes'] = gdf.Nscenes.sum()
        else:
            outline2 = gdf.dissolve(by='region')
            # outline2[['xmax', 'xmin', 'ymax', 'ymin']] = gdf.xmax.max(), gdf.xmin.min(), gdf.ymax.max(), gdf.ymin.min()
            # outline2['Nscenes'] = gdf.Nscenes.sum()
            outline_gdf = pd.concat([outline_gdf, outline2])  # Concatenate
        count += 1
    # Remove Columns that does not have meaning in merged product
    outline_gdf = outline_gdf.drop(['version', 'seg_id', 'rmse', 'strip'], axis=1)
    # Add more attributes after merge Area, Perimeter
    outline_gdf['Area'] = np.round(outline_gdf.area/10**6)  # round t0 1 or zero decimals
    outline_gdf['Perim'] = np.round(outline_gdf.length/10**3)
    outline_gdf = outline_gdf[['Area', 'Perim', 'geometry']]
    outline_gdf.to_file(f'{metadata_folder}/{prefix}_regions.shp')
    # outline_gdf.to_file(f'{metadata_folder}/{prefix}_regions.geojson', driver='GeoJSON')

    # ------------------------------------------------------------------------------------------------------------
    # Part 3: Merge all strips into one Pan-strip-dem
    count = 0
    for shp_file in shp_files:
        if count == 0:
            gdf = gpd.read_file(f'{metadata_folder}/{shp_file}')
            gdf['region'] = shp_file.split('.')[0]  # Attach region column
        else:
            gdf2 = gpd.read_file(f'{metadata_folder}/{shp_file}')
            gdf2['region'] = shp_file.split('.')[0]  # Attach region column
            gdf = pd.concat([gdf, gdf2])  # Concatenate
        count += 1
    # Save shapefile
    gdf.to_file(f'{metadata_folder}/{prefix}_strips.shp')
    # gdf.to_file(f'{metadata_folder}/{prefix}_strips.geojson', driver='GeoJSON')

# ------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
    print('Finished Script')

"""
OLD Codes: Prior to March 31, 2022
================================== 
def fill_missing_time_x(row, region):
    ''' Fill in the time1 and time2 for older tiles 
        ALERT: This is hard-coded to REMA (see meta_file path below)
        After March 31, 2022: this is not being used
    '''
    # region = row['region']
    strip = row['strip']
    version = row['version']
    seg_id = row['seg_id']
    # meta_file = f'{dem_folder}/{region}/strips_v4/2m/{strip}_{version}/{strip}_seg{seg_id}_meta.txt'
    meta_file = f'{dem_folder}/{region}/strips_v4.1/2m/{strip}_{version}/{strip}_seg{seg_id}_meta.txt'
    dt1 = []
    dt2 = []
    with open(meta_file, 'r') as infile:
        for line in infile.readlines():
            if 'Image 1=' in line:
                file_parts = line.split('/')
                # region1 = '/'.join(file_parts[1:-3]) #file_parts[-3]
                img1 = file_parts[-1]
                img1 = img1.split('.tif')[0] #to remove extension (.tif) from name
                # From this img1 (ie part of filename, extract the datetime of image acquisition)
                dt = img1.split('_')[1] # sample '20091120125025' [YYYYMMDDhhmmss]
                dt = pd.to_datetime(dt)
                dt1.append(dt)
                # this can be readily converted datetime, eg pd.to_datetime('20091120125025') => Timestamp('2009-11-20 12:50:25')
            if 'Image 2=' in line:
                file_parts = line.split('/')
                img2 = file_parts[-1]
                img2 = img2.split('.tif')[0]
                dt = img2.split('_')[1]
                dt = pd.to_datetime(dt)
                dt2.append(dt)
    # Cacluate average acquisition times for images
    if len(dt1) > 0:
        ser1 = pd.Series(dt1)
        time1 = ser1.mean()
        time1 = time1.strftime("%Y-%m-%d %H:%M:%S")
    if len(dt2) > 0:
        ser2 = pd.Series(dt2)
        time2 = ser2.mean()
        time2 = time2.strftime("%Y-%m-%d %H:%M:%S")
    return time1, time2


def get_rema_strip_polygon_x(strip_folder, option=None):
    ''' Use this only for older DEMS (Mar 22, 2022)
        After March 31, 2022: this is not being used

        Parses a REMA folder to reveal if one or more tiles are present in  folder
        strip_folder: Full path to 2m new REMA strips
        option: to match DEMs that do not match the default criteria of DEM ; eg hma folder structure was different
        Assumptions:
            strip folder name always end with _v0xxxx
            Eg: W1W1_20081215_1020010005DB7A00_1020010005DE1800_2m_lsf_v040000
    '''
    strip = strip_folder.split('/')[-1]
    # strip the trailing version number from folder names because actual DEM, 
    # matchtag etc inside the folder does not have verison info
    strip = strip.split('_')
    strip_version = strip[-1]
    strip = '_'.join(strip[:-1])

    files = os.listdir(strip_folder)
    metas = [f for f in files if f.endswith('_meta.txt')]
    if len(metas) == 0:
        # To Guard aginst empty/missing data inside the strip folder; but revise again; can cause other error downstream
        logging.error(f'Missing metadata: {strip_folder}')
        return
    for seg_id in range(1, len(metas)+1):
        # f'{dir_prefix}/REMA/{region}/strips_unf/2m/{strip}/{strip}_seg{seg_id}_meta.txt'
        meta_file = f'{strip_folder}/{strip}_seg{seg_id}_meta.txt'  # This may be wrong; though overwritten used use to if/else condition below
        # Quick fix for HMA older files
        if option == 'hma':
            meta_file = f'{strip_folder}/{metas[seg_id-1]}'
        if seg_id == 1:
            gdf = extract_polygon(meta_file)
            gdf['seg_id'] = seg_id
        else:
            gdf_temp = extract_polygon(meta_file)
            gdf_temp['seg_id'] = seg_id
            gdf = pd.concat([gdf, gdf_temp])
    # Perhaps dissolve to simplify the geometry if more than one tiles present within a strip
    # if len(gdf)>1:
    #    gdf = gdf.dissolve(by='name')
    gdf['strip'] = strip
    gdf['version'] = strip_version # Ian wanted in newer version for metadata [Aug 05, 2020]
    gdf = gdf[['strip', 'version', 'seg_id', 'time1', 'time2', 'rmse', 'geometry']] #Reorganize for presentation
    return gdf


"""
