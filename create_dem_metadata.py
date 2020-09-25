# -*- coding: utf-8 -*-
"""
Extract polygon boundary from the metadata textfile of REMA, ArcticDEM, and EarthDEM
Sep 17, 2020 [Modified for EarthDEM]
@author: yadav.111

Change this script to choose folder location, dem type, resolution before running
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
# from rema_dem import get_rema_strip_polygon
from shapely.geometry import Polygon

dem_type = sys.argv[1] # pass REMA or ArcticDEM as an argument from the calling function

if dem_type == 'REMA':
    dir_prefix = '/fs/project/howat.4' # Choose the main location where dems are staged by region
    prefix = 'rema' # prefix for a given dem
    dem_folder = f'{dir_prefix}/{dem_type}' # Path to DEM parent folder containing strips 
    metadata_folder = f'{dem_folder}/metadata_{prefix}' # Shapefile metadata will be saved inside this folder
    # List regions: subfolders in each region contain different segments of a DEM
    regions = os.listdir(dem_folder)
    regions = [f for f in regions if f.startswith(f'{prefix}_')]
    regions = [f for f in regions if not f.endswith('.mat')] # right now not required for arcticDEM
elif dem_type == 'ArcticDEM':
    dir_prefix = '/fs/byo/howat-data5' # Example for ArcticDEM: /fs/byo/howat-data5/ArcticDEM/arcticdem_02_greenland_southeast/strips_v4/2m/WV01_20200416_10200100992F7600_10200100961F1200_2m_lsf_v040203
    prefix = 'arcticdem'
    dem_folder = f'{dir_prefix}/{dem_type}' # Each strip inside this folder
    metadata_folder = f'{dir_prefix}/metadata_{prefix}' # 2. only for ArciticDEM because I do not have write priviledge inside the folder
    # List regions: subfolders in each region contain different segments of a DEM
    regions = os.listdir(dem_folder)
    regions = [f for f in regions if f.startswith(f'{prefix}_')]
    regions = [f for f in regions if not f.endswith('.mat')] # right now not required for arcticDEM
elif dem_type == 'EarthDEM':
    print('Creating metadata for EarthDEM')
    dir_prefix = '/fs/project/howat.4' # Choose the main location where dems are staged by region
    prefix = 'EarthDEM' # prefix for a given dem
    dem_folder = f'{dir_prefix}/{dem_type}' # Path to DEM parent folder containing strips 
    #dem_folder = f'{dir_prefix}/{dem_type}/{region}/strips_unf/2m' # EarthDEM
    metadata_folder = f'{dem_folder}/metadata_{prefix}' # Shapefile metadata will be saved inside this folder
    #.../howat.4/EarthDEM/region_31_alaska_south/strips_v4/2m
    regions = ['region_01_iceland', 'region_31_alaska_south', 'region_34_alaska_north'] # directly
else:
    # Raise Error
    assert(False, 'Must pass either REMA or ArcticDEM')
    sys.exit(1)

if not os.path.exists(metadata_folder):
    # Create Output folder to save metadata
    os.makedirs(metadata_folder)
#------------------------------------------------------------------------------------------------------------

def extract_rema_polygon(meta_file):
    ''' Extract coordinates, creation date, rmse etc from txt file for each REMA DEM
        
        There some some missing strips for 2m DEM, hence, this should be used only for 8m DEMs
        Inputs:
            dir_prefix, txt_file
            proj = rema or arctic
    '''
    # With newer files we will pass the full path to meta_file directly
    with open(meta_file, 'r') as infile:
        data_dict = {} #use dictionary to organize everything
        line = '' # just a starting point for variable to use in while loop [infile.readline()]
        line = infile.readline()
        dt1 = []
        dt2 = []
        while line:
        #while 'Scene Metadata' not in line:
            # if "Creation Date" in line: #perhaps not required, this is just the day algorithm was run - has no meaning for metadata
            #     creationDate = line
            if 'Strip projection (proj4)' in line:
                proj4 = line.split(':')[-1]
                # Strip new line character, ', then extra space at end
                proj4 = proj4.strip().strip("""'""").strip()

            if 'X:' in line:
                X = line.split(':')[1].strip() # second item is a string of X's
                X = np.fromstring(X, sep=' ')  # convert string to numpy array
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
                        break #break the while loop
                # These list comprehension works for empty list as well
                rmse_list = [float(v) for v in rmse_list if v!='nan'] #remove nans and convert to float
                rmse_list = [v for v in rmse_list if v>0] # remove all zeros
                if len(rmse_list) > 0:
                    rmse = np.array(rmse_list)
                    rmse = np.mean(rmse)
                else:
                    rmse = -1 #
            # Append new version of mean time1 and time2 (Ian Howat Aug 08 2020)
            # for line in lines:
            if 'Image_1_Acquisition_time' in line:
                dt = line.split('=')[-1]
                dt = pd.to_datetime(dt)
                dt1.append(dt)
            if 'Image_2_Acquisition_time' in line:
                dt = line.split('=')[-1]
                dt = pd.to_datetime(dt)
                dt2.append(dt)
            line = infile.readline() #Read the next line
        # Cacluate average acquisition times for images
        if len(dt1)>0:
            ser1 = pd.Series(dt1)
            time1 = ser1.mean()
            time1 = time1.strftime("%Y-%m-%d %H:%M:%S")
        else:
            time1 = ''
        if len(dt2)>0:
            ser2 = pd.Series(dt2)
            time2 = ser2.mean()
            time2 = time2.strftime("%Y-%m-%d %H:%M:%S")
        else:
            time2 = ''
        # Collect data in a dictionary
        # To use numpy array as a cell inside a dataframe, Use a wrapper around the numpy array i.e. pass the numpy array as list
        # https://stackoverflow.com/questions/45548426/store-numpy-array-in-cells-of-a-pandas-dataframe
        data_dict['x'] = list(X)
        data_dict['y'] = list(Y) # for pandas to export to csv
        data_dict['time1'] = time1
        data_dict['time2'] = time2
        data_dict['rmse'] = rmse

        # And creat dataframe and geodataframe
        df = pd.DataFrame([data_dict])
        geom = []
        for i in range(len(df)):
            # This can be made easier or more legible
            geom.append(Polygon(zip(df.x.values[i], df.y.values[i])))
        df['geometry'] = geom
        gdf = gpd.GeoDataFrame(df[['rmse', 'time1', 'time2', 'geometry']], geometry='geometry')
        gdf.crs = proj4 #This is same for all tiles
        #gdf = gdf.to_crs('EPSG:4326') #OLD: ({'init': 'epsg:4326'}) # convert to change to lat/longitude #TODO either pass parameter to convert to lat/lon for icesat-2; or do conversion in the calling function
        # check why future warning started coming now
        #shp_json = gdf.to_json()
        #gdf.to_file('D:/wspace/rema_2m/{}.shp'.format(region))
        return gdf

#------------------------------------------------------------------------------------------------------------
def get_rema_strip_polygon(strip_folder, option=None):
    ''' Parses a REMA folder to reveal if one or more than one rema tiles are present in the folder
        strip_folder: Full path to 2m new REMA strips
        option: to match DEMs that do not match the default criteria of DEM ; eg hma folder structure was different
        Assumptions:
            strip folder name always end with _v0xxxx
            Eg: W1W1_20081215_1020010005DB7A00_1020010005DE1800_2m_lsf_v040000
    '''
    strip = strip_folder.split('/')[-1]
    # Separate strip name and strip version; because all the files inside this folder is named by strip only, no version number; we also want to strip version number to shapefile metadata file
    # strip = strip.split('_v0')[0] #OLD delete after a while [Aug 08, 2020]
    strip = strip.split('_')
    strip_version = strip[-1]
    strip = '_'.join(strip[:-1])

    files = os.listdir(strip_folder)
    metas = [f for f in files if f.endswith('_meta.txt')]
    if len(metas) == 0:
        # To Guard aginst empty/missing data inside the strip folder; but revise again; can cause other error downstream
        #print(f'Missing metadata: {strip_folder}')
        return
    for seg_id in range(1, len(metas)+1):
        #f'{dir_prefix}/REMA/{region}/strips_unf/2m/{strip}/{strip}_seg{seg_id}_meta.txt'
        meta_file = f'{strip_folder}/{strip}_seg{seg_id}_meta.txt' #This maybe be wrong; though overwritten used use to if/else condition below
        # Quick fix for HMA older files
        if option == 'hma':
            meta_file = f'{strip_folder}/{metas[seg_id-1]}'
        if seg_id == 1:
            gdf = extract_rema_polygon(meta_file)
            gdf['seg_id'] = seg_id
        else:
            gdf_temp = extract_rema_polygon(meta_file)
            gdf_temp['seg_id'] = seg_id
            gdf = pd.concat([gdf, gdf_temp])
    # Perhaps dissolve to simplify the geometry if more than one tiles present within a strip
    #if len(gdf)>1:
    #    gdf = gdf.dissolve(by='name')
    gdf['strip'] = strip # moved from another function to here because 
    gdf['version'] = strip_version # Ian wanted in newer version for metadata [Aut 05, 2020]
    gdf = gdf[['strip', 'version', 'seg_id', 'time1', 'time2', 'rmse', 'geometry']] #Reorganize for presentation
    return gdf

def fill_missing_time(row, region):
    ''' Fill in the time1 and time2 for older tiles that do not conform to new format 
        ALERT: This is hard-coded to REMA (see meta_file path below)
    '''
    # region = row['region']
    strip = row['strip']
    version = row['version']
    seg_id = row['seg_id']
    meta_file = f'{dem_folder}/{region}/strips_v4/2m/{strip}_{version}/{strip}_seg{seg_id}_meta.txt'
    dt1 = []
    dt2 = []
    with open(meta_file, 'r') as infile:
        for line in infile.readlines():
            if 'Image 1=' in line:
                file_parts = line.split('/')
                #region1 = '/'.join(file_parts[1:-3]) #file_parts[-3]
                img1 = file_parts[-1]
                img1 = img1.split('.tif')[0] #to remove extension (.tif) from name
                # From this img1 (ie part of filename, extract the datetime of image acquisition)
                dt = img1.split('_')[1] # sample '20091120125025' [YYYYMMDDhhmmss]
                dt = pd.to_datetime(dt)
                dt1.append(dt)
                #this can be readily converted datetime, eg pd.to_datetime('20091120125025') => Timestamp('2009-11-20 12:50:25')
            if 'Image 2=' in line:
                file_parts = line.split('/')
                img2 = file_parts[-1]
                img2 = img2.split('.tif')[0]
                dt = img2.split('_')[1]
                dt = pd.to_datetime(dt)
                dt2.append(dt)
    # Cacluate average acquisition times for images
    if len(dt1)>0:
        ser1 = pd.Series(dt1)
        time1 = ser1.mean()
        time1 = time1.strftime("%Y-%m-%d %H:%M:%S")
    if len(dt2)>0:
        ser2 = pd.Series(dt2)
        time2 = ser2.mean()
        time2 = time2.strftime("%Y-%m-%d %H:%M:%S")
    return time1, time2
#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
for region in regions:
    #First list all the dem strips
    strips_folder = f'{dem_folder}/{region}/strips_v4/2m' # Each strip inside this folder
    strips = os.listdir(f'{strips_folder}')
    print(f'Number of Strips for {region} = {len(strips)}')
    count = 0
    for strip in strips:
        #strip_folder = f'{dir_prefix}/{dem_type}/{region}/strips_unf/2m/{strip}'
        strip_folder = f'{dem_folder}/{region}/strips_v4/2m/{strip}'
        if count == 0:
            gdf = get_rema_strip_polygon(strip_folder) #, option='hma'
        else:
            gdf_temp = get_rema_strip_polygon(strip_folder) #, option='hma'
            #Concatenate
            gdf = pd.concat([gdf, gdf_temp])
        count += 1
    gdf.index = range(len(gdf)) #assingn index, else subsetting for missing wont work
    # Check if any row is missing time1 time2 due to strips being in old format
    # subset = gdf[gdf.time1.isnull()]
    #subset = gdf[gdf.time1.isnull()] #when saved shp file is open, empty string to read as null
    subset = gdf[gdf.time1==''] #Since the shp is not yet saved, it still contains empty string
    if len(subset)>0:
        '''This ad-hoc code is only for old files. Remove once only new files are processed PGC'''
        print(f'Missing time {len(subset)}')
        t = subset.apply(fill_missing_time, args=(region,), axis=1)
        gdf.loc[t.index, 'time1'] = t.apply(lambda x: x[0])
        gdf.loc[t.index, 'time2'] = t.apply(lambda x: x[1])
        # gdf.loc[t.index, 'flag'] = 'old'
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
# Extract Outlines of the shapefile (by dissolving) and concatenating all the outlines
# NB: Try merge instead of concat ()
count = 0
for shp_file in shp_files:
    gdf = gpd.read_file(f'{metadata_folder}/{shp_file}')
    # Dissolve Polygons within a shapefile to get only the Outline
    # Need to add some attribute based on which to dissolve
    gdf['region'] = shp_file.split('.')[0]#shp_file[:9]
    if count == 0:
        #No concatenation
        outline_gdf = gdf.dissolve(by='region')  ## This is a geodataframe with single row
        # Attach more metadata from the original shapefile
        #outline_gdf[['xmax', 'xmin', 'ymax', 'ymin']] = gdf.xmax.max(), gdf.xmin.min(), gdf.ymax.max(), gdf.ymin.min()
        #outline_gdf['Nscenes'] = gdf.Nscenes.sum()
    else:
        outline2 = gdf.dissolve(by='region')
        #outline2[['xmax', 'xmin', 'ymax', 'ymin']] = gdf.xmax.max(), gdf.xmin.min(), gdf.ymax.max(), gdf.ymin.min()
        #outline2['Nscenes'] = gdf.Nscenes.sum()
        #Concatenate
        outline_gdf = pd.concat([outline_gdf, outline2])
    count += 1
# Remove Columns that does not have meaning in merged product
outline_gdf = outline_gdf.drop(['version','seg_id','rmse', 'strip'], axis=1)
# Add more attributes after merge Area, Perimeter
outline_gdf['Area'] = np.round(outline_gdf.area/10**6) #round t0 1 or zero decimals
outline_gdf['Perim'] = np.round(outline_gdf.length/10**3)
outline_gdf = outline_gdf[['Area', 'Perim', 'geometry']]
outline_gdf.to_file(f'{metadata_folder}/{prefix}_regions.shp')
outline_gdf.to_file(f'{metadata_folder}/{prefix}_regions.geojson', driver='GeoJSON')

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# Part 3: Merge all strips into one Pan-strip-dem
count = 0
for shp_file in shp_files:
    if count == 0:
        gdf = gpd.read_file(f'{metadata_folder}/{shp_file}')
        gdf['region'] = shp_file.split('.')[0] # Attach region column
    else:
        gdf2 = gpd.read_file(f'{metadata_folder}/{shp_file}')
        gdf2['region'] = shp_file.split('.')[0] # Attach region column
        #Concatenate
        gdf = pd.concat([gdf, gdf2])
    count += 1
# Save shapefile
gdf.to_file(f'{metadata_folder}/{prefix}_strips.shp')
gdf.to_file(f'{metadata_folder}/{prefix}_strips.geojson', driver='GeoJSON')

# ------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print('Finished Script')
