# scripts for fetching and preparing external data (e.g. data cube)
import os
import warnings
import uuid
import datacube
import numpy as np
import geopandas as gpd
import xarray as xr
import rasterio
from PIL import Image
from base64 import b64encode
from io import BytesIO
from matplotlib import cm
from gdv_scripts import helpers
from dea_scripts import DEADataHandling  # old
from dea_scripts import dea_datahandling # new

# nci
#import sys
#sys.path.append('/home/573/lt5065/dea-notebooks/10_Scripts')  # todo remove
#import DEADataHandling                                        # todo remove


## HELPER FUNCTIONS ##
# drop original bands if they exist, fail silently if they dont
def drop_original_bands(ds):
    if ds:
        print('Dropping original blue, green, red, nir, swir1, swir2 bands.')
        ds = ds.drop(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'], errors='ignore')
        return ds
    else:
        print('No dataste provided during drop_original_bands. Aborting.')
        return None
    
# add empty data for first date for consistent resampling
def add_empty_first_date_array(ds, time_slice, first_query_date):
    
    def create_new_date_string(first_query_date, time_slice):
        yr = first_query_date.split('-')[0]
        mt = '{:02d}'.format(time_slice)
        dy = first_query_date.split('-')[2]
        return yr + '-' + mt + '-' + dy
    
    # tell user
    print('Inserting empty first date, if required, to improve resampling. Please wait.')
    
    # for specific month or season, correct month
    if time_slice == 'q1':
        first_query_date = np.datetime64(create_new_date_string(first_query_date, 1))
    elif time_slice == 'q2':
        first_query_date = np.datetime64(create_new_date_string(first_query_date, 4))
    elif time_slice == 'q3':
        first_query_date = np.datetime64(create_new_date_string(first_query_date, 7))
    elif time_slice == 'q4':
        first_query_date = np.datetime64(create_new_date_string(first_query_date, 10))
    elif isinstance(time_slice, int):
        first_query_date = np.datetime64(create_new_date_string(first_query_date, time_slice))
    else:
        first_query_date = np.datetime64(first_query_date)
        
    # check if already exists
    if ds['time'].isel(time=0) == first_query_date:
        return ds
        
    # generate empty and modify date
    empty_date = xr.full_like(ds.isel(time=0), fill_value=np.nan)
    empty_date['time'] = first_query_date
    
    # combine and sort
    ds = xr.concat([ds, empty_date], dim='time').sortby('time')
    
    return ds

# create reprojected tif image for interactive map
def create_and_reproject_image_from_array(da, out_folder, filename, filetype, mask_value):
    
    if da.any() and out_folder and filename and filetype:

        # create temp outfile locations and names
        filename = uuid.uuid1().hex + '_' + filename
        
        # prepare temp file locations
        temp_file_orig = os.path.join(out_folder, filename + filetype)
        temp_file_proj = os.path.join(out_folder, filename + '_' + 'proj' + filetype)
        
        # prepare output crs on da
        da.attrs['crs'] = 'EPSG:3577'
        datacube.helpers.write_geotiff(temp_file_orig, da.to_dataset())

        # set destination crs (wgs84)
        dst_crs = 'EPSG:4326' 

        # get original tif information
        with rasterio.open(temp_file_orig) as src:
            transform, width, height = rasterio.warp.calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update({'crs': dst_crs, 'transform': transform, 'width': width, 'height': height})

            # use rasterio to write out the new projected raster
            with rasterio.open(temp_file_proj, 'w', **kwargs) as dst:
                rasterio.warp.reproject(source=rasterio.band(src, 1), destination=rasterio.band(dst, 1), 
                                        src_transform=src.transform, src_crs=src.crs, dst_transform=transform, 
                                        dst_crs=dst_crs, resampling=rasterio.warp.Resampling.nearest)

        # use rasterio to import the reprojected data as img
        with rasterio.open(temp_file_proj) as src:
            boundary = src.bounds
            img = src.read()
            nodata = src.nodata

        # clena up img for zeros, cmaps
        if mask_value is not None:
            img = np.ma.masked_where(img[0] == mask_value, img[0])
        else:
            img = img[0]

        try:
            # delete temp files
            os.remove(temp_file_orig)
            os.remove(temp_file_proj)
        except:
            print('Problem during temp file deletion, skipping deletion.')
            return None, None

        return img, boundary

    else:
        print('Did not provide a data array, index, or file output location during image generation.')

# create png for ipyleaflet mapping
def generate_png_url(in_img, color_map):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
   
        try:
            # normalise image color range
            img_norm = in_img - np.nanpercentile(in_img, q=0.5)
            img_norm = img_norm / np.nanmax(img_norm)

        except:
            print('Problem with image normalisation.')
            return None

    try:
        # get colormap from user
        cmap = cm.get_cmap(color_map)

        # create pil image in memory
        pil_img = Image.fromarray(np.uint8(cmap(img_norm) * 255))
        mask = np.where(np.isfinite(in_img), 255, 0)
        pil_mask = Image.fromarray(np.uint8(mask), mode='L')
        
    except:
        print('Problem with pil image creation and masking.')
        raise
        return None
    try:
        # get image size
        size = in_img.shape[::-1]

        # create new img in memory
        out_img = Image.new(mode='RGBA', size=size, color=None)
        out_img.paste(pil_img, mask=pil_mask)        
        
        # resize to remove blur
        size = (size[0] * 3, size[1] * 3)
    except:
        print('Problem with pill image copy and sizing.')
        return None

    # open in mem
    with BytesIO() as out:
        out_img = out_img.resize(size, resample=Image.NEAREST)
        out_img.save(out, 'png')
        contents = out.getvalue()
        contents = b64encode(contents)
        contents = contents.decode('ascii')
        img_url = 'data:image/png;base64,' + contents
    
    return img_url
        

## REDUCER FUNCTIONS ##
# reduce a given dataset down into wet (djf) and dry (son) months
def reduce_ds_to_wet_and_dry(ds):
    
    def is_wet_or_dry(month):
        return ((month >= 1) & (month <= 3)) | ((month >= 9) & (month <= 11))
    
    if not ds:
        print('No dataset provided during reduction to wet and dry season. Aborting.')
        return None  
    
    # tell user and reduce dataset to months in wet (JFM) or dry (SON) seasons only
    print('Reducing dataset down into wet (JFM) and dry (SON) seasons only.')
    ds = ds.sel(time=is_wet_or_dry(ds['time.month']))
         
    return ds

# reduce to specific quarter/month (original mannkendall or change analysis only)
def reduce_to_specific_month_or_season(ds, time_slice):
    try:
        if ds and time_slice:              
            if time_slice == 'q1':
                print('> Dataset being reduced to Quarter 1 (DJF) months. Proceeding.')
                return ds.sel(time=((ds['time.month'] == 12) | (ds['time.month'] == 1)) | ((ds['time.month'] == 2)))
            elif time_slice == 'q2':
                print('> Dataset being reduced to Quarter 2 (MAM) months. Proceeding.')
                return ds.sel(time=((ds['time.month'] == 3) | (ds['time.month'] == 4)) | ((ds['time.month'] == 5)))
            elif time_slice == 'q3':
                print('> Dataset being reduced to Quarter 3 (JJA) months. Proceeding.')
                return ds.sel(time=((ds['time.month'] == 6) | (ds['time.month'] == 7)) | ((ds['time.month'] == 8)))
            elif time_slice == 'q4':
                print('> Dataset being reduced to Quarter 4 (SON) months. Proceeding.')
                return ds.sel(time=((ds['time.month'] == 9) | (ds['time.month'] == 10)) | ((ds['time.month'] == 11)))
            elif isinstance(time_slice, int):
                print('> Dataset being reduced to month: {0}. Proceeding.'.format(time_slice))
                return ds.sel(time=((ds['time.month'] == time_slice)))
            else:
                print('> No adequate trend time slice. No reduction done.')
                return None
        else:
            print('> No dataset or trend dimension provided. Please check.')
            return None 
    except:
        print('Something went wrong. Check your trend dimension.')
        return None

    
## MAIN FUNCTIONS ##
# fetch data cube data for landsat and sentinel
def fetch_cube_data(platform, product, q_dates, cloud_cover, cloud_mask, coords_nw, 
                    coords_se, in_epsg='EPSG:3577', out_epsg='EPSG:3577'):
    
    # general clean up function (drop data perc, normalize time, sort
    def general_clean_up(ds):
        
        # tell user
        print('Doing some remaining dataset clean up (e.g. dropping mask var, normalising time, sorting)')
    
        # if data_perc exists, drop
        ds = ds.drop('data_perc', errors='ignore')
            
        # normalise datetime to remove seconds, minutes
        ds['time'] = ds.indexes['time'].normalize()

        # sort by date in case anything moved
        ds = ds.sortby('time')
        
        return ds
    
    # tell user
    print('Beginning satellite data fetch. This can take awhile.\n')
    
    # create message regarding params
    print('=== Your Chosen Parameters ===')
    print('Platform: {0}.'.format(platform))
    print('Product: {0}.'.format(product))
    print('Date Range: {0} to {1}.'.format(q_dates[0], q_dates[1]))
    print('Cloud Cover: {0}.'.format(str(100 - (cloud_cover * 100))))
    print('Mask Clouds: {0}.'.format(str(cloud_mask)))
    print('NW Corner: {0} | SE Corner {1}.'.format(str(coords_nw), str(coords_se)))
    print('Input EPSG: {0} | Output EPSG {1}.'.format(in_epsg, out_epsg))
    print('=== Your Chosen Parameters ===\n')
    
    try:
        # build query
        x = (coords_nw[0], coords_se[0])
        y = (coords_nw[1], coords_se[1])
        query = {'x': x, 'y': y, 'time': q_dates, 'crs': in_epsg, 'output_crs': out_epsg}

        # get datacube object
        dc = datacube.Datacube(app='gdv_tool')
        
        # check platform, get bands, proceed
        if platform == 'ls_old':

            # set bands, we reduce later
            bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']

            # add resolution key val pair appropriate for ls
            query.update(resolution = (-25, 25))

            # pull from dea load_clearlandsat script TODO OLD METHOD CHANGE TO dea_datahandling
            ds = DEADataHandling.load_clearlandsat(dc=dc, query=query, bands_of_interest=bands, 
                                                   sensors=['ls5', 'ls7', 'ls8'], product=product, 
                                                   masked_prop=cloud_cover, ls7_slc_off=False, 
                                                   mask_invalid_data=True, mask_pixel_quality=cloud_mask,
                                                   satellite_metadata=False, dask_chunks={'time': 1}, 
                                                   lazy_load=True)

            if ds:
                # tell user and return
                print('\nSuccessfully fetched Landsat satellite data!')

                # clean up
                ds = general_clean_up(ds=ds)

                # return
                return ds
            else:
                print('Aborting data fetch due to no valid data for query. Please change date range or platform.')
                return None

        elif platform == 's2_old':

            # set bands depending on product selection, we reduce later
            bands = [product + '_blue', product + '_green', product + '_red', product + '_nir_1', 
                     product + '_swir_2', product + '_swir_3']

            # add resolution key val pair appropriate for sent2
            query.update(resolution = (-10, 10))

            # pull from dea load_clearsentinel script modified for lazy load
            ds = DEADataHandling.load_clearsentinel2_mod(dc=dc, query=query, sensors=['s2a', 's2b'],
                                                         bands_of_interest=bands, masked_prop=cloud_cover,
                                                         mask_invalid_data=True, mask_pixel_quality=cloud_mask, 
                                                         satellite_metadata=False, dask_chunks={'time': 1}, 
                                                         lazy_load=True)

            if ds:
                # rename bands for output
                print('Renaming Sentinel 2 bands to match landsat standards.')
                for var in ds.data_vars:
                    if var == product + '_blue':
                        ds = ds.rename({var: 'blue'})
                    elif var == product + '_green':
                        ds = ds.rename({var: 'green'})
                    elif var == product + '_red':
                        ds = ds.rename({var: 'red'})
                    elif var == product + '_nir_1':
                        ds = ds.rename({var: 'nir'})                
                    elif var == product + '_swir_2':
                        ds = ds.rename({var: 'swir1'})                  
                    elif var == product + '_swir_3':
                        ds = ds.rename({var: 'swir2'}) 

                # tell user and return
                print('\nSuccessfully fetched Sentinel 2 satellite data!')

                # clean up
                ds = general_clean_up(ds=ds)

                # return
                return ds
            else:
                print('Aborting data fetch due to no valid data for query. Please change date range or platform.')
                return None

        elif platform == 'ls':

            # set bands depending on product selection, we reduce later
            bands = [product + '_blue', product + '_green', product + '_red', product + '_nir', 
                     product + '_swir_1', product + '_swir_2']

            # add resolution key val pair appropriate for ls
            query.update(measurements = bands)
            query.update(resolution = (-30, 30))

            # move this over once on sandbox
            ds = dea_datahandling.load_ard(dc=dc, products=['ga_ls5t_ard_3', 'ga_ls7e_ard_3', 'ga_ls8c_ard_3'],
                                           min_gooddata=cloud_cover, mask_pixel_quality=cloud_mask, 
                                           mask_invalid_data=True, ls7_slc_off=False, 
                                           dask_chunks={'time': 1}, lazy_load=True, **query)

            if ds:
                # tell user and return
                print('\nSuccessfully fetched Landsat satellite data!')
                
                # rename bands for output
                print('Renaming Landsat bands for clarity.')
                for var in ds.data_vars:
                    if var == product + '_blue':
                        ds = ds.rename({var: 'blue'})
                    elif var == product + '_green':
                        ds = ds.rename({var: 'green'})
                    elif var == product + '_red':
                        ds = ds.rename({var: 'red'})
                    elif var == product + '_nir':
                        ds = ds.rename({var: 'nir'})                
                    elif var == product + '_swir_1':
                        ds = ds.rename({var: 'swir1'})                  
                    elif var == product + '_swir_2':
                        ds = ds.rename({var: 'swir2'}) 

                # clean up
                ds = general_clean_up(ds=ds)

                # return
                return ds
            else:
                print('Aborting data fetch due to no valid data for query. Please change date range or platform.')
                return None
        
        elif platform == 's2':

            # set bands depending on product selection, we reduce later
            bands = [product + '_blue', product + '_green', product + '_red', product + '_nir_1', 
                     product + '_swir_2', product + '_swir_3']

            # add resolution key val pair appropriate for ls
            query.update(measurements = bands)
            query.update(resolution = (-10, 10))

            # move this over once on sandbox
            ds = dea_datahandling.load_ard(dc=dc, products=['s2a_ard_granule', 's2b_ard_granule'],
                                           min_gooddata=cloud_cover, mask_pixel_quality=cloud_mask, 
                                           mask_invalid_data=True, ls7_slc_off=False, 
                                           dask_chunks={'time': 1}, lazy_load=True, **query)

            if ds:
                # tell user and return
                print('\nSuccessfully fetched Sentinel satellite data!')
                
                # rename bands for output
                print('Renaming Sentinel bands for clarity.')
                for var in ds.data_vars:
                    if var == product + '_blue':
                        ds = ds.rename({var: 'blue'})
                    elif var == product + '_green':
                        ds = ds.rename({var: 'green'})
                    elif var == product + '_red':
                        ds = ds.rename({var: 'red'})
                    elif var == product + '_nir_1':
                        ds = ds.rename({var: 'nir'})                
                    elif var == product + '_swir_2':
                        ds = ds.rename({var: 'swir1'})                  
                    elif var == product + '_swir_3':
                        ds = ds.rename({var: 'swir2'}) 

                # clean up
                ds = general_clean_up(ds=ds)

                # return
                return ds
            else:
                print('Aborting data fetch due to no valid data for query. Please change date range or platform.')
                return None
        
    except Exception as e:
            print('Error occurred during fetch_cube_data of type {0}. Stopping.'.format(e))
            raise e    
            
# fetch roc shapefile and prepare into actual, prediction lists
def open_and_prepare_shapefile_for_roc(da, shapefile_path):

    if not shapefile_path.endswith('.shp'):
        print('Did not select a shapefile. Make sure the file you selected is of type .shp. Aborting.')
        return None
    
    # load shapefile
    shp = gpd.read_file(shapefile_path)
        
    # checks
    if shp.crs['init'] != 'epsg:3577':
        print('> Incorrect CRS detected. Needs to be Australia Albers, EPSG:3577. Please reproject in a GIS and try again.')
        return
    elif shp.geom_type.all() != 'Point':
        print('> Incorrect shapefile type detected (Point). Please ensure shapefile is of type Point and try again.')
        return
    elif 'gdv_act' not in [col.lower() for col in shp.columns]:
        print('> Incorrect groundwater dependent species location column (GDV_ACT) found in data. Please ensure observations are in column called GDV_ACT.')
        return
    elif len(shp['GDV_ACT'].unique().tolist()) != 2:
        print('> More than 0 and 1 in GDV_ACT column detected. Please ensure only 0 (absent) or 1 (present) exist.')
        return
    elif not all(val in shp['GDV_ACT'].unique().tolist() for val in [0, 1]):
        print('> Incorrect values in GDV_ACT column detected. Please ensure only 0 (absent) or 1 (present) exist.')
        return
    
    # tell user
    print('Groundtruth shapefile opened successfully. Beginning to extract data for ROC analysis. Please hold.')
    
    # add new column for predictions
    shp['GDV_PRED'] = None

    # intersect values from likelihood to points
    for idx, row in shp.iterrows():
                
        x = row['geometry'].x
        y = row['geometry'].y

        # get value from xarray at x and y loc
        try:
            val = da.sel(x=x, y=y, method='nearest', tolerance=30)
            val = val.values
            val = val.item()
        except:
            val = None

        # update row value at index
        if val:
            shp.loc[idx, 'GDV_PRED'] = val
            
    # drop where prediction is na
    shp = shp.dropna(subset=['GDV_PRED'])

    if len(shp) > 0:
        print('Sucessfully prepared shapefile for ROC analysis. Num of ground sites within selected study area: {0}'.format(len(shp)))
        return shp

    else:
        print('No groundtruthed points intersected with selected study area. Aborting.')
        return None

# export likelihood as netcdf
def write_like_data_path(da_raw, da_thresh, out_path):
    try:
        msg = 'Exporting groundwater dependent vegetation likelihood maps (raw and thresholded). '
        msg += 'Files being written as netcdf files for ease of input/output.\nLocation of export: {0}'      
        print(msg.format(out_path))
        
        # export
        da_raw.to_netcdf(out_path + '_raw.nc', mode='w', engine='scipy')
        da_thresh.to_netcdf(out_path + '_thresh.nc', mode='w', engine='scipy')
        
    except:
        print('Error during write_like_data_path. Aborting.')

# import likelihood from netcdf
def read_like_data_path(in_path):
    try:
        msg = 'Importing groundwater dependent vegetation likelihood data (should be thresholded). '
        msg += 'Files being read as netcdf file.\nLocation of import: {0}'      
        print(msg.format(in_path))
        
        # check if nc type
        if in_path.endswith('.nc'):
            ds = xr.open_dataset(in_path).compute()
            print('File successully read. Proceeding.\n')
            return ds
        else:
            print('Selected likelihood file is not a netcdf file type (.nc). Please select an .nc file.')
            return None
    except:
        print('Error during write_like_data_path. Aborting.')
        return None

# load imported netcdf
def load_and_prepare_like_netcdf(params):

    # load data, if not valid, abort
    if params.read_data_path:
        
        # tell user
        print('Loading external netcdf file. Please wait.')
        
        # read data and set crs
        scene_like_threshed = read_like_data_path(params.read_data_path)
        scene_like_threshed.attrs['crs'] = 'EPSG:3577'
        
        # get coords nw and se
        n, w = float(scene_like_threshed['y'].max()), float(scene_like_threshed['x'].min())
        s, e = float(scene_like_threshed['y'].min()), float(scene_like_threshed['x'].max())
        
        # pop crs, affine and albers in param
        params.affine, params.crs = scene_like_threshed.affine, scene_like_threshed.crs
        params.alb_nw, params.alb_se = (w, n),  (e, s)
        
        # project to wgs84 and add to params
        params.gcs_nw, params.gcs_se = (helpers.transform_point(e, s, 'epsg:3577', 'epsg:4326')), (helpers.transform_point(w, n, 'epsg:3577', 'epsg:4326'))

        # tell user and return
        print('Loaded external netcdf file Sucessfully. Proceeding.')
        return scene_like_threshed, params
    else:
        print('No file path selected in above control. Please select a netcdf file (.nc) and try again.')
        return None, params