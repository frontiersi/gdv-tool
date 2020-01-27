# scripts for calculating various vegetation and moisture indices
import warnings
import xarray as xr
import pandas as pd
import numpy as np

from gdv_scripts import data
from dea_scripts.BandIndices import tasseled_cap

#import sys                                                    # todo remove
#sys.path.append('/home/573/lt5065/dea-notebooks/10_Scripts')  # todo remove
#from BandIndices import tasseled_cap

## HELPER FUNCTIONS ##
# check if blue, green, red, nir, swir1, swir2 bands exist, else False is returned
def check_required_bands_exist(ds):
    
    # basic checks
    if not ds:
        print('No dataset provided. Please ensure all previous cells have been run.')
        return False
    elif not len(ds.data_vars) > 0:
        print('No bands in dataset. Try fetching your satellite data again.')
        return False
    
    try:
        existing_bands = [b for b in ds.data_vars]
        required_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
        result = all(b in existing_bands for b in required_bands)
        
        return result
   
    except Exception as e:
        print('Error occurred during check_required_bands_exist of type {0}. Stopping.'.format(e))
        raise e  

# used in apply function, loop years, count months in each, populate given list with missing wet/dry dates
def get_missing_dates_wet_dry(da, miss_seasons):
    if len(da['time.month']) > 0 and len(da['time.month']) < 2:
        if da['time.month'] == 1:
            missing_month = 9
        elif da['time.month'] == 9:
            missing_month = 1
        else:
            print('Error occured. Months are not 1 or 9. Aborting.')
            return da

        # get scene, date, and timestamp
        scene = da.isel(time=0)
        date = pd.to_datetime(scene['time'].values)
        ts = pd.Timestamp(date.year, missing_month, date.day)
        date = np.datetime64(ts)

        miss_seasons.append(date)
        
    return da


## MAIN FUNCTIONS ##
# calculate one of the nine vegetation indices at user request, rescale if desired
def calc_veg_index(ds, index, rescale, l_value):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please fetch your datacube again.')
        return None
    elif not index:
        print('No vegetation index set. Please ensure they are set above above.')
        return None
    
    # tell user
    print('Calculating vegetation index: {0}. Please wait.'.format(index))
    
    try:
        # ndvi
        if index == 'ndvi':
            ds['veg_idx'] = (ds.nir - ds.red) / (ds.nir + ds.red)
            if rescale:
                ds['veg_idx'] = ds['veg_idx'] + 1

        # savi
        elif index == 'savi':
            ds['veg_idx'] = ((ds.nir - ds.red) / (ds.nir + ds.red + l)) * (1 + l_value)
            if rescale:
                ds['veg_idx'] = ds['veg_idx'] + 1

        # stvi3
        elif index == 'stvi3':
            ds['veg_idx'] = ds.nir / (ds.red + ds.swir1)

        # slavi
        elif index == 'slavi':
            ds['veg_idx'] = ds.nir / (ds.red + ds.swir2)

        # msvi1
        elif index == 'msvi1':
            ds['veg_idx'] = ds.nir / ds.swir1

        # msvi2
        elif index == 'msvi2':
            ds['veg_idx'] = ds.nir / ds.swir2

        # msvi3
        elif index == 'msvi3':
            ds['veg_idx'] = ds.nir / (ds.swir1 + ds.swir2)

        # mavi
        elif index == 'mavi':
            ds['veg_idx'] = (ds.nir - ds.red) / (ds.nir + ds.red + ds.swir1)
            if rescale:
                ds['veg_idx'] = ds['veg_idx'] + 1

        # tasseled cap greenness (via dea tool), rename bands
        elif index == 'tcap':
            ds = tasseled_cap(ds[['blue', 'green', 'red', 'nir', 'swir1', 'swir2']], 
                              tc_bands=['greenness', 'brightness'], drop=False)
            ds = ds.rename({'greenness': 'veg_idx', 'brightness': 'brt_idx'})
            
            if rescale:
                ds['veg_idx'] = ds['veg_idx'] / 10000
                ds['brt_idx'] = ds['brt_idx'] / 10000

        # tell user
        print('Successfully generated vegetation index: {0}. Proceeding.'.format(index))
        
        return ds
    
    except Exception as e:
        print('Error occurred during calculate_veg_index of type {0}. Stopping.'.format(e))
        raise e  

# calculate one of the two moisture indices at user request, rescale if desired
def calc_mst_index(ds, index, rescale):
     
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please fetch your datacube again.')
        return None
    elif not index:
        print('No moisture index set. Please ensure they are set above above.')
        return None
    
    # tell user
    print('Calculating moisture index: {0}. Please wait.'.format(index))
    
    try:
        # ndmi
        if index == 'ndmi':
            ds['mst_idx'] = (ds.nir - ds.swir1) / (ds.nir + ds.swir1)
            if rescale:
                ds['mst_idx'] = ds['mst_idx'] + 1

        # gvmi    
        elif index == 'gvmi':
            ds['mst_idx'] = ((ds.nir + 0.1) - (ds.swir2 + 0.02)) / ((ds.nir + 0.1) + (ds.swir2 + 0.02))
            if rescale:
                ds['mst_idx'] = ds['mst_idx'] + 1

        # tell user
        print('Successfully generated moisture index: {0}. Proceeding.'.format(index))
        
        return ds
    
    except Exception as e:
        print('Error occurred during calculate_mst_index of type {0}. Stopping.'.format(e))
        raise e  

# calculate the one of one fire indices at user request
def calc_nbr_index(ds, index):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please fetch your datacube again.')
        return None
    elif not index:
        print('No bush fire index set. Please ensure they are set above above.')
        return None
    
    # tell user
    print('Calculating fire index: {0}. Please wait.'.format(index))
    
    try:
        # nbr
        if index == 'nbr':
            ds['nbr_idx'] = (ds.nir - ds.swir2) / (ds.nir + ds.swir2)

        # tell user
        print('Successfully generated fire index: {0}. Proceeding.'.format(index))
        
        return ds
    
    except Exception as e:
        print('Error occurred during calc_nbr_index of type {0}. Stopping.'.format(e))
        raise e  

# calculate all indices, drop bands, rename, resample for dry and wet seasons
def create_resample_idxs_wet_and_dry(ds, params):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please re-run the tool again.')
        return None
    elif not params.check_raw_idx_params():
        print('No vegetation or moisture index set. Please ensure they are set above above.')
        return None
    
    # generate vegetation and moisture indices based on selection, add as vars to existing ds
    ds = calc_veg_index(ds=ds, index=params.veg_idx, rescale=params.rescale, l_value=params.l_value)
    ds = calc_mst_index(ds=ds, index=params.mst_idx, rescale=params.rescale)
    
    # drop original bands
    ds = data.drop_original_bands(ds)
    
    # rename and drop tcap brightness band if tcap requested - dont need yet
    if params.veg_idx == 'tcap':
        ds = ds.drop(['brt_idx'], errors='ignore')
            
    # basic check
    if not ds:
        print('No dataset remaining after drop original bands. Aborting.')
        return None   
    
    # resample into two seasons. always going to be 01-01 and 09-01.
    print('\nResampling data down to annual wet and dry seasonal medians.')
    ds_wet = ds.where(ds['time.month'] <= 3, drop=True).resample(time='AS-JAN')
    ds_dry = ds.where(ds['time.month'] >= 9, drop=True).resample(time='AS-SEP')
       
    # compute dataset via reduction (.median() doesnt work on dask yet...). disable nan warnings (i want them!)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            print('Computing dask during reduction. This can take quite awhile. Please wait.')
            ds = ds_wet.reduce(np.nanmedian, dim='time', keep_attrs=True)
            ds = xr.concat([ds, ds_dry.reduce(np.nanmedian, dim='time', keep_attrs=True)], dim='time').sortby('time')
            ds = ds.compute() ### ADDED
            
            print('Computated dataset successfully! Proceeding.\n')
    except Exception as e:
        print('Error occurred during compute of type {0}. Stopping.'.format(e))
        print('Please reduce the number of years and/or study area size.')
        raise e
        
    # get dates of missing wet/dry dates, append on as all nan dates
    miss_seasons = []
    ds.groupby('time.year').apply(get_missing_dates_wet_dry, miss_seasons=miss_seasons)
    if len(miss_seasons) > 0:
        for date in np.array(miss_seasons, dtype='datetime64[ns]'):
            new_scene = xr.full_like(ds.isel(time=0), np.nan)
            new_scene.coords['time'] = date
            ds = xr.concat([ds, new_scene], dim='time').sortby('time')

    # basic check
    if not ds:
        print('No dataset remaining after reduction. Aborting.')
        return None
    else:
        return ds

# calculate nbr index
def create_resample_nbr_annual(ds):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please re-run the tool again.')
        return None
    
    # drop original bands
    ds = data.drop_original_bands(ds)
            
    # basic check
    if not ds:
        print('No dataset remaining after drop original bands. Aborting.')
        return None   
    
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)

            # resample into two seasons. always going to be 01-01 and 09-01.
            print('\nResampling data down to annual burn areas. Please wait.')
            ds = ds.resample(time='AS').mean()

            # reverse order and difference each date from next to subtract post from pre. sort again back to normal afterwards
            ds = ds.sortby('time', ascending=False).diff('time', label='lower')
            ds = ds.sortby('time')

            # compute
            print('Computing dask during reduction. This can take quite awhile. Please wait.')
            ds = ds.compute()
    
    except Exception as e:
        print('Error occurred during compute of type {0}. Stopping.'.format(e))
        print('Please reduce the number of years and/or study area size.')
        raise e    

    # basic check
    if not ds:
        print('No dataset remaining after reduction. Aborting.')
        return None
    else:
        return ds
    
# calculate all indices, drop bands, rename, resample for trend time slice
def create_resample_idxs_trend_time_slice(ds, params):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please re-run the tool again.')
        return None
    elif not params.check_trend_idx_params():
        print('No vegetation index set. Please ensure they are set above above.')
        return None
    
    # generate vegetation and moisture indices based on selection, add as vars to existing ds
    ds = calc_veg_index(ds=ds, index=params.veg_idx, rescale=params.rescale, l_value=params.l_value)
    
    # drop original bands
    ds = data.drop_original_bands(ds)
    
    # rename and drop tcap brightness band if tcap requested - dont need yet
    if params.veg_idx == 'tcap':
        ds = ds.drop(['brt_idx'], errors='ignore')
            
    # basic check
    if not ds:
        print('No dataset remaining after drop original bands. Aborting.')
        return None   
    
    # resample into months/seasons
    print('Resampling data down to quarterly/monthly medians. Please hold, this can take awhile.')

    # compute dataset via reduction (.median() doesnt work on dask yet...). disable nan warnings (i want them!)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            
            if params.time_slice in ['qt', 'q1', 'q2', 'q3', 'q4']:
                ds = ds.resample(time='QS-DEC', label='right')
                ds = ds.reduce(np.nanmedian, dim='time', keep_attrs=True)
                ds = ds.compute() ### ADDED
                                
                if params.time_slice == 'qt':
                    print('No need to reduce further - all quarters required (DJF, MAM, JJA, SON). Proceeding.')
                elif params.time_slice == 'q1':
                    print('Reducing data further to Quarter 1 (DJF), please wait...')
                    ds = ds.where(ds['time.month'] == 3, drop=True)
                elif params.time_slice == 'q2':
                    print('Reducing data further to Quarter 2 (MAM), please wait...')
                    ds = ds.where(ds['time.month'] == 6, drop=True)
                elif params.time_slice == 'q3':
                    print('Reducing data further to Quarter 3 (JJA), please wait...')
                    ds = ds.where(ds['time.month'] == 9, drop=True)
                elif params.time_slice == 'q4':
                    print('Reducing data further to Quarter 4 (SON), please wait...')
                    ds = ds.where(ds['time.month'] == 12, drop=True)

            elif params.time_slice == 'mt' or isinstance(params.time_slice, int):
                ds = ds.resample(time='1MS')
                ds = ds.reduce(np.nanmedian, dim='time', keep_attrs=True)
                ds = ds.compute() ### ADDED
                
                if params.time_slice == 'mt':
                    print('No need to reduce further - all months required (JAN-DEC). Proceeding.')
                elif isinstance(params.time_slice, int):
                    print('Reducing data further to specific month: {0}, pleade wait...'.format(params.time_slice))
                    ds = ds.where(ds['time.month'] == params.time_slice, drop=True)

    except Exception as e:
        print('Error occurred during compute of type {0}. Stopping.'.format(e))
        print('Please reduce the number of years and/or study area size.')
        raise e

    # basic check
    if not ds:
        print('No dataset remaining after reduction. Aborting.')
        return None
    else:
        print('Computated dataset successfully! Proceeding.\n')
        return ds
    
# calculate all indices, drop bands, rename, resample for change time slice
def create_resample_idxs_change_time_slice(ds, params):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif not check_required_bands_exist(ds):
        print('Not all required bands available in dataset. Please re-run the tool again.')
        return None
    elif not params.check_change_idx_params():
        print('No vegetation index set. Please ensure they are set above above.')
        return None
    
    # generate vegetation and moisture indices based on selection, add as vars to existing ds
    ds = calc_veg_index(ds=ds, index=params.veg_idx, rescale=params.rescale, l_value=0.5)
    
    # drop original bands
    ds = data.drop_original_bands(ds)
                
    # basic check
    if not ds:
        print('No dataset remaining after drop original bands. Aborting.')
        return None   
    
    # resample into months/seasons
    print('Resampling data down to quarterly/monthly medians. Please hold, this can take awhile.')

    # compute dataset via reduction (.median() doesnt work on dask yet...). disable nan warnings (i want them!)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            
            if params.time_slice in ['q1', 'q2', 'q3', 'q4']:
                ds = ds.resample(time='QS-DEC', label='right')
                ds = ds.reduce(np.nanmedian, dim='time', keep_attrs=True)
                ds = ds.compute() ### ADDED
                if params.time_slice == 'q1':
                    print('Reducing data further to Quarter 1 (DJF), please wait...')
                    ds = ds.where(ds['time.month'] == 3, drop=True)
                elif params.time_slice == 'q2':
                    print('Reducing data further to Quarter 2 (MAM), please wait...')
                    ds = ds.where(ds['time.month'] == 6, drop=True)
                elif params.time_slice == 'q3':
                    print('Reducing data further to Quarter 3 (JJA), please wait...')
                    ds = ds.where(ds['time.month'] == 9, drop=True)
                elif params.time_slice == 'q4':
                    print('Reducing data further to Quarter 4 (SON), please wait...')
                    ds = ds.where(ds['time.month'] == 12, drop=True)

            elif isinstance(params.time_slice, int):
                ds = ds.resample(time='1MS')
                ds = ds.reduce(np.nanmedian, dim='time', keep_attrs=True)
                ds = ds.compute() ### ADDED
                if isinstance(params.time_slice, int):
                    print('Reducing data further to specific month: {0}, pleade wait...'.format(params.time_slice))
                    ds = ds.where(ds['time.month'] == params.time_slice, drop=True)

    except Exception as e:
        print('Error occurred during compute of type {0}. Stopping.'.format(e))
        print('Please reduce the number of years and/or study area size.')
        raise e

    # basic check
    if not ds:
        print('No dataset remaining after reduction. Aborting.')
        return None
    else:
        print('Computated dataset successfully! Proceeding.\n')
        return ds