# scripts for plotting various interpolation methods
import warnings
import numpy as np
import pandas as pd
import xarray as xr

## HELPER METHODS ##
# drop final year if any months within are missing or all nan, else skip
def drop_incomplete_first_year_wet_dry(ds):
    try:
        # tell user
        print('Checking if first year has missing seasons, and dropping if so. Please wait.')

        # get last year and months
        year_first, ds_first = list(ds.groupby('time.year'))[0]
        if year_first and len(ds_first['time.month']) != 2:
            ds = ds.where(ds['time.year'] != year_first, drop=True)
            print('Not enough data in first year: {0}. Had to drop the whole year.'.format(year_first))
        else:
            print('Sufficient data in first year: {0}. No data removed. Proceeding.'.format(year_first))

        return ds

    except Exception as e:
        print('Error occurred during drop_incomplete_first_year_wet_dry of type {0}. Stopping.'.format(e))
        raise e

# drop final year if any months within are missing or all nan, else skip
def drop_incomplete_final_year_wet_dry(ds):
    try:
        # tell user
        print('Checking if final year has missing seasons, and dropping if so. Please wait.')

        # get last year and months, drop if not == 2 after drop
        year_last, ds_last = list(ds.groupby('time.year'))[-1]
        if year_last and len(ds_last.dropna(dim='time', how='all')['time']) != 2:
            ds = ds.where(ds['time.year'] != year_last, drop=True)
            print('Not enough data in final year: {0}. Had to drop the whole year.'.format(year_last))
        else:
            print('Sufficient data in final year: {0}. No data removed. Proceeding.'.format(year_last))

        return ds

    except Exception as e:
        print('Error occurred during drop_incomplete_final_year_wet_dry of type {0}. Stopping.'.format(e))
        raise e
        
# flag whole years if any all nan months exist
def flag_years_with_missing_wet_or_dry(ds):
    print('\nChecking for years with wet/dry season data gaps, will flag and remove if found. Please wait.')
    try:
        # set full year to nan if any one season nan in it
        nan_dates = np.setdiff1d(ds['time'], ds.dropna(dim='time', how='all')['time'])
        nan_years = pd.to_datetime(nan_dates).year.drop_duplicates()

        if len(nan_years) > 0:
            print('Following years contain wet/dry season gaps: {0}'.format(list(nan_years)))

            ds = ds.where(ds['time.year'].isin(nan_years) == False)
            ds = ds.drop(['year_outlier'], errors='ignore')
        else:
            print('No years contain wet/dry season data gaps. Proceeding.')
                
        return ds
    
    except Exception as e:
        print('Error occurred during flag_years_with_missing_wet_or_dry of type {0}. Stopping.'.format(e))
        raise e


## MAIN FUNCTIONS ##
# take ds, interpolate linear along wet then dry seasons seperately
def interp_wet_dry(ds):

    # get missing dates for wet, dry seperate
    print('Beginning linear interpolation of missing data. Please wait.')
    nan_dates_wet = np.setdiff1d(ds.where(ds['time.month'] < 6, drop=True)['time'], 
                                 ds.where(ds['time.month'] < 6, drop=True).dropna(dim='time', how='all')['time'])
    nan_dates_dry = np.setdiff1d(ds.where(ds['time.month'] > 6, drop=True)['time'], 
                                 ds.where(ds['time.month'] > 6, drop=True).dropna(dim='time', how='all')['time'])

    # drop any scenes where all nan
    ds = ds.dropna('time', how='all')

    if len(nan_dates_wet) > 0:
        # interp along time for wet, dry seperately and extrapolate first year if all nan (not last year)
        interp_wet = ds.where(ds['time.month'] < 6, drop=True).interp(time=nan_dates_wet, kwargs={'fill_value': 'extrapolate'})
        ds = xr.concat([ds, interp_wet], dim='time').sortby('time')
    else:
        print('No need to interpolate wet season data. Skipping and proceeding.')
                
    if len(nan_dates_dry) > 0:
        interp_dry = ds.where(ds['time.month'] > 6, drop=True).interp(time=nan_dates_dry, kwargs={'fill_value': 'extrapolate'})
        ds = xr.concat([ds, interp_dry], dim='time').sortby('time')
    else:
        print('No need to interpolate dry season data. Skipping and proceeding.')

    # tell user and return
    print('Successfully interpolated missing data. Proceeding.')
    return ds

# take ds, interpolate (fill) missing pixels within each singular image
def fill_pixels_wet_dry(ds):
    try:
        print('\nFilling in missing pixels. Warning, this can take considerable time to process. Please hold.')
        ds = xr.concat([ds.where(ds['time.month'] < 6, drop=True).interpolate_na(dim='time'),
                        ds.where(ds['time.month'] > 6, drop=True).interpolate_na(dim='time')], dim='time').sortby('time')
        
        # tell user and return
        print('Missing pixels successfully filled. Proceeding.\n')
        return ds

    except Exception as e:
        print('Error occurred during flag_years_with_missing_wet_or_dry of type {0}. Stopping.'.format(e))
        raise e
        
# take ds, interpolate linear across years for original mk
def interp_original_trend(ds):

    # get missing dates
    print('Beginning linear interpolation of missing data for original mannkendall. Please wait.')
    nan_dates = np.setdiff1d(ds['time'], ds.dropna(dim='time', how='all')['time'])

    # drop any scenes where all nan
    ds = ds.dropna('time', how='all')

    if len(nan_dates) > 0:
        interp = ds.interp(time=nan_dates, kwargs={'fill_value': 'extrapolate'})
        ds = xr.concat([ds, interp], dim='time').sortby('time')
    else:
        print('No need to interpolate data. Skipping and proceeding.')
                
    # tell user and return
    print('Successfully interpolated missing data. Proceeding.')
    return ds

# take ds, interpolate linear across years, seasons/ months for seasonal mk
def interp_seasonal_trend(ds):
    
    def apply_interp(ds, nan_date, slice_month, interped_data):
        if ds['time.month'].values[0] == slice_month:
            interped_data.append(ds.interp(time=nan_date))

        return ds

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
    
        # get missing dates for wet, dry seperate
        print('\nBeginning linear interpolation of missing data for seasonal mannkendall. Please wait.')
        nan_years_interp, nan_dates_interp = [], []

        # determine if any whole years are nan, add to interp list, and remove
        nan_years_list = np.setdiff1d(ds.resample(time='1Y').median()['time.year'], 
                                      ds.resample(time='1Y').median().dropna(dim='time', how='all')['time.year'])
        nan_years_interp = ds.where(ds['time.year'].isin(nan_years_list), drop=True)['time']
        ds = ds.where(~ds['time.year'].isin(nan_years_list), drop=True)
        
        # determine if any remaining dates are nan, add to interp list, remove
        nan_dates_interp = np.setdiff1d(ds['time'], ds.dropna(dim='time', how='all')['time'])
        ds = ds.dropna(dim='time', how='all')      
        
        # interp across dates e.g. 2011 jan-dec, 2012 jan-dec
        if len(nan_dates_interp) > 0:
            print('Interpolating linearly within seasons/months, please wait.')
            ds = xr.concat([ds, ds.interp(time=nan_dates_interp)], dim='time').sortby('time')

        # loop interp missing year dates, loop ds, if month match, interp at date
        print('Interpolating linearly across years, please wait.')
        interped_data = []
        for nan_date in nan_years_interp:
            slice_month = nan_date.dt.month.values
            ds.groupby('time.month').apply(apply_interp, nan_date=nan_date, slice_month=slice_month, interped_data=interped_data)

        # put it back together
        if len(interped_data) > 0:
            ds = xr.concat([ds, xr.concat(interped_data, dim='time')], dim='time').sortby('time')

        # do a final check if anything was missed
        nan_final_interp = np.setdiff1d(ds['time'], ds.dropna(dim='time', how='all')['time'])
        if len(nan_final_interp):
            print('Interpolating missed dates. Please wait.')
            ds = ds.dropna(dim='time', how='all')
            ds = xr.concat([ds, ds.interp(time=nan_final_interp, kwargs={'fill_value': 'extrapolate'})], 
                           dim='time').sortby('time')

        # final clean up
        ds = ds.sortby('time')
        if 'month' in list(ds.data_vars):
            ds = ds.drop('month')
            
        # tell user and return
        print('Interpolation successful. Proceeding.\n')
        return ds

# take ds, interpolate (fill) missing pixels within each singular image for trend
def fill_pixels_trend(ds):
    try:
        print('Filling in missing pixels. Warning, this can take considerable time to process. Please hold.')
        ds = ds.interpolate_na(dim='time').sortby('time')
        
        # tell user and return
        print('Missing pixels successfully filled. Proceeding.\n')
        return ds

    except Exception as e:
        print('Error occurred during fill_pixels_trend of type {0}. Stopping.'.format(e))
        raise e
        
# take ds, interpolate linear across years for change vector analysis
def interp_change_vector(ds):

    # get missing dates
    print('Beginning linear interpolation of missing data for change vector analysis. Please wait.')
    nan_dates = np.setdiff1d(ds['time'], ds.dropna(dim='time', how='all')['time'])

    # drop any scenes where all nan
    ds = ds.dropna('time', how='all')

    if len(nan_dates) > 0:
        interp = ds.interp(time=nan_dates, kwargs={'fill_value': 'extrapolate'})
        ds = xr.concat([ds, interp], dim='time').sortby('time')
    else:
        print('No need to interpolate data. Skipping and proceeding.')
                
    # tell user and return
    print('Successfully interpolated missing data. Proceeding.')
    return ds      

# take ds, interpolate (fill) missing pixels within each singular image for change
def fill_pixels_change(ds):
    try:
        print('Filling in missing pixels. Warning, this can take considerable time to process. Please hold.')
        ds = ds.interpolate_na(dim='time').sortby('time')
        
        # tell user and return
        print('Missing pixels successfully filled. Proceeding.\n')
        return ds

    except Exception as e:
        print('Error occurred during fill_pixels_trend of type {0}. Stopping.'.format(e))
        raise e