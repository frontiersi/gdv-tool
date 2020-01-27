# scripts for statstical processors and summaries
import warnings
import math
import numpy as np
import pandas as pd
import xarray as xr
import sklearn.metrics as metrics
import pymannkendall as mk           # need pip to install
import ruptures as rpt               # need pip to install
from scipy.ndimage import convolve
from gdv_scripts import data

## HELPER FUNCTIONS ##
# check if any da in ds are all nan
def is_all_nan(da, all_nan_list):
    try:
        if da.isnull().all().values:
            datetime = pd.to_datetime(str(da['time'].values))
            datetime = datetime.strftime('%Y-%m-%d')
            all_nan_list.append(datetime)

        return da

    except Exception as e:
        print('Error occurred during is_all_nan of type {0}. Stopping.'.format(e))
        raise e      

# check if any da in ds are any nan but not all nan
def is_any_nan(da, any_nan_list):
    try:
        if not da.isnull().all().values:
            if da.isnull().any().values:
                datetime = pd.to_datetime(str(da['time'].values))
                datetime = datetime.strftime('%Y-%m-%d')
                any_nan_list.append(datetime)

        return da

    except Exception as e:
        print('Error occurred during is_any_nan of type {0}. Stopping.'.format(e))
        raise e
        
# do a majority filter on 3x3 window for 3d dataset
def apply_major_filter(ds):
    
    # tell user
    print('Applying majority filter. Please wait.')
    
    # get 3x3 window, convolve it and apply it
    kernel = np.ones((1, 3, 3))
    conv = lambda x: convolve(x, kernel, mode="wrap")
    ds = xr.apply_ufunc(conv, ds)
    
    return ds

# take a da, ds, and z crit, calc zscore, compare to crit, and flag in var outliers for wet/dry
def do_zscore_and_flag_wet_dry(da, ds, z_val):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        if da['time.month'] < 6:
            ds_season = ds.where(ds['time.month'] < 6, drop=True).groupby('time').mean(...)
        elif da['time.month'] > 6:
            ds_season = ds.where(ds['time.month'] > 6, drop=True).groupby('time').mean(...)
        else:
            return da

        # get test median, pop means, pop stdv (veg and mst index in each)
        medians = da.mean()
        p_means = ds_season.mean()
        p_stdvs = ds_season.std()

        # do veg z-score and check if > or < than crit val
        if p_means['veg_idx'] > 0:
            z = (medians['veg_idx'].values - p_means['veg_idx'].values) / p_stdvs['veg_idx'].values
            if (z > z_val) or (z < z_val * -1):                   
                da['veg_outlier'] = True
            else:
                da['veg_outlier'] = False

        # do mst z-score and check if > or < than crit val
        if p_means['mst_idx'] > 0:
            z = (medians['mst_idx'].values - p_means['mst_idx'].values) / p_stdvs['mst_idx'].values
            if (z > z_val) or (z < z_val * -1):
                da['mst_outlier'] = True
            else:
                da['mst_outlier'] = False    

        return da  

# take a da, ds, and z crit, calc zscore, compare to crit, and flag in var outliers for trend
def do_zscore_and_flag_trend(da, ds, z_val):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        # get mean of all times
        da_trend = ds.groupby('time').mean(...)

        # get test median, pop means, pop stdv (veg and mst index in each)
        medians = da.mean()
        p_means = da_trend.mean()
        p_stdvs = da_trend.std()

        # do veg z-score and check if > or < than crit val
        if p_means['veg_idx'] > 0:
            z = (medians['veg_idx'].values - p_means['veg_idx'].values) / p_stdvs['veg_idx'].values
            if (z > z_val) or (z < z_val * -1):                   
                da['veg_outlier'] = True
            else:
                da['veg_outlier'] = False

        return da  
    
# take a da, ds, and z crit, calc zscore, compare to crit, and flag in var outliers for change
def do_zscore_and_flag_change(da, ds, z_val):
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        # get mean of all times
        da_trend = ds.groupby('time').mean(...)

        # get test median, pop means, pop stdv (veg and mst index in each)
        medians = da.mean()
        p_means = da_trend.mean()
        p_stdvs = da_trend.std()

        # do veg z-score and check if > or < than crit val
        if p_means['veg_idx'] > 0:
            z = (medians['veg_idx'].values - p_means['veg_idx'].values) / p_stdvs['veg_idx'].values
            if (z > z_val) or (z < z_val * -1):                   
                da['veg_outlier'] = True
            else:
                da['veg_outlier'] = False

        # do mst z-score and check if > or < than crit val
        if p_means['brt_idx'] > 0:
            z = (medians['brt_idx'].values - p_means['brt_idx'].values) / p_stdvs['brt_idx'].values
            if (z > z_val) or (z < z_val * -1):
                da['brt_outlier'] = True
            else:
                da['brt_outlier'] = False    

        return da  
    
# look up function to find orthogonal polynomial coefficient ss, constant and list of coefficients
def build_opc_params(num_scenes):
    if num_scenes == 3:
        coeff_list = [-1, 0, 1]
    elif num_scenes == 4:
        coeff_list = [-3, -1, 1, 3]
    elif num_scenes == 5:
        coeff_list = [-2, -1, 0, 1, 2]
    elif num_scenes == 6:
        coeff_list = [-5, -3, -1, 1, 3, 5]
    elif num_scenes == 7:
        coeff_list = [-3, -2, -1, 0, 1, 2, 3]
    elif num_scenes == 8:
        coeff_list = [-7, -5, -3, -1, 1, 3, 5, 7]
    elif num_scenes == 9:
        coeff_list = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    elif num_scenes == 10:
        coeff_list = [-9, -7, -5, -3, -1, 1, 3, 5, 7, 9] 
    elif num_scenes == 11:
        coeff_list = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
    elif num_scenes == 12:
        coeff_list = [-11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11] 
    elif num_scenes == 13:
        coeff_list = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
    elif num_scenes == 14:
        coeff_list = [-13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13]
    elif num_scenes == 15:
        coeff_list = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
    elif num_scenes == 16:
        coeff_list = [-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15]
    elif num_scenes == 17:
        coeff_list = [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]           
    elif num_scenes == 18:
        coeff_list = [-17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17]   
    elif num_scenes == 19:
        coeff_list = [-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    elif num_scenes == 20:
        coeff_list = [-19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
    elif num_scenes == 21:
        coeff_list = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    elif num_scenes == 22:
        coeff_list = [-21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    elif num_scenes == 23:
        coeff_list = [-11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    elif num_scenes == 24:
        coeff_list = [-23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]
    elif num_scenes == 25:
        coeff_list = [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    elif num_scenes == 26:
        coeff_list = [-25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]           
    elif num_scenes == 27:
        coeff_list = [-13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    elif num_scenes == 28:
        coeff_list = [-27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27]    
    elif num_scenes == 29:
        coeff_list = [-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    elif num_scenes == 30:
        coeff_list = [-29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]
    elif num_scenes == 31:
        coeff_list = [-15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    elif num_scenes == 32:
        coeff_list = [-31, -29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
    elif num_scenes == 33:
        coeff_list = [-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    elif num_scenes == 34:
        coeff_list = [-33, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33]
    elif num_scenes == 35:
        coeff_list = [-17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    elif num_scenes == 36:
        coeff_list = [-35, -33, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35]
    elif num_scenes == 37:
        coeff_list = [-18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    elif num_scenes == 38:
        coeff_list = [-37, -35, -33, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37]
    elif num_scenes == 39:
        coeff_list = [-19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            
    # calculate sum of squares and constant
    if len(coeff_list) > 0:
        ss = sum(map(lambda i : i * i, coeff_list))   # sum of squares
        constant = 2 if num_scenes % 2 == 0 else 1    # constant, if num scenes even 2, else 1

        return coeff_list, ss, constant
        
    else:
        print('Cannot create opc parameters, less than 3 dates of data in dataset. Please fetch more data.')
        return [], None, None
    
# get n% percentile of veg, mst, create igt masks, give user basic info and display
def mask_via_percentile(da, percentile_val):
    
    # build mask from percentile
    mask = xr.where(da > np.nanpercentile(da, q=percentile_val), 1, 0)
    
    return mask

# get steadiest of slope pixels via percentile close to 0 either side
def get_steadiest_opc_slope(da, slope_steadiness_perecentile):
    try:
        neg = np.nanpercentile(da.where(da <= 0), q=(100 - slope_steadiness_perecentile))
        pos = np.nanpercentile(da.where(da >= 0), q=slope_steadiness_perecentile)
    except:
        print('Problem during get_steadiest_opc_slope. Aborting.')
        return None
    
    if neg and pos:
        return xr.where((da > neg) & (da < pos), 1, 0)
    else:
        return None
    
# mask invariant greennest, moistest sites
def make_invariant_sites(mask_highest, mask_steadiest):
    
    # mask where greenest (for example) and steadiesy
    mask = xr.where(mask_highest & mask_steadiest, 1, 0)
    
    return mask
 
# reduce invariant sites to user defined if greater than
def reduce_invariant_sites_mask(da, mask, max_num_sites):

    # get lowest highest index veg or mst value within invariant sites
    def get_lowest_value_in_invariant_sites(da, max_num_sites):
        vals = da.stack(v=('x', 'y')).values            # stack into 1d
        vals = vals[np.isfinite(vals)]                  # remove nan
        vals = vals[np.argsort(vals)[-max_num_sites:]]  # sort low to high, reduce to user param

        return np.nanmin(vals)

    # check if number of user chosen sites < than what has come back
    if max_num_sites <= mask.sum():

        # tell user
        msg = '> Number of invariant sites ({0}) greater than maximum allowed ({1}). Reducing down.'
        print(msg.format(mask.sum().values, max_num_sites))

        # get lowest highest index veg or mst value within invariant sites
        min_max_val = get_lowest_value_in_invariant_sites(da=da.where(mask, drop=True), max_num_sites=max_num_sites)

        # reduce invariant target mask to meet max num sites, getting highest values
        return xr.where(mask & (da >= min_max_val), 1, 0)
        
    else:
        print('> Number of invariant sites lower than user defined. Skipping reduction.')
        return mask
    
# create empty mannkendall dataset for holding mk results
def prepare_empty_mk_dataset(scene_like_threshed):

    if not scene_like_threshed:
        print('Need to provide a likelihood dataset. Could not run. Try again.')
        return None
        
    try:
        # create empty dataset for seasonal mk
        print('> Preparing empty mannkendall dataset, hold on...')
        empty = xr.full_like(scene_like_threshed['likelihood'], np.nan)
        mk_result = scene_like_threshed.drop(scene_like_threshed.data_vars)

        # fill with new empty nan arrays
        mk_result['tau'] = empty
        mk_result['p'] = empty
        mk_result['z'] = empty
        mk_result['slope'] = empty
        mk_result['s'] = empty

        # tell user and return
        print('> Prepared mannkendall dataset successfully. Proceeding.')
        return mk_result

    except:
        print('Problem during mannkendall dataset creation. Please check.')
        return None

# mask out non-gdv areas from dataset to speed up processing
def mask_non_gdv_data(ds, scene_like_threshed):

    if not ds or not scene_like_threshed:
        print('No dataset or likelihood map dataset provided. Cannot run. Please check.')
        return None
        
    # tell user
    print('> Masking out non-groundwater dependent vegetation pixels to speed up processing. Please hold.')

    # mask out non-gdv pixels
    ds = ds.where(np.isfinite(scene_like_threshed['likelihood']), np.nan)
        
    # tell user and return
    print('> Dataset sucessfully masked to GDV pixels. Proceeding.')
    return ds

# stack pixels into pixel vectors
def stack_and_prepare_pixels(ds):
        
    # basic checks
    if not ds:
        print('No vegetation dataset provided. Cannot run. Please check.')
        return None
            
    # tell user and stack
    print('> Beginning vector stacking process, please hold...')
    stack = ds.stack(v=('x', 'y'))                # stack into z-index
    stack = stack.dropna(dim='v', how='all')      # drop non-gdv (all nan vectors)

    # setup vars
    stack['p'] = np.nan        # p-value of the significance test
    stack['tau'] = np.nan      # kendall tau value, rank correlation
    stack['z'] = np.nan        # normalized test statistics
    stack['slope'] = np.nan    # sen's slope value
    stack['s'] = np.nan        # mannkendalls score

    # tell user and return
    print('> Data stacked into vectors successfully! Mann-Kendall trend analysis now ready.')
    return stack

# helper apply func to do mannkendall on per vector basis
def apply_mk_to_vector(vector, mk_type, sig_only, trend_type, period):

    # fill a vector
    def fill_a_vector(vector, result):
        vector['tau'] = result.Tau        # kendall tau value, rank correlation
        vector['p'] = result.p            # p-value of the significance test
        vector['z'] = result.z            # normalized test statistics
        vector['slope'] = result.slope    # sen's slope value
        vector['s'] = result.s            # mannkendalls score

        return vector

    if not np.isnan(vector['veg_idx'].data).all():
        try:
            # check if original or seasonal mannkendall
            if mk_type == 'original':
                result = mk.original_test(vector['veg_idx'])
            elif mk_type == 'seasonal':
                result = mk.seasonal_test(vector['veg_idx'], period=period)
            else:
                return vector

            # return mk only if meets user params (increasing/decreasing/sig/notsig)
            if sig_only:
                if result.h and (result.trend == trend_type or trend_type == 'all'):
                    return fill_a_vector(vector, result)
            elif not sig_only:
                if result.trend == trend_type or trend_type == 'all':
                    return fill_a_vector(vector, result)
            else:
                return vector    
        except:
            print('Problem occured at a pixel during mannkendall. Skipping current pixel.')
            return vector

    return vector       

# perform original, seasonal mannkendall
def perform_mk(stacked, mk_type, sig_only, trend_type, period):
        
    # basic checks
    if not stacked:
        print('No stacked vegetation dataset provided. Cannot run. Please check.')
        return None

    # tell user
    print('Beginning mannkendall analysis. This can take awhile. Go get a coffee!')

    # perform our mannkendall, depending on user params
    if mk_type == 'original':
        print('> Original mannkendall analysis being applied, hold tight.')
    elif mk_type == 'seasonal':
        print('> Seasonal mannkendall analysis with period of {0} being applied, hold tight...'.format(period))

    # do the mk!
    stacked = stacked.groupby('v').progress_apply(apply_mk_to_vector, mk_type=mk_type, sig_only=sig_only, 
                                                  trend_type=trend_type, period=period)
    # tell user
    print('> Successfully completed the mannkendall!')

    # tell user and clean up
    print('> Cleaning up some bits and pieces, please hold.')
    unstacked = stacked.unstack('v')
    unstacked = unstacked.transpose('time', 'y', 'x')

    # tell user and return
    print('> Finished mannkendal and cleaning procedure! Proceeding.')
    return unstacked

# do final transfer and clean up of vector values over to mk_result dataset
def transfer_unstack_to_mk_result(mk_result, unstacked):
        
    # basic checks
    if not mk_result or not unstacked:
        print('No mannkendall result or stacked vegetation dataset provided. Cannot run. Please check.')
        return None

    # tell user
    print('> Transfering mannkendall result over to final dataset.')
        
    # do update from vector, remove uneeded vars, dims
    try:
        # update empty mk result in place (via coords) for each vector
        print('> Updating empty mannkendall result in-place using vectors, hold tight.')
        mk_result = mk_result.update(unstacked)

        # drop any of the forbidden vars
        mk_result = mk_result.drop(['time', 'veg_idx'], errors='ignore')
        
        return mk_result

    except:
        print('Problem during final mk result clean. Please check.')
        return None
    
# seperate change dataset into a baseline (all time median) and ds of comparison dates
def seperate_baseline_and_compare_ds(ds, base_dates, comp_dates):
    
    # basic check
    if not ds:
        print('No dataset provided. Aborting.')
        return None
    elif not base_dates or not comp_dates:
        print('No baseline or comparison years provided. Aborting.')
        return ds
    elif len(ds['time.year']) == 0:
        print('No years of data in dataset. Please include more years in your query above.')
        return ds
    
    try:
        # create all time median for selected baseline dates
        print('> Generating baseline dataset for date range: {0} to {1}'.format(base_dates[0], base_dates[1]))
        base_ds = ds.sel(time=slice(base_dates[0], base_dates[1])).median('time')

        # create dataset of all dates in compare range
        print('> Generating comparison data for date range: {0} to {1}'.format(comp_dates[0], comp_dates[1]))
        comp_ds = ds.sel(time=slice(comp_dates[0], comp_dates[1]))

        return base_ds, comp_ds
    
    except:
        print('Something happened during seperate_baseline_and_compare_ds. Please check.')

# prepare comparison dataset for holding cva results
def prepare_cva_comparison_dataset(comp_ds):
    try:
        # create empty da
        empty = xr.full_like(comp_ds['brt_idx'], np.nan)

        # make new vars for cva
        comp_ds['angout'] = empty
        comp_ds['magout'] = empty
        comp_ds['magnitude'] = empty
        comp_ds['brt_diff'] = empty
        comp_ds['veg_diff'] = empty
        comp_ds['sel'] = empty
            
        return comp_ds
    except:
        print('Something went wrong during cva preparation. Please check.')
        return comp_ds

# restrict cva result to specific angle range set by user
def reduce_cva_to_specific_angles(cva_ds, min_ang, max_ang):
    
    if not cva_ds:
        print('No cva dataset provided. Aborting.')
        return None
    elif not min_ang or not max_ang:
        print('No minimum/maximum angle provided. Aborting.')
        return cva_ds
    
    # restrict to specific angle range e.g. >= 90 to <=180 is vegetation decline
    print('> Reducing cva angles from {0} to {1} degrees'.format(min_ang, max_ang))
    cva_ds['angout'] = cva_ds['angout'].where(((cva_ds['angout'] >= min_ang) & (cva_ds['angout'] <= max_ang)))

    # reduce magnitude now to only those angles that exist
    print('> Restricting magnitudes to only those within selected angles')
    cva_ds['magout'] = cva_ds['magout'].where(np.isfinite(cva_ds['angout']))
        
    return cva_ds

# mask cva results to gdv likelihood thresholded layer
def mask_cva_to_gdv_likelihood(cva_ds, likelihood):
    try:
        mask_gdv = xr.where(np.isfinite(likelihood), True, False)
        cva_ds['angout'] = cva_ds['angout'].where(mask_gdv)
        cva_ds['magout'] = cva_ds['magout'].where(mask_gdv)
        
        return cva_ds
    
    except:
        print('Problem during mask_cva_to_gdv_likelihood. Skipping mask process.')
        return cva_ds
    
## MAIN FUNCTIONS ##
# generate summary statistics for wet and dry season-type dataset
def calc_summary_stats_wet_and_dry(ds):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) <= 2:
        print('Less than three years of data in dataset. Please include more years in your query above.')
        return None
    elif len(ds.groupby('time.month')) != 2:
        print('Num of months of data in dataset not equal to 2. Please include more years in your query above.')
        return None
    
    try:
        # number of dates, years num wet, num dry
        total_dates = len(ds['time'])
        total_years = len(ds['time'].groupby('time.year'))
        num_wet_dates = len(ds['time'].where(ds['time.month'] < 6, drop=True))
        num_dry_dates = len(ds['time'].where(ds['time.month'] > 6, drop=True))

        # get dates that are all nan
        all_nan_dates_veg, all_nan_dates_mst = [], []
        ds['veg_idx'].groupby('time').apply(is_all_nan, all_nan_list=all_nan_dates_veg)
        ds['mst_idx'].groupby('time').apply(is_all_nan, all_nan_list=all_nan_dates_mst)

        # get any other dates that have any (not all) nan
        any_nan_dates_veg, any_nan_dates_mst = [], []
        ds['veg_idx'].groupby('time').apply(is_any_nan, any_nan_list=any_nan_dates_veg)
        ds['mst_idx'].groupby('time').apply(is_any_nan, any_nan_list=any_nan_dates_mst)
        
        # get years missing a season or containing more than 2 seasons (wet and dry)
        num_seasons = ds['time'].groupby('time.year').count('time')
        years_miss_season = num_seasons['year'].where(num_seasons < 2, drop=True).values
        years_extra_season = num_seasons['year'].where(num_seasons > 2, drop=True).values

        # get mean of wet, dry, veg, mst
        da_mean_wet = ds.where(ds['time.month'] < 6).mean()
        da_mean_dry = ds.where(ds['time.month'] > 6).mean()

        # get median of wet, dry, veg, mst
        da_medn_wet = ds.where(ds['time.month'] < 6).median()
        da_medn_dry = ds.where(ds['time.month'] > 6).median()

        # get std dev of wet, dry, veg, mst
        da_stdv_wet = ds.where(ds['time.month'] < 6).std()
        da_stdv_dry = ds.where(ds['time.month'] > 6).std()

        # relay time summary to user
        print('=== Summary dataset statistics ===')
        print('> Total number of dates in dataset: {0}'.format(total_dates))
        print('> Total number of years in dataset: {0}'.format(total_years))
        print('> Number of wet (DJF) dates in dataset: {0}'.format(num_wet_dates))
        print('> Number of dry (SON) dates in dataset: {0}'.format(num_dry_dates))

        # relay missing data summary to user
        print('> Years with missing season(s): {0}'.format(years_miss_season))
        print('> Years with extra season(s): {0}'.format(years_extra_season))
        print('> Completely missing vegetation data on dates: {0}'.format(all_nan_dates_veg))
        print('> Completely missing moisture data on dates: {0}'.format(all_nan_dates_mst))
        print('> Partially missing vegetation data on dates: {0}'.format(any_nan_dates_veg))
        print('> Partially missing moisture data on dates: {0}'.format(any_nan_dates_mst))

        # relay mean stats to user
        msg = '> Mean values for entire wet season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_mean_wet['veg_idx'].values, 3), np.round(da_mean_wet['mst_idx'].values, 3)))
        msg = '> Mean values for entire dry season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_mean_dry['veg_idx'].values, 3), np.round(da_mean_dry['mst_idx'].values, 3)))

        # relay median stats to user
        msg = '> Median values for entire wet season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_medn_wet['veg_idx'].values, 3), np.round(da_medn_wet['mst_idx'].values, 3)))
        msg = '> Median values for entire dry season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_medn_dry['veg_idx'].values, 3), np.round(da_medn_dry['mst_idx'].values, 3)))

        # relay std dev stats to user
        msg = '> Standard deviation values for entire wet season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_stdv_wet['veg_idx'].values, 3), np.round(da_stdv_wet['mst_idx'].values, 3)))
        msg = '> Standard deviation values for entire dry season: vege: {0}, moist: {1}'
        print(msg.format(np.round(da_stdv_dry['veg_idx'].values, 3), np.round(da_stdv_dry['mst_idx'].values, 3)))

        # close up
        print('=== Summary dataset statistics ===\n')

    except Exception as e:
        print('Error occurred during calc_summary_stats_wet_and_dry of type {0}. Skipping.'.format(e))

# calc if any veg, mst are z-score sig for wet, dry season dataset, flag, show stats, reset to all nan if outlier
def zscore_flag_show_reset_wet_dry(ds, z_val, reset):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) <= 2:
        print('Less than three years of data in dataset. Please include more years in your query above.')
        return None
    elif len(ds.groupby('time.month')) != 2:
        print('Num of months of data in dataset not equal to 2. Please include more years in your query above.')
        return None
    elif z_val == 0:
        print('Z-score critical value is invalid. Please ensure it is not 0.')
        return None
    
    try:
        # tell user
        print('Performing Z-score test to determine data outliers. Please wait.')  
        
        # apply z-score and flag outliers
        ds = ds.groupby('time').apply(do_zscore_and_flag_wet_dry, ds=ds, z_val=z_val)
    
        # see if vars in ds and set all flagged to nan
        if all(elem in list(ds.data_vars) for elem in ['veg_outlier', 'mst_outlier']):
            
            # get dates for printing
            veg_outs = [str(e['time'].dt.strftime('%Y-%m-%d').values) for e in ds['time'].where(ds['veg_outlier'], drop=True)]
            mst_outs = [str(e['time'].dt.strftime('%Y-%m-%d').values) for e in ds['time'].where(ds['mst_outlier'], drop=True)]
            
            # show summary
            print('=== Z-Score Results using score of {0} ==='.format(z_val))
            print('> Vege Outlier Dates: {0}'.format(veg_outs))
            print('> Moisture Outlier Dates: {0}'.format(mst_outs))
            print('=== Z-Score Results using score of {0} ===\n'.format(z_val))
            
            # if reset, set outlier scenes to all nan
            if reset:
                ds = ds.where((ds['veg_outlier'] == False) & (ds['mst_outlier'] == False))
                
        else:
            print('No vegetation and/or moisture outlier mask. Skipping and proceeding.')
                  

        # drop outliet mask vars
        ds = ds.drop(['veg_outlier', 'mst_outlier'], errors='ignore')
        return ds

    except Exception as e:
        print('Error occurred during zscore_flag_show_reset_wet_dry of type {0}. Stopping.'.format(e))
        raise e   

# generate opc slope function, take ds and var name (e.g. veg_idx). dry period only.
def build_opc_slope_dry(ds, var):

    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) <= 2:
        print('Less than three years of data in dataset. Please include more years in your query above.')
        return None
    elif len(ds[var].where(ds['time.month'] > 6, drop=True)) <= 2:
        print('Less than three scenes of data in dataset. Please include more years in your query above.')
        return None
    elif len(ds.groupby('time.month')) != 2:
        print('Num of months of data in dataset not equal to 1. Please include more years in your query above.')
        return None

    # get opc coefficient parameters
    num_scenes = len(ds[var].where(ds['time.month'] > 6, drop=True))
    coeff_list, ss, constant = build_opc_params(num_scenes=num_scenes)
    
    # tell user
    msg = '> Generating dry season slope for var: {0}. Number of coefficients: {1}. Sum of squares: {2}. Constant: {3}.'
    print(msg.format(var, len(coeff_list), ss, constant))

    # get veg opc slope
    if not coeff_list or not ss or not constant:
        print('No coefficient list, sum of squares, or constant. Aborting.')
        return None

    try:
         # multiply each date in ds by coefficient list (3d array x 1d array), divide ss and * constant, return
        slope = (ds[var].where(ds['time.month'] > 6, drop=True) * np.array(coeff_list)[:, None, None]).sum('time') / ss * constant
        return slope
    
    except Exception as e:
        print('Error occurred during build_opc_slope_dry of type {0}. Stopping.'.format(e))
        return None
    
# standardise veg, mst da to dry season using invariant targets, low and high inflection. used in apply func
def perform_increase_sigmoid(da, mask_veg, mask_mst, high_percentile, low_inflect_veg, low_inflect_mst):
    try:
        # get high inflection point via percentile for veg, mst
        high_inflect_veg = np.nanpercentile(da['veg_idx'].where(mask_veg), q=high_percentile)
        high_inflect_mst = np.nanpercentile(da['mst_idx'].where(mask_mst), q=high_percentile)

        # make masks for veg, mst
        low_mask_veg, high_mask_veg = xr.where(da['veg_idx'] > low_inflect_veg, 1, 0), xr.where(da['veg_idx'] < high_inflect_veg, 1, 0)
        low_mask_mst, high_mask_mst = xr.where(da['mst_idx'] > low_inflect_mst, 1, 0), xr.where(da['mst_idx'] < high_inflect_mst, 1, 0)   

        # do increasing fuzzy sigmoidal for veg, mst
        stand_veg = np.square(np.cos((1 - ((da['veg_idx'] - low_inflect_veg) / (high_inflect_veg - low_inflect_veg))) * (math.pi / 2)))
        stand_mst = np.square(np.cos((1 - ((da['mst_idx'] - low_inflect_mst) / (high_inflect_mst - low_inflect_mst))) * (math.pi / 2)))

        # keep veg in inflection mask, else 0 (low) or 1 (high)
        stand_veg = stand_veg.where(low_mask_veg, 0.0)
        stand_veg = stand_veg.where(high_mask_veg, 1.0)

        # keep mst in inflection mask, else 0 (low) or 1 (high)
        stand_mst = stand_mst.where(low_mask_mst, 0.0)
        stand_mst = stand_mst.where(high_mask_mst, 1.0)

        # add to da
        da['veg_idx'] = stand_veg
        da['mst_idx'] = stand_mst
        return da
    
    except:
        print('Problem during perform_increase_sigmoid. Returning original data array.')
        return da
    
# standardise veg, mst stability da. used in apply func
def perform_inc_dec_sigmoid(da):
    try:
        # make left and right curve masks, veg and mst
        mask_left_veg, mask_right_veg = xr.where(da['veg_stb'] <= 1, 1, 0), xr.where(da['veg_stb'] > 1, 1, 0)
        mask_left_mst, mask_right_mst = xr.where(da['mst_stb'] <= 1, 1, 0), xr.where(da['mst_stb'] > 1, 1, 0)

        # do left sigmoidal (increasing) for veg, mst
        left_veg = np.square(np.cos((1 - ((da['veg_stb'] - 0) / (1 - 0))) * (math.pi / 2)))      
        left_mst = np.square(np.cos((1 - ((da['mst_stb'] - 0) / (1 - 0))) * (math.pi / 2)))
        
        # now right sigmoidal (decreasing) for veg, mst
        right_veg = np.square(np.cos(((da['veg_stb'] - 1) / (2 - 1)) * (math.pi / 2)))        
        right_mst = np.square(np.cos(((da['mst_stb'] - 1) / (2 - 1)) * (math.pi / 2)))
        
        # mask veg areas outside of left, right range
        left_veg = left_veg.where(mask_left_veg, 0.0)
        right_veg = right_veg.where(mask_right_veg, 0.0)
        
        # mask mst areas outside of left, right range
        left_mst = left_mst.where(mask_left_mst, 0.0)
        right_mst = right_mst.where(mask_right_mst, 0.0)       
        
        # combine left and right sides
        da['veg_stb'] = left_veg + right_veg
        da['mst_stb'] = left_mst + right_mst

        return da
    
    except:
        print('Problem during perform_inc_dec_sigmoid. Returning original data array.')
        return da
     
# used in apply, will divide by percentile for every scene
def perform_divide_by_percentile(da, mask_greenest, high_percentile):
    
    max_value = np.nanpercentile(da['veg_idx'].where(mask_greenest), q=high_percentile)
    da['veg_idx'] = da['veg_idx'] / max_value
    da['veg_idx'] = xr.where(da['veg_idx'] > 1, 1, da['veg_idx'])
    da['veg_idx'] = xr.where(da['veg_idx'] < 0, 0, da['veg_idx'])
    
    return da

# used in apply, will divide by percentile for every scene (for change vector)
def perform_divide_by_percentile_change_vector(da, mask_greenest, mask_brightest, high_percentile):
    
    # get max veg, brt values using percentile
    max_veg_value = np.nanpercentile(da['veg_idx'].where(mask_greenest), q=high_percentile)
    max_brt_value = np.nanpercentile(da['brt_idx'].where(mask_brightest), q=high_percentile)
    
    # standardise using basic division for veg
    da['veg_idx'] = da['veg_idx'] / max_veg_value
    da['veg_idx'] = xr.where(da['veg_idx'] > 1, 1, da['veg_idx'])
    da['veg_idx'] = xr.where(da['veg_idx'] < 0, 0, da['veg_idx'])
    
    # standardise using basic division for brt
    da['brt_idx'] = da['brt_idx'] / max_brt_value
    da['brt_idx'] = xr.where(da['brt_idx'] > 1, 1, da['brt_idx'])
    da['brt_idx'] = xr.where(da['brt_idx'] < 0, 0, da['brt_idx'])
    
    return da
        
# calculate stability between pixels between wet and dry seasons
def calc_seasonal_stability(da, stability_list):
    
    if len(da['time'] == 2):
        
        try:
            # get year (int) of current da
            group = list(da.groupby('time.year').groups)
            year = str(group[0])
        except:
            return da

        if year:
            wet = da.where(da['time.month'] < 6, drop=True).squeeze(drop=True)
            dry = da.where(da['time.month'] > 6, drop=True).squeeze(drop=True)

            # generate stability veg, mst
            stability = xr.merge([(wet['veg_idx'] - dry['veg_idx']), (wet['mst_idx'] - dry['mst_idx'])])
            stability = stability.assign_coords(time=year)
            stability_list.append(stability)
        
        return da
    
# mask stability maps to veg, mst > 75
def mask_seasonal_stability(stability, ds, percentile):
    
    def do_mask(da, ds, percentile):
        
        try:
            year = str(da['time'].values)

            if year:
                # get mean of two months in year
                scene_mean = ds.sel(time=year).groupby('time.year').mean()

                # mask where > percentile
                mask = xr.where(scene_mean > np.nanpercentile(scene_mean, q=percentile), 1, 0)
                mask = mask.squeeze('year', drop=True)
                da = da.where(mask, 0.0)   
            else:
                return da
        except:
            return da

        return da
    
    # do masking process for veg, mst
    stability['veg_stb'] = stability['veg_stb'].groupby('time').apply(do_mask, ds=ds['veg_idx'], percentile=percentile)
    stability['mst_stb'] = stability['mst_stb'].groupby('time').apply(do_mask, ds=ds['mst_idx'], percentile=percentile)
    
    
    return stability

# generate likelihood maps using ahp weights
def generate_likelihood_map(da, stability, likelihood_list, params):
        
    if len(da['time'] == 2):

        # get year (int) of current da
        group = list(da.groupby('time.year').groups)
        year = str(group[0])

        if year:
                       
            # get stability for current year, split da into wet, dry
            stb = stability.sel(time=year)
            wet = da.where(da['time.month'] < 6, drop=True).squeeze(drop=True)
            dry = da.where(da['time.month'] > 6, drop=True).squeeze(drop=True)
            
            # weight differently depending on index with swir or not
            if params.veg_idx in ['ndvi', 'savi']:
                wet['veg_idx'] = wet['veg_idx'] * 0.162613224   # orig 0.104621642
                wet['mst_idx'] = wet['mst_idx'] * 0.104621642   # orig 0.162613224
                dry['veg_idx'] = dry['veg_idx'] * 0.362481346
                dry['mst_idx'] = dry['mst_idx'] * 0.284883443
                stb['veg_stb'] = stb['veg_stb'] * 0.057825868
                stb['mst_stb'] = stb['mst_stb'] * 0.027574478
            
            else:
                wet['veg_idx'] = wet['veg_idx'] * 0.212613224 
                wet['mst_idx'] = wet['mst_idx'] * 0.054621642
                dry['veg_idx'] = dry['veg_idx'] * 0.462481346
                dry['mst_idx'] = dry['mst_idx'] * 0.184883443
                stb['veg_stb'] = stb['veg_stb'] * 0.057825868
                stb['mst_stb'] = stb['mst_stb'] * 0.027574478
                        
            # sum weights
            like = wet['veg_idx'] + dry['veg_idx'] + stb['veg_stb'] + wet['mst_idx'] + dry['mst_idx'] + stb['mst_stb']
            
            # append date (year) on, append to likelihood list
            like = like.assign_coords(time=year)
            likelihood_list.append(like)
            
        else:
            print('Warning, not enough years in current generate_likelihood_map iteration. Skipping.')

        return da
    
# perform roc analysis, return  tpr, fpr, thresholds, and auc
def perform_roc(actual_list, prediction_list):
   
    if len(actual_list) == len(prediction_list):
        print('Beginning ROC analysis. Please hold.')
        
        fpr, tpr, thresholds = metrics.roc_curve(actual_list, prediction_list)
        auc = metrics.roc_auc_score(actual_list, prediction_list)
        cut_off = thresholds[np.argmax(tpr - fpr)]
        
        # tell user and return
        print('Successfully completed ROC analysis. Proceeding.')
        return fpr, tpr, thresholds, auc, cut_off
    else:
        return None, None, None, None, None
    
# threshold da using standard deviation size
def thresh_via_stdv(da, num_of_std_devs):
    
    # get standard dev and mean
    mean = np.nanmean(da)
    stdv = np.nanstd(da)

    # correct standard dev via user input
    thresh = mean + (stdv * num_of_std_devs)

    # now, threshold likelihood median to std dev
    da = da.where(da > thresh)
    
    # tell user and return
    print('> Mean: {0}. Standard Deviation: {1}. Final threshold: {2}'.format(round(mean, 3), round(stdv, 3), round(thresh, 3)))
    return da

# threshold da using groundtruth shapefile
def thresh_via_shape(da, shp_path):
    
    # vars
    thresh, shp, fpr, tpr, auc, cutoff = None, None, None, None, None, None
    
    try:
        print('\nThresholding groundwater dependent vegetation likelihood using groundtruthed shapefile. Please wait.')
        shp = data.open_and_prepare_shapefile_for_roc(da=da, shapefile_path=shp_path)  
        actual_list, prediction_list = shp['GDV_ACT'].ravel(), shp['GDV_PRED'].ravel()

        if len(actual_list) > 0 and len(prediction_list) > 0:
            fpr, tpr, thresholds, auc, cutoff = perform_roc(actual_list=actual_list, prediction_list=prediction_list)
                    
            if cutoff:
                thresh = da.where(da > cutoff)
            else:
                print('Error with groundtruthing data actual vs prediction count. Aborting.')
                
        # tell user and return        
        print('Successfully thresholded groundwater dependent vegetation likelihood. Proceeding.')
        return thresh, shp, fpr, tpr, auc, cutoff
            
    except:
        print('Error thresholding shapefile-based thresholding. Aborting.')
        return thresh, shp, fpr, tpr, auc, cutoff
            
# generate summary statistics for trend dataset
def calc_summary_stats_trend(ds):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) < 5:
        print('Less than five years of data in dataset. Please include more years in your query above.')
        return None
    
    try:
        # number of dates, years num wet, num dry
        total_dates = len(ds['time'])
        total_years = len(ds['time'].groupby('time.year'))

        # get dates that are all nan
        all_nan_dates_veg = []
        ds['veg_idx'].groupby('time').apply(is_all_nan, all_nan_list=all_nan_dates_veg)

        # get any other dates that have any (not all) nan
        any_nan_dates_veg = []
        ds['veg_idx'].groupby('time').apply(is_any_nan, any_nan_list=any_nan_dates_veg)
        
        # relay time summary to user
        print('=== Summary dataset statistics ===')
        print('> Total number of dates in dataset: {0}'.format(total_dates))
        print('> Total number of years in dataset: {0}'.format(total_years))

        # relay missing data summary to user
        print('> Completely missing vegetation data on dates: {0}'.format(all_nan_dates_veg))
        print('> Partially missing vegetation data on dates: {0}'.format(any_nan_dates_veg))

        # close up
        print('=== Summary dataset statistics ===\n')

    except Exception as e:
        print('Error occurred during calc_summary_stats_trend of type {0}. Skipping.'.format(e))
        
# calc if any veg, mst are z-score sig for trend dataset, flag, show stats, reset to all nan if outlier
def zscore_flag_show_reset_trend(ds, z_val, reset):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) < 5:
        print('Less than five years of data in dataset. Please include more years in your query above.')
        return None
    elif z_val == 0:
        print('Z-score critical value is invalid. Please ensure it is not 0.')
        return None
    
    try:
        # tell user
        print('Performing Z-score test to determine data outliers. Please wait.')  
        
        # apply z-score and flag outliers
        ds = ds.groupby('time').apply(do_zscore_and_flag_trend, ds=ds, z_val=z_val)
            
        # see if vars in ds and set all flagged to nan
        if all(elem in list(ds.data_vars) for elem in ['veg_outlier']):
            
            # get dates for printing
            veg_outs = [str(e['time'].dt.strftime('%Y-%m-%d').values) for e in ds['time'].where(ds['veg_outlier'], drop=True)]
            
            # show summary
            print('=== Z-Score Results using score of {0} ==='.format(z_val))
            print('> Vege Outlier Dates: {0}'.format(veg_outs))
            print('=== Z-Score Results using score of {0} ===\n'.format(z_val))
            
            # if reset, set outlier scenes to all nan
            if reset:
                ds = ds.where(ds['veg_outlier'] == False)
                
        else:
            print('No vegetation and/or moisture outlier mask. Skipping and proceeding.')
                  
        # drop outliet mask vars
        ds = ds.drop(['veg_outlier'], errors='ignore')
        return ds

    except Exception as e:
        print('Error occurred during zscore_flag_show_reset_wet_dry of type {0}. Stopping.'.format(e))
        raise e   
        
# do basic green site standardisation process sites for trends
def do_standardisation_indices_trend(ds, params):
    
    # basic check
    if not ds:
        print('Error occured after compute, no dataset exists. Aborting.')
        return None
    elif not params.check_trend_standardisation_params():
        print('No standardisation parameters set. Please ensure they are set above.')
        return ds
    elif len(ds['time.year']) < 5:
        print('Less than five years of data in dataset. Please include more years in your query above.')
        return ds
            
    # tell user
    print('\nGenerating standardisation targets for vegetation for original mannkendall. Please wait.')
    
    # get median all time for veg
    med_alltime = ds.median('time')
    
    # get n% percentile of veg, create igt masks, give user basic info and display
    print('> Getting highest valued vegetation pixels and reducing to max number of sites. Please wait.') 
    mask_green = mask_via_percentile(da=med_alltime['veg_idx'], percentile_val=params.green_percentile)
    mask_greenest = reduce_invariant_sites_mask(da=med_alltime['veg_idx'], mask=mask_green, max_num_sites=params.max_num_sites)
    
    # basic check
    if not mask_greenest.any():
        print('> No invariant target masks exist or all empty. Aborting.')
        return ds
  
    # standardise every scene to upper percentile within greenest mask
    print('> Standardising raw vegetation indices using targets and basic / max function. Please wait.')
    ds = ds.groupby('time').apply(perform_divide_by_percentile, mask_greenest=mask_greenest, high_percentile=99.5)

    # tell user and return
    print('Basic standardisation for original mannkendall applied successfully. Proceeding.')
    return ds

# get mk info for selected pixel from mk_result var
def get_mk_pixel_info(mk_result, x, y):
    
    # default message
    msg = 'TAU: NA | PVALUE: NA | SENS SLOPE: NA'
    
    # check if x and y exist, else create default string
    if not x or not y:
        return msg

    try:       
        # find closest possible pixel, pass if not found
        vec = mk_result.sel(x=x, y=y, method='nearest', tolerance=30)
        
        # get values if exist, else na
        tau_val = '{0:.3f}'.format(float(vec['tau'])) if np.isfinite(vec['tau']) else 'NA'
        p_val = '{0:.3f}'.format(float(vec['p'])) if np.isfinite(vec['p']) else 'NA'
        z_val = '{0:.3f}'.format(float(vec['z'])) if np.isfinite(vec['z']) else 'NA'
        slope_val = '{0:.3f}'.format(float(vec['slope'])) if np.isfinite(vec['slope']) else 'NA'
        s_val = '{0:.3f}'.format(float(vec['s'])) if np.isfinite(vec['s']) else 'NA'
                
        # populate string
        msg = 'TAU: {0} | PVALUE: {1} | SENS SLOPE: {2}'.format(tau_val, p_val, slope_val)
        return msg

    except:
        print('Warning, something went wrong during get_mk_pixel_info. Reverting to default message')
        return msg
    
# perform change point detection on vector of vegetation values
def do_change_pd(ds, x, y, model, min_size, jump, pelt_penelty, dynp_n_bkps):
    
    # check if x and y exist, else create default string
    if not ds or not x or not y:
        return None, None, None, None
    
    try:       
        # find closest possible pixel, get dates, 
        vec = ds['veg_idx'].sel(x=x, y=y, method='nearest', tolerance=30)
        dates = vec['time']
                        
        # get data as basic array
        vec = vec.data
        
    except:
        return None, None, None, None
                        
    # perform change point detection via pelt
    try:
        result_pelt = rpt.Pelt(model=model, min_size=min_size, jump=jump).fit_predict(signal=vec, pen=pelt_penelty)
    except:
        print('No break points for pelt detected, skipping.')
        result_pelt = None
        
    # perform change point detection via pelt
    try:
        result_dynp = rpt.Dynp(model=model, min_size=min_size, jump=jump).fit_predict(signal=vec, n_bkps=dynp_n_bkps)
    except:
        print('No break points for dynamic programming detected, skipping.')
        result_dynp = None
    
    return vec, dates, result_pelt, result_dynp

# generate summary statistics for change dataset
def calc_summary_stats_change(ds):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) == 0:
        print('No years of data in dataset. Please include more years in your query above.')
        return None

    try:
        # number of dates, years
        total_dates = len(ds['time'])
        total_years = len(ds['time'].groupby('time.year'))

        # get dates that are all nan
        all_nan_dates_veg, all_nan_dates_brt = [], []
        ds['veg_idx'].groupby('time').apply(is_all_nan, all_nan_list=all_nan_dates_veg)
        ds['brt_idx'].groupby('time').apply(is_all_nan, all_nan_list=all_nan_dates_brt)
        
        # get any other dates that have any (not all) nan
        any_nan_dates_veg, any_nan_dates_brt = [], []
        ds['veg_idx'].groupby('time').apply(is_any_nan, any_nan_list=any_nan_dates_veg)
        ds['brt_idx'].groupby('time').apply(is_any_nan, any_nan_list=any_nan_dates_brt)
        
        # get mean of veg, brt
        da_mean_veg = ds['veg_idx'].mean()
        da_mean_brt = ds['brt_idx'].mean()

        # get median of veg, brt
        da_medn_veg = ds['veg_idx'].median()
        da_medn_brt = ds['brt_idx'].median()

        # get std dev of veg, brt
        da_stdv_veg = ds['veg_idx'].std()
        da_stdv_brt = ds['brt_idx'].std()
        
        # relay time summary to user
        print('=== Summary dataset statistics ===')
        print('> Total number of dates in dataset: {0}'.format(total_dates))
        print('> Total number of years in dataset: {0}'.format(total_years))
        
        # relay missing data summary to user
        print('> Completely missing greeness data on dates: {0}'.format(all_nan_dates_veg))
        print('> Completely missing brightness data on dates: {0}'.format(all_nan_dates_brt))
        print('> Partially missing greeness data on dates: {0}'.format(any_nan_dates_veg))
        print('> Partially missing brightness data on dates: {0}'.format(any_nan_dates_brt))
        
        # relay mean stats to user
        msg = '> Mean values for entire dataset : greeness: {0}, brightness: {1}'
        print(msg.format(np.round(da_mean_veg.values, 3), np.round(da_mean_brt.values, 3)))

        # relay median stats to user
        msg = '> Median values for entire dataset: greeness: {0}, brightness: {1}'
        print(msg.format(np.round(da_medn_veg.values, 3), np.round(da_medn_brt.values, 3)))

        # relay std dev stats to user
        msg = '> Standard deviation values for entire wet season: greeness: {0}, brightness: {1}'
        print(msg.format(np.round(da_stdv_veg.values, 3), np.round(da_stdv_brt.values, 3)))

        # close up
        print('=== Summary dataset statistics ===\n')

    except Exception as e:
        print('Error occurred during calc_summary_stats_change of type {0}. Skipping.'.format(e))
        
# calc if any veg, brt are z-score sig for trend dataset, flag, show stats, reset to all nan if outlier
def zscore_flag_show_reset_change(ds, z_val, reset):
    
    # basic check
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None
    elif len(ds['time.year']) == 0:
        print('No years of data in dataset. Please include more years in your query above.')
        return None
    elif z_val == 0:
        print('Z-score critical value is invalid. Please ensure it is not 0.')
        return None
    
    try:
        # tell user
        print('Performing Z-score test to determine data outliers. Please wait.')  
        
        # apply z-score and flag outliers
        ds = ds.groupby('time').apply(do_zscore_and_flag_change, ds=ds, z_val=z_val)
            
        # see if vars in ds and set all flagged to nan
        if all(elem in list(ds.data_vars) for elem in ['veg_outlier', 'brt_outlier']):
            
            # get dates for printing
            veg_outs = [str(e['time'].dt.strftime('%Y-%m-%d').values) for e in ds['time'].where(ds['veg_outlier'], drop=True)]
            brt_outs = [str(e['time'].dt.strftime('%Y-%m-%d').values) for e in ds['time'].where(ds['brt_outlier'], drop=True)]
            
            # show summary
            print('=== Z-Score Results using score of {0} ==='.format(z_val))
            print('> Greeness Outlier Dates: {0}'.format(veg_outs))
            print('> Brightness Outlier Dates: {0}'.format(brt_outs))
            print('=== Z-Score Results using score of {0} ===\n'.format(z_val))
            
            # if reset, set outlier scenes to all nan
            if reset:
                ds = ds.where((ds['veg_outlier'] == False) & (ds['brt_outlier'] == False))

        else:
            print('No vegetation and/or moisture outlier mask. Skipping and proceeding.')
                  
        # drop outliet mask vars
        ds = ds.drop(['veg_outlier', 'brt_outlier'], errors='ignore')
        return ds

    except Exception as e:
        print('Error occurred during zscore_flag_show_reset_change of type {0}. Stopping.'.format(e))
        raise e   
        
# do basic green site standardisation process sites for change
def do_standardisation_indices_change(ds, params):
    
    # basic check
    if not ds:
        print('Error occured after compute, no dataset exists. Aborting.')
        return None
    elif not params.check_change_standardisation_params():
        print('No standardisation parameters set. Please ensure they are set above.')
        return ds
    elif len(ds['time.year']) == 0:
        print('Less than five years of data in dataset. Please include more years in your query above.')
        return ds
            
    # tell user
    print('\nGenerating standardisation targets for vegetation for change vector analysis. Please wait.')
    
    # get median all time for veg
    med_alltime = ds.median('time')
    
    # get n% percentile of veg, create igt masks, give user basic info and display
    print('> Getting highest valued greeness pixels and reducing to max number of sites. Please wait.') 
    mask_green = mask_via_percentile(da=med_alltime['veg_idx'], percentile_val=params.green_percentile)
    mask_greenest = reduce_invariant_sites_mask(da=med_alltime['veg_idx'], mask=mask_green, max_num_sites=params.max_num_sites)
    
    # get n% percentile of brt, create igt masks, give user basic info and display
    print('> Getting highest valued brightness pixels and reducing to max number of sites. Please wait.') 
    mask_bright = mask_via_percentile(da=med_alltime['brt_idx'], percentile_val=params.green_percentile)
    mask_brightest = reduce_invariant_sites_mask(da=med_alltime['brt_idx'], mask=mask_bright, max_num_sites=params.max_num_sites)
    
    # basic check
    if not mask_greenest.any() or not mask_brightest.any():
        print('> No invariant target masks exist or all empty. Aborting.')
        return ds
  
    # standardise every scene to upper percentile within greenest mask
    print('> Standardising raw vegetation and brightness indices using targets and basic / max function. Please wait.')
    ds = ds.groupby('time').apply(perform_divide_by_percentile_change_vector, mask_greenest=mask_greenest, mask_brightest=mask_brightest, 
                                  high_percentile=99.5)

    # tell user and return
    print('Basic standardisation for change vector analysis applied successfully. Proceeding.')
    return ds

# perform change vector analysis on one baseline vs one comparison year
def apply_cva(comp_da, base_da, tmf):
            
    # get differences between bands between start and end dates
    comp_da['brt_diff'] = comp_da['brt_idx'] - base_da['brt_idx']
    comp_da['veg_diff'] = comp_da['veg_idx'] - base_da['veg_idx']
    
    # get magnitude and threshold
    comp_da['magnitude'] = np.sqrt((np.square(comp_da['brt_diff'])) + np.square(comp_da['veg_diff']))
    threshold = tmf * comp_da['magnitude'].where(comp_da['magnitude'] > 0).median(skipna=True).values

    # create mask less than threshold and not nan
    comp_da['sel'] = xr.where((comp_da['magnitude'] > threshold) & (np.isfinite(comp_da['magnitude'])), True, False)

    # if anything came back from mask, calculate angles and magnitudes
    if comp_da['sel'].any():
        comp_da['angout'] = np.arctan2(comp_da['brt_diff'].where(comp_da['sel']), comp_da['veg_diff'].where(comp_da['sel'])) / math.pi * 180
        comp_da['magout'] = comp_da['magnitude'].where(comp_da['sel'])

    # if angle is negative then + 360, followed by masking
    comp_da['angout'] = xr.where(comp_da['angout'] < 0, comp_da['angout'] + 360, comp_da['angout'])
    comp_da['angout'] = comp_da['angout'].where(comp_da['sel'])

    return comp_da

# get cva info for selected pixel from mk_result var
def get_cva_pixel_info(cva_da, x, y):
    
    # default message
    msg = 'ANGLE: NA | MAGNITUDE: NA'
    
    # check if x and y exist, else create default string
    if not x or not y:
        return msg
    
    try:       
        # find closest possible pixel, pass if not found
        vec = cva_da.sel(x=x, y=y, method='nearest', tolerance=30)
        
        # get values if exist, else na
        ang_val = '{0:.3f}'.format(float(vec['angout'])) if np.isfinite(vec['angout']) else 'NA'
        mag_val = '{0:.3f}'.format(float(vec['magout'])) if np.isfinite(vec['magout']) else 'NA'
                        
        # populate string
        msg = 'ANGLE: {0} | MAGNITUDE: {1}'.format(ang_val, mag_val)
        
        return msg

    except:
        print('Warning, something went wrong during get_cva_pixel_info. Reverting to default message')
        return msg
    
# perform change point detection on vector of vegetation and brightness values for cva
def do_cva_change_pd(ds, index, x, y, model, min_size, jump, pelt_penelty, dynp_n_bkps):
    
    # check if x and y exist, else create default string
    if not ds or not x or not y:
        return None, None, None, None
    
    try:       
        # find closest possible pixel, get dates, 
        vec = ds[index].sel(x=x, y=y, method='nearest', tolerance=30)
        dates = vec['time']
                        
        # get data as basic array
        vec = vec.data
        
    except:
        return None, None, None, None
                        
    # perform change point detection via pelt
    try:
        result_pelt = rpt.Pelt(model=model, min_size=min_size, jump=jump).fit_predict(signal=vec, pen=pelt_penelty)
    except:
        result_pelt = None
        
    # perform change point detection via pelt
    try:
        result_dynp = rpt.Dynp(model=model, min_size=min_size, jump=jump).fit_predict(signal=vec, n_bkps=dynp_n_bkps)
    except:
        result_dynp = None
    
    return vec, dates, result_pelt, result_dynp