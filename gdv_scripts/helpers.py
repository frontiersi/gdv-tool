# general helper functions for various repeated processes
from datetime import datetime
from pyproj import Proj as projection
from pyproj import transform
from tqdm.notebook import tqdm
import xarray as xr

# add a progress bar to xarray for apply function (based on https://github.com/tqdm/tqdm/issues/653)
def enable_xr_progressbar():
    def inner_generator(df_function='apply'):
        
        def inner(df, func, *args, **kwargs):
            t = tqdm(total = len(df))
            
            def wrapper(*args, **kwargs):
                t.update( n=1 if not t.total or t.n < t.total else 0 )
                return func( *args, **kwargs )

            result = getattr(df, df_function)(wrapper, **kwargs)

            t.close()
            
            return result
        
        return inner

    xr.core.groupby.DataArrayGroupBy.progress_apply = inner_generator()
    xr.core.groupby.DatasetGroupBy.progress_apply = inner_generator()

# convert tuple of years into cube-friendly query dates
def prepare_dates(date_tuple):
    if date_tuple:
        try:
            start = '{0}-{1}'.format(date_tuple[0], '01-01')
            end = '{0}-{1}'.format(date_tuple[1], '12-31')
            
            return (start, end)
        except Exception as e:
            print('Error occurred during prepare_dates of type {0}. Stopping.'.format(e))
            raise e
    else:
        return ('1980-01-01', '2029-01-01')

# check if selected range of years is valid
def valid_year_range(date_tuple):
    if date_tuple:
        try:
            start = datetime.strptime(date_tuple[0], '%Y-%m-%d').year
            end = datetime.strptime(date_tuple[1], '%Y-%m-%d').year
            diff = (end + 1) - start
            
            return diff
        except Exception as e:
            print('Error occurred during valid_year_range of type {0}. Stopping.'.format(e))
            raise e
            
# get earliest and latest dates from two tuples (for change vector)
def get_earliest_and_latest_date_range(base_dates, compare_dates):
    try:       
        #base dates
        base_s_date = datetime.strptime(base_dates[0], '%Y-%m-%d')
        base_e_date = datetime.strptime(base_dates[1], '%Y-%m-%d')
        
        # comparison dates
        comp_s_date = datetime.strptime(compare_dates[0], '%Y-%m-%d')
        comp_e_date = datetime.strptime(compare_dates[1], '%Y-%m-%d')

        # get earliest date or bust
        if base_s_date.year <= comp_s_date.year:
            earliest_date = str(base_s_date)
            print('Earliest date is {0}. Proceeding.'.format(earliest_date))
        else:
            print('Comparison date earlier than baseline. Please change.')
            return None
        
        # get latest date or bust
        if base_e_date.year >= comp_e_date.year:
            latest_date = str(base_e_date)
            print('Latest date is {0}. Proceeding.'.format(latest_date))
        else:
            latest_date = str(comp_e_date)
            print('Latest date is {0}. Proceeding.'.format(latest_date))
        
        return earliest_date, latest_date
        
    except:
        print('Something went wrong. Check your dates.')

# convert cloud cover percentage 0-100% to inverted float e.g. 5 is 0.95 in cube query
def prepare_cloud_percentage(raw_cover_val):
    if (raw_cover_val >= 0) and (raw_cover_val <= 1):
        try:
            clean_cover_val = round(1.0 - raw_cover_val, 2)

            return clean_cover_val
        except Exception as e:
            print('Error occurred during prepare_cloud_percentage of type {0}. Stopping.'.format(e))
            raise e
    else:
        return 0.98

# reproject x, y coordinate from wgs84 to albers or vice versa 'epsg:4326' or 'epsg:3577'
def transform_point(lon, lat, in_epsg, out_epsg):
    try:
        in_prj = projection(init=in_epsg)
        out_prj = projection(init=out_epsg)
        x, y = transform(in_prj, out_prj, lon, lat)
            
        return x, y
    except Exception as e:
        print('Error occurred during transform_point_wgs_to_alb of type {0}. Stopping.'.format(e))
        raise e

# find platform via pixel resolution
def identify_platform(da):
    try:
        if abs(da.x[0] - da.x[1]) > 20:
            return 'ls'
        else:
            return 's2'
    except:
        return ''