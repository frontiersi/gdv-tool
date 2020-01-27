# classes for holding user parameters from various interface controls
from pathlib import Path
from gdv_scripts import helpers


# user parameter class for storing user params
class GeneralParams:
    def __init__(self):
        self.project = 'unknown'                         # project name (default project is unknown)
        self.study_center = (-22.925, 119.645)           # default pilbara study area center
        self.out_folder = '/home'                        # working folder location for temp outputs
        
        # try to get home folder
        self.get_home_folder()
    
    def get_home_folder(self):
        try:
            self.out_folder = str(Path.home())
        except:
            self.out_folder = None

# gdv likelihood class for storing gdv params
class LikelihoodParams:
    def __init__(self):
        self.platform = 'ls'                               # platform name (default landsat)
        self.product = 'nbart'                             # product name (default nbart)
        self.query_dates = ('1999-01-01', '2019-12-31')    # query date range
        self.cloud_cover = 0.98                            # maximum cloud cover %
        self.cloud_mask = True                             # mask any existing clouds as nan
        self.crs = 'EPSG:3577'                             # crs of cube dataset. bug removes it sometimes
        self.affine = None                                 # affine coords of cube dataset
        self.gcs_nw = (None, None)                         # nw study area corner lon and lat in albers
        self.gcs_se = (None, None)                         # se study area corner lon and lat in albers
        self.alb_nw = (None, None)                         # nw study area corner x and y in albers
        self.alb_se = (None, None)                         # se study area corner x and y in albers
        self.veg_idx = 'ndvi'                              # vegetation index (default ndvi)
        self.mst_idx = 'ndmi'                              # moisture index (default ndmi)
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.l_value = 0.5                                 # soil l value for savi index
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.plot_raw_idx = False                          # plot raw indices if selected
        self.miss_data_handle = 'interpolate'              # how to handle missing data
        self.fill_pixels = False                           # interpolate pixels within each scene
        self.plot_interp_idx = False                       # plot interpolated indices if selected
        self.green_moist_percentile = 99.5                 # upper percentile for highest veg/mst values during standardisation
        self.slope_steady_percentile = 5                   # lowest percentile steady areas of slope for standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.plot_standard_idx = False                     # plot standardised indices if selected
        self.like_thresh_type = 'stdv'                     # boolean for thresholding likelihood. standard dev, or shape
        self.ground_sites_path = None                      # if thresholding via groundtruth shapefile, provide path
        self.like_thresh_date = 'all'                      # likelihood map to display on map, defaults to all (median all time)
        self.like_thresh_stdv = 1.5                        # standard deviation value to threshold by if no shapefile provided, default to 1.5
        self.map_likelihood = True                         # put results of threshold onto map for inspection
        self.write_data_path = None                        # write data path
        
    def check_general_params(self):
        if not self.platform:
            return False
        if not self.product:
            return False
        if not self.query_dates[0] or not self.query_dates[1]:
            return False
        if self.cloud_cover < 0.0 or self.cloud_cover > 1.0:
            return False
        if not all(self.gcs_nw) or not all(self.gcs_se):
            return False
        if not all(self.alb_nw) or not all(self.alb_se):
            return False
        return True
    
    def check_raw_idx_params(self):
        if not self.veg_idx:
            return False
        if not self.mst_idx:
            return False
        if self.l_value < 0 or self.l_value > 1:
            return False
        return True
    
    def check_interpolate_idx_params(self):
        if self.zscore_value < 1 or self.zscore_value > 3:
            return False
        if not self.miss_data_handle:
            return False
        return True

    def check_prepare_standardisation_params(self):
        if self.green_moist_percentile <= 0 or self.green_moist_percentile >= 100:
            return False
        if self.slope_steady_percentile <= 0 or self.slope_steady_percentile >= 100:
            return False
        if self.max_num_sites <= 0 or self.max_num_sites > 100:
            return False
        return True

    def reset_raw_idx_params(self):
        self.veg_idx = 'ndvi'                              # vegetation index (default ndvi)
        self.mst_idx = 'ndmi'                              # moisture index (default ndmi)
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.l_value = 0.5                                 # soil l value for savi index
        self.plot_raw_idx = False                          # plot raw indices if selected
        
    def reset_interpolate_idx_params(self):
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.miss_data_handle = 'interpolate'              # how to handle missing data 
        self.plot_interp_idx = False                       # plot interpolated indices if selected
        
    def reset_prepare_standardisation_params(self):
        self.green_moist_percentile = 99.5                 # upper percentile for highest veg/mst values during standardisation
        self.slope_steady_percentile = 5                   # lowest percentile steady areas of slope for standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.plot_standard_idx = False                     # plot standardised indices if selected
        
    def reset_likelihood_threshold_params(self):
        self.like_thresh_type = 'stdv'                     # boolean for thresholding likelihood. standard dev, or shape
        self.ground_sites_path = None                      # if thresholding via groundtruth shapefile, provide path
        self.like_thresh_date = 'all'                      # likelihood map to display on map, defaults to all (median all time)
        self.like_thresh_stdv = 1.5                        # standard deviation value to threshold by if no shapefile provided, default to 1.5
        self.map_likelihood = True                         # put results of threshold onto map for inspection
        self.write_data_path = None                        # write data path
        
# trend class for storing trend params     
class TrendParams:
    def __init__(self):
        self.platform = 'ls'                               # platform name (default landsat)
        self.product = 'nbart'                             # product name (default nbart)
        self.query_dates = ('2009-01-01', '2019-12-31')    # query date range
        self.cloud_cover = 0.98                            # maximum cloud cover %
        self.cloud_mask = True                             # mask any existing clouds as nan
        self.crs = 'EPSG:3577'                             # crs of cube dataset. bug removes it sometimes
        self.affine = None                                 # affine coords of cube dataset
        self.gcs_nw = (None, None)                         # nw study area corner lon and lat in albers
        self.gcs_se = (None, None)                         # se study area corner lon and lat in albers
        self.alb_nw = (None, None)                         # nw study area corner x and y in albers
        self.alb_se = (None, None)                         # se study area corner x and y in albers
        self.mk_type = 'seasonal'                          # mannkendall type (original or seasonal)
        self.time_slice = 'qt'                             # time slice, quarters or specific months
        self.mk_sig_only = False                           # return pixels where pvalue is sig < 0.05 only, or not
        self.trend_type = 'all'                            # show all trends, or increasing, or decreasing only
        self.model = 'rbf'                                 # change point detection model type. l1, l2 and rbf offered
        self.min_size = 2                                  # minimum distance between change points
        self.jump = 1                                      # change point indexes multiple of this particular value to consider
        self.dynp_n_bkps = 2                               # number of break points to predict for dynamic programming method
        self.pelt_penelty = 3                              # penelty number for pelt method
        self.veg_idx = 'ndvi'                              # vegetation index (default ndvi)
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.l_value = 0.5                                 # soil l value for savi index
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.miss_data_handle = 'interpolate'              # how to handle missing data
        self.fill_pixels = False                           # interpolate pixels within each scene
        self.green_percentile = 99.5                       # upper percentile for highest veg/mst values during standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.read_data_path = None                         # path to netcdf of likelihood, if exists

    def check_general_params(self):
        if not self.platform:
            return False
        if not self.product:
            return False
        if not self.query_dates[0] or not self.query_dates[1]:
            return False
        if self.cloud_cover < 0.0 or self.cloud_cover > 1.0:
            return False
        if not all(self.gcs_nw) or not all(self.gcs_se):
            return False
        if not all(self.alb_nw) or not all(self.alb_se):
            return False
        return True
    
    def check_trend_idx_params(self):
        if not self.veg_idx:
            return False
        if self.l_value < 0 or self.l_value > 1:
            return False
        return True
    
    def check_trend_interpolate_idx_params(self):
        if self.zscore_value < 1 or self.zscore_value > 3:
            return False
        if not self.miss_data_handle:
            return False
        return True
    
    def check_trend_standardisation_params(self):
        if self.green_percentile <= 0 or self.green_percentile >= 100:
            return False
        if self.max_num_sites <= 0 or self.max_num_sites > 100:
            return False
        return True

    def check_mk_trend_params(self):
        if self.mk_type not in ['original', 'seasonal']:
            return False
        if self.trend_type not in ['all', 'increasing', 'decreasing']:
            return False
        return True        
        
    def reset_trend_params(self):
        self.platform = 'ls'                               # platform name (default landsat)
        self.product = 'nbart'                             # product name (default nbart)
        self.query_dates = ('2009-01-01', '2019-12-31')    # query date range
        self.cloud_cover = 0.98                            # maximum cloud cover %
        self.cloud_mask = True                             # mask any existing clouds as nan
        self.mk_type = 'seasonal'                          # mannkendall type (original or seasonal)
        self.time_slice = 'qt'                             # time slice, quarters or specific months
        self.mk_sig_only = False                           # return pixels where pvalue is sig < 0.05 only, or not
        self.trend_type = 'all'                            # show all trends, or increasing, or decreasing only
        self.veg_idx = 'ndvi'                              # vegetation index (default ndvi)
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.l_value = 0.5                                 # soil l value for savi index
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.miss_data_handle = 'interpolate'              # how to handle missing data
        self.fill_pixels = False                           # interpolate pixels within each scene
        self.green_percentile = 99.5                       # upper percentile for highest veg/mst values during standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.read_data_path = None                         # path to netcdf of likelihood, if exists     

    def reset_cpd_params(self):
        self.model = 'rbf'                                 # change point detection model type. l1, l2 and rbf offered
        self.min_size = 2                                  # minimum distance between change points
        self.jump = 1                                      # change point indexes multiple of this particular value to consider
        self.pelt_penelty = 3                              # penelty number for pelt method
        self.dynp_n_bkps = 2                               # number of break points to predict for dynamic programming method
        
    def correct_sat_launch_date(self):
        if self.platform == 'ls':
            if int(self.query_dates[0].split('-')[0]) < 1986:
                self.query_dates = ('1986-01-01', self.query_dates[1])
                print('> Corrected query date range to 1986-01-01 based on landsat satellite launch date.')
        elif self.platform == 's2':
            if int(self.query_dates[0].split('-')[0]) < 2015:
                self.query_dates = ('2015-01-01', self.query_dates[1])
                print('> Corrected query date range to beginning 2015-01-01 based on sentinel 2 satellite launch date.')

# change class for storing change params
class ChangeParams:
    def __init__(self):
        self.platform = 'ls'                               # platform name (default landsat)
        self.product = 'nbart'                             # product name (default nbart)
        self.base_dates = ('1989-01-01', '2019-12-31')     # baseline date range
        self.comp_dates = ('2009-01-01', '2019-12-31')     # comparison date range
        self.cloud_cover = 0.98                            # maximum cloud cover %
        self.cloud_mask = True                             # mask any existing clouds as nan
        self.crs = 'EPSG:3577'                             # crs of cube dataset. bug removes it sometimes
        self.affine = None                                 # affine coords of cube dataset
        self.gcs_nw = (None, None)                         # nw study area corner lon and lat in albers
        self.gcs_se = (None, None)                         # se study area corner lon and lat in albers
        self.alb_nw = (None, None)                         # nw study area corner x and y in albers
        self.alb_se = (None, None)                         # se study area corner x and y in albers
        self.time_slice = 'q1'                             # time slice, quarters or specific months
        self.veg_idx = 'tcap'                              # vegetation index (default tcap)
        self.change_angles_range = (90, 180)               # cva angle range 90-180 is veg decline
        self.change_threshold = 1.5                        # threshold for cva magnitude
        self.model = 'rbf'                                 # change point detection model type. l1, l2 and rbf offered
        self.min_size = 2                                  # minimum distance between change points
        self.jump = 1                                      # change point indexes multiple of this particular value to consider
        self.dynp_n_bkps = 2                               # number of break points to predict for dynamic programming method
        self.pelt_penelty = 3                              # penelty number for pelt method
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.miss_data_handle = 'interpolate'              # how to handle missing data
        self.fill_pixels = False                           # interpolate pixels within each scene
        self.green_percentile = 99.5                       # upper percentile for highest veg/mst values during standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.read_data_path = None                         # path to netcdf of likelihood, if exists

    def check_general_params(self):
        if not self.platform:
            return False
        if not self.product:
            return False
        if not self.base_dates[0] or not self.base_dates[1]:
            return False
        if not self.comp_dates[0] or not self.comp_dates[1]:
            return False
        if self.cloud_cover < 0.0 or self.cloud_cover > 1.0:
            return False
        if not all(self.gcs_nw) or not all(self.gcs_se):
            return False
        if not all(self.alb_nw) or not all(self.alb_se):
            return False
        return True 
    
    def check_change_idx_params(self):
        if not self.veg_idx:
            return False
        return True
    
    def check_change_interpolate_idx_params(self):
        if self.zscore_value < 0.5 or self.zscore_value > 3:
            return False
        if not self.miss_data_handle:
            return False
        return True
    
    def check_change_standardisation_params(self):
        if self.green_percentile <= 0 or self.green_percentile >= 100:
            return False
        if self.max_num_sites <= 0 or self.max_num_sites > 100:
            return False
        return True
    
    def check_cva_params(self):
        if not self.change_angles_range:
            return False
        if self.change_threshold < 0.5 or self.change_threshold > 3:
            return False
        return True   
    
    def reset_change_params(self):
        self.platform = 'ls'                               # platform name (default landsat)
        self.product = 'nbart'                             # product name (default nbart)
        self.base_dates = ('1989-01-01', '2019-12-31')     # baseline date range
        self.comp_dates = ('2009-01-01', '2019-12-31')     # comparison date range
        self.cloud_cover = 0.98                            # maximum cloud cover %
        self.cloud_mask = True                             # mask any existing clouds as nan
        self.time_slice = 'q1'                             # time slice, quarters or specific months
        self.veg_idx = 'tcap'                              # vegetation index (default tcap)
        self.change_angles_range = (90, 180)               # cva angle range 90-180 is veg decline
        self.change_threshold = 1.5                        # threshold for cva magnitude
        self.rescale = True                                # rescale indices if -1 to 1 (go 0 to 2)
        self.zscore = True                                 # remove outlier scenes via zscore
        self.zscore_value = 1.96                           # z-score value (default 1.96, which is p 0.05)
        self.miss_data_handle = 'interpolate'              # how to handle missing data
        self.fill_pixels = False                           # interpolate pixels within each scene
        self.green_percentile = 99.5                       # upper percentile for highest veg/mst values during standardisation
        self.max_num_sites = 50                            # maximum number of invariant sites
        self.read_data_path = None                         # path to netcdf of likelihood, if exists
        
    def reset_cpd_params(self):
        self.model = 'rbf'                                 # change point detection model type. l1, l2 and rbf offered
        self.min_size = 2                                  # minimum distance between change points
        self.jump = 1                                      # change point indexes multiple of this particular value to consider
        self.pelt_penelty = 3                              # penelty number for pelt method
        self.dynp_n_bkps = 2                               # number of break points to predict for dynamic programming method
        
    def correct_sat_launch_date(self):
        if self.platform == 'ls':
            if int(self.base_dates[0].split('-')[0]) < 1986:
                self.base_dates = ('1986-01-01', self.base_dates[1])
                print('> Corrected base date range to 1986-01-01 based on landsat satellite launch date.')
            if int(self.comp_dates[0].split('-')[0]) < 1986:
                self.comp_dates = ('1986-01-01', self.comp_dates[1])
                print('> Corrected comparison date range to 1986-01-01 based on landsat satellite launch date.')
        elif self.platform == 's2':
            if int(self.base_dates[0].split('-')[0]) < 2015:
                self.base_dates = ('2015-01-01', self.base_dates[1])
                print('> Corrected base date range to beginning 2015-01-01 based on sentinel 2 satellite launch date.')
            if int(self.comp_dates[0].split('-')[0]) < 2015:
                self.comp_dates = ('2015-01-01', self.comp_dates[1])
                print('> Corrected comparison date range to beginning 2015-01-01 based on sentinel 2 satellite launch date.')