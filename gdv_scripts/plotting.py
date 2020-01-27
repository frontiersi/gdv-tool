# scripts for plotting various figures
import warnings
import numpy as np
import pandas as pd
import ruptures as rpt
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from itertools import cycle
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# use ggplot style
plt.style.use('ggplot')

## HELPER FUNCTIONS ##
def normalise_values(da, stretch_min, stretch_max):
    try:
        vmin, vmax = da.quantile([stretch_min, stretch_max]).values
        norm = clrs.Normalize(vmin=vmin, vmax=vmax, clip=True)
        
        return norm
    
    except Exception as e:
        print('Error occurred during normalise_values of type {0}. Stopping.'.format(e))
        raise e   
    
    
## MAIN FUNCTIONS ##
# plot summary stats for raw vegetation and moisture index data
def plot_raw_idx_stats(ds, fignum):

    if ds and fignum:
        # tell user
        print('Generating summary statistics and graphs, please hold.')
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            # basic checks
            if not ds:
                print('No cube dataset provided. Please ensure all cells above have been run.')
                return None, None

            # get annual means of veg and mst indices for wet and dry seasons
            ds_wet_means = ds.where(ds['time.month'] < 6).groupby('time.year').mean(...)
            ds_dry_means = ds.where(ds['time.month'] > 6).groupby('time.year').mean(...)
            
            # generate fig and axes
            width = 0.35
            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)
            
            # get years and list of positions
            x = ds.groupby('time.year').mean()['year'].values
            x_pos = np.arange(len(ds_wet_means.year.values))

            # get vegetation data
            wet_means = ds_wet_means['veg_idx'].data
            dry_means = ds_dry_means['veg_idx'].data

            # vegetation plot
            axes[0].bar(x_pos, wet_means, width, label='Wet', color='forestgreen', align='center', alpha=0.75)
            axes[0].bar(x_pos + width, dry_means, width, label='Dry', color='limegreen', align='center', alpha=0.75)
            axes[0].axhline(np.nanmean(wet_means), label='Wet Mean', linestyle='--', color='red')
            axes[0].axhline(np.nanmean(dry_means), label='Dry Mean', linestyle='--', color='orange')
            axes[0].set_title('Annual Vegetation Means for Wet/Dry Seasons')
            axes[0].set_xlabel('Year')
            axes[0].set_ylabel('Vegetation Index')
            axes[0].set_ylim(np.nanmin(dry_means) - (np.nanmin(dry_means) * 0.1), 
                             np.nanmax(wet_means) + (np.nanmax(wet_means) * 0.1))
            axes[0].set_xticks(x_pos)
            axes[0].set_xticklabels(x, rotation=90)
            axes[0].legend(loc='best')

            # get moisture data
            wet_means = ds_wet_means['mst_idx'].data
            dry_means = ds_dry_means['mst_idx'].data

            # moisture plot
            axes[1].bar(x_pos, wet_means, width, label='Wet', color='royalblue', align='center', alpha=0.75)
            axes[1].bar(x_pos + width, dry_means, width, label='Dry', color='deepskyblue', align='center', alpha=0.75)
            axes[1].axhline(np.nanmean(wet_means), label='Wet Mean', linestyle='--', color='red')
            axes[1].axhline(np.nanmean(dry_means), label='Dry Mean', linestyle='--', color='orange')
            axes[1].set_title('Annual Moisture Means for Wet/Dry Seasons')
            axes[1].set_xlabel('Year')
            axes[1].set_ylabel('Moisture Index')
            axes[1].set_ylim(np.nanmin(dry_means) - (np.nanmin(dry_means) * 0.1), 
                             np.nanmax(wet_means) + (np.nanmax(wet_means) * 0.1))
            axes[1].set_xticks(x_pos)
            axes[1].set_xticklabels(x, rotation=90)
            axes[1].legend(loc='best')

            return fig, axes
    else:
        print('Did not provide a valid dataset.')

# plot raw vegetation and moisture data
def plot_raw_idx(ds, stretch_min, stretch_max):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # num fo columns
    ncols = 4
    
    # get aspect and size
    try:
        aspect =  ds.dims['x'] / ds.dims['y']
        size = 3 / aspect
    except:
        aspect = 1
        size = 3

    # create vegetation indices figure
    norm = normalise_values(da=ds['veg_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    fig1 = ds['veg_idx'].plot.imshow(col='time', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[], 
                                     norm=norm, cmap='RdYlGn', add_labels=False, add_colorbar=False)
    fig1.set_xlabels('')
    fig1.set_ylabels('')
    fig1.fig.suptitle(t='Annual Median Vegetation Index for Wet (DJF) and Dry (SON) Seasons', 
                     fontsize='x-large', fontweight='bold', y=1.015)

    # create moisture indices figure
    norm = normalise_values(da=ds['mst_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    fig2 = ds['mst_idx'].plot.imshow(col='time', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[], 
                                     norm=norm, cmap='RdYlBu', add_labels=False, add_colorbar=False)
    fig2.set_xlabels('')
    fig2.set_ylabels('')
    fig2.fig.suptitle(t='Annual Median Moisture Index for Wet (DJF) and Dry (SON) Seasons', 
                      fontsize='x-large', fontweight='bold', y=1.015)
    
    return fig1, fig2

# plot vegetation and moisture medians all time
def plot_median_alltime_idx_dry(ds, stretch_min, stretch_max, fignum):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # generate fig and axes
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)

    # plot vegetation median all time
    norm = normalise_values(da=ds['veg_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    im_veg = ds['veg_idx'].plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, norm=norm, cmap='RdYlGn', add_colorbar=False)
    axes[0].set_title('All-time Vegetation Median for Dry Season', fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal')

    # plot moisture median all time
    norm = normalise_values(da=ds['mst_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    im_mst = ds['mst_idx'].plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, norm=norm, cmap='RdYlBu', add_colorbar=False)
    axes[1].set_title('All-time Moisture Median for Dry Season', fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal')
    
    return fig, axes

# plot masks of greennest, moistest vegetation and moisture medians all time
def plot_mask_alltime_idx_dry(mask_green, mask_moist, fignum):
    
    # basic checks
    if not mask_green.any() or not mask_moist.any():
        print('No data arrays provided. Please ensure all cells above have been run.')
        return None, None
    
    # for some reason matplot renders binary data weird... convert 0 to nan temp here
    mask_green = mask_green.where(mask_green == 1)
    mask_moist = mask_moist.where(mask_moist == 1)

    # generate fig and axes and cbar classes
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)
    class_bins = np.array([1])
    cmap_veg = clrs.ListedColormap(['seagreen'])
    cmap_mst = clrs.ListedColormap(['royalblue'])

    # plot vegetation median all time mask
    im_veg = mask_green.plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_veg, add_colorbar=False)
    axes[0].set_title('All-time Greenest Pixels (Dry)', fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal', ticks=class_bins)

    # plot moisture median all time mask
    im_mst = mask_moist.plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_mst, add_colorbar=False)
    axes[1].set_title('All-time Moistest Pixels Median (Dry)', fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal', ticks=class_bins)

    return fig, axes

# plot vegetation and moisture orthogonal polynomial coeff slopes
def plot_opc_slope_idx_dry(slope_veg, slope_mst, fignum):
    
    # basic checks
    if not slope_veg.any() or not slope_mst.any():
        print('No data arrays provided. Please ensure all cells above have been run.')
        return None, None
    
    # generate fig and axes
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)

    # plot vegetation median all time
    im_veg = slope_veg.plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=True, cmap='twilight_shifted_r', add_colorbar=False)
    axes[0].set_title('All-time Vegetation Trend (Dry)', fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal')

    # plot moisture median all time
    im_mst = slope_mst.plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=True, cmap='twilight_shifted_r', add_colorbar=False)
    axes[1].set_title('All-time Moisture Trend (Dry)', fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal')
    
    return fig, axes

# plot masks of steadiest vegetation and moisture opc slopes
def plot_mask_steady_slope_dry(mask_steady_veg, mask_steady_mst, fignum):
    
    # basic checks
    if not mask_steady_veg.any() or not mask_steady_mst.any():
        print('No data arrays provided. Please ensure all cells above have been run.')
        return None, None

    # generate fig and axes and cbar classes and properties
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)
    class_bins = np.array([0, 1])
    cmap = clrs.ListedColormap(['whitesmoke', 'black'])

    # plot steady veg slope median all time mask
    im_veg = mask_steady_veg.plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap, add_colorbar=False)
    axes[0].set_title('Steadiest Vegetation Pixels (Dry)', fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal', ticks=class_bins)

    # plot moisture median all time mask
    im_mst = mask_steady_mst.plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap, add_colorbar=False)
    axes[1].set_title('Steadiest Moisture Pixels (Dry)', fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal', ticks=class_bins)

    return fig, axes

# plot masks for vegetation and moisture invariant targets 
def plot_invariant_targets_dry(invariant_targets_veg, invariant_targets_mst, fignum):
    
    # basic checks
    if not invariant_targets_veg.any() or not invariant_targets_mst.any():
        print('No data arrays provided. Please ensure all cells above have been run.')
        return None, None
    
    # get number of sites for title
    try:
        num_sites_veg, num_sites_mst = int(invariant_targets_veg.sum()), int(invariant_targets_mst.sum())
    except:
        num_sites_veg, num_sites_mst = 0, 0
        
    # for some reason matplot renders binary data weird... convert 0 to nan temp here
    invariant_targets_veg = invariant_targets_veg.where(invariant_targets_veg == 1)
    invariant_targets_mst = invariant_targets_mst.where(invariant_targets_mst == 1)

    # generate fig and axes and cbar classes and properties
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)
    class_bins = np.array([1])
    cmap_veg = clrs.ListedColormap(['seagreen'])
    cmap_mst = clrs.ListedColormap(['royalblue'])
    
    # plot steady veg slope median all time mask
    im_veg = invariant_targets_veg.plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_veg, add_colorbar=False)
    axes[0].set_title('All ({0}) Invariant Vegetation Targets (Dry)'.format(num_sites_veg), fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal', ticks=class_bins)

    # plot moisture median all time mask
    im_mst = invariant_targets_mst.plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_mst, add_colorbar=False)
    axes[1].set_title('All ({0}) Invariant Moisture Targets (Dry)'.format(num_sites_mst), fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal', ticks=class_bins)

    return fig, axes

# plot reduced vegetation and moisture invariant targets
def plot_reduced_invariant_targets_dry(invariant_targets_veg, invariant_targets_mst, fignum):
    
    # basic checks
    if not invariant_targets_veg.any() or not invariant_targets_mst.any():
        print('No data arrays provided. Please ensure all cells above have been run.')
        return None, None
    
    # for some reason matplot renders binary data weird... convert 0 to nan temp here
    invariant_targets_veg = invariant_targets_veg.where(invariant_targets_veg == 1)
    invariant_targets_mst = invariant_targets_mst.where(invariant_targets_mst == 1)
    
    # get number of sites for title
    try:
        num_sites_veg, num_sites_mst = int(invariant_targets_veg.sum()), int(invariant_targets_mst.sum())
    except:
        num_sites_veg, num_sites_mst = 0, 0

    # generate fig and axes and cbar classes and properties
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)
    class_bins = np.array([1])
    cmap_veg = clrs.ListedColormap(['seagreen'])
    cmap_mst = clrs.ListedColormap(['royalblue'])
    
    # plot steady veg slope median all time mask
    im_veg = invariant_targets_veg.plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_veg, add_colorbar=False)
    axes[0].set_title('Reduced ({0}) Vegetation Targets (Dry)'.format(num_sites_veg), fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal', ticks=class_bins)

    # plot moisture median all time mask
    im_mst = invariant_targets_mst.plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, cmap=cmap_mst, add_colorbar=False)
    axes[1].set_title('Reduced ({0}) Moisture Targets (Dry)'.format(num_sites_mst), fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal', ticks=class_bins)

    return fig, axes

# plot standardised vegetation and moisture data
def plot_standard_idx(ds, stretch_min, stretch_max):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # num fo columns
    ncols = 4
    
    # get aspect and size
    try:
        aspect =  ds.dims['x'] / ds.dims['y']
        size = 3 / aspect
    except:
        aspect = 1
        size = 3

    # create vegetation indices figure
    norm = normalise_values(da=ds['veg_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    fig1 = ds['veg_idx'].plot.imshow(col='time', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[],
                                     norm=norm, cmap='RdYlGn', add_labels=False, add_colorbar=False)
    fig1.set_xlabels('')
    fig1.set_ylabels('')
    fig1.fig.suptitle(t='Standardised Vegetation Index for Wet (DJF) and Dry (SON) Seasons', 
                     fontsize='x-large', fontweight='bold', y=1.015)

    # create moisture indices figure
    norm = normalise_values(da=ds['mst_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    fig2 = ds['mst_idx'].plot.imshow(col='time', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[],
                                     norm=norm, cmap='RdYlBu', add_labels=False, add_colorbar=False)
    fig2.set_xlabels('')
    fig2.set_ylabels('')
    fig2.fig.suptitle(t='Standardised Moisture Index for Wet (DJF) and Dry (SON) Seasons', 
                      fontsize='x-large', fontweight='bold', y=1.015)
    
    return fig1, fig2

# plot raw vegetation and moisture stability data
def plot_raw_stability_idx(ds):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # num fo columns
    ncols = 5
    
    # get aspect and size
    try:
        aspect =  ds.dims['x'] / ds.dims['y']
        size = 2.4 / aspect
    except:
        aspect = 1
        size = 2.4      

    # create vegetation indices figure
    fig1 = ds['veg_stb'].plot.imshow(col='time', col_wrap=5, size=size, aspect=aspect, xticks=[], yticks=[],
                                     robust=False, cmap='seismic', add_labels=False, add_colorbar=False)
    fig1.set_xlabels('')
    fig1.set_ylabels('')
    fig1.fig.suptitle(t='Vegetation Stability between Wet (DJF) and Dry (SON) Seasons', 
                     fontsize='x-large', fontweight='bold', y=1.015)

    # create moisture indices figure
    fig2 = ds['mst_stb'].plot.imshow(col='time', col_wrap=5, size=size, aspect=aspect, xticks=[], yticks=[],
                                     robust=False, cmap='seismic', add_labels=False, add_colorbar=False)
    fig2.set_xlabels('')
    fig2.set_ylabels('')
    fig2.fig.suptitle(t='Moisture Stability between Wet (DJF) and Dry (SON) Seasons', 
                      fontsize='x-large', fontweight='bold', y=1.015)
    
    return fig1, fig2

# plot raw vegetation and moisture stability data
def plot_stand_stable_idx(ds):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # num fo columns
    ncols = 5
    
    # get aspect and size
    try:
        aspect =  ds.dims['x'] / ds.dims['y']
        size = 2.4 / aspect
    except:
        aspect = 1
        size = 2.4 

    # create vegetation indices figure
    fig1 = ds['veg_stb'].plot.imshow(col='time', col_wrap=5, size=size, aspect=aspect, xticks=[], yticks=[],
                                     robust=False, cmap='Greens', add_labels=False, add_colorbar=False)
    fig1.set_xlabels('')
    fig1.set_ylabels('')
    fig1.fig.suptitle(t='Standardised Vegetation Stability per Year', 
                     fontsize='x-large', fontweight='bold', y=1.015)

    # create moisture indices figure
    fig2 = ds['mst_stb'].plot.imshow(col='time', col_wrap=5, size=size, aspect=aspect, xticks=[], yticks=[],
                                     robust=False, cmap='Blues', add_labels=False, add_colorbar=False)
    fig2.set_xlabels('')
    fig2.set_ylabels('')
    fig2.fig.suptitle(t='Standardised Moisture Stability per Year', 
                      fontsize='x-large', fontweight='bold', y=1.015)
    
    return fig1, fig2

# plot likelihood data
def plot_likelihood(ds, stretch_min, stretch_max):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # num fo columns
    ncols = 4
    
    # get aspect and size
    try:
        aspect =  ds.dims['x'] / ds.dims['y']
        size = 3 / aspect
    except:
        aspect = 1
        size = 3

    # create vegetation indices figure
    norm = normalise_values(da=ds['likelihood'], stretch_min=stretch_min, stretch_max=stretch_max)
    fig1 = ds['likelihood'].plot.imshow(col='time', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[], 
                                        norm=norm, cmap='jet_r', add_labels=False, add_colorbar=False)
    fig1.set_xlabels('')
    fig1.set_ylabels('')
    fig1.fig.suptitle(t='Raw Groundwater Dependent Vegetation per Year', 
                     fontsize='x-large', fontweight='bold', y=1.015)

    
    return fig1

# plot bush fire over likelihood
def plot_nbr_idx(like_years, nbr_years):
    
    # basic checks
    if not like_years.any() or not nbr_years.any():
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None

    # num fo columns
    ncols = 4

    # get aspect and size
    try:
        aspect =  like_years.dims['x'] / like_years.dims['y']
        size = 3 / aspect
    except:
        aspect = 1
        size = 3

    # create vegetation indices facetgrid 
    fig = like_years.plot.imshow(col='year', col_wrap=4, size=size, aspect=aspect, xticks=[], yticks=[], cmap='jet_r', 
                                     add_labels=False, add_colorbar=False)

    # iterate axes and add burn and likelihood above it
    for ax in fig.axes.flat:
        title = ax.get_title().split(' ')[-1]
        if title:
            year = int(title)
            scene_like = like_years.sel(year=year)
            ax.imshow(scene_like, alpha=0.25, cmap='jet_r', zorder=0)
            ax.set_title('{0} - {1}'.format(year - 1, year))

            try:
                scene_nbr = nbr_years.sel(year=year)
                ax.imshow(scene_nbr, alpha=0.90, cmap='binary_r', zorder=1)
            except:
                pass

    # clean up labels and sup title
    fig.set_xlabels('')
    fig.set_ylabels('')
    fig.fig.suptitle(t='Potential Burn Areas', fontsize='x-large', fontweight='bold', y=1.015)

    return fig

# plot shapefile points and roc analysis result
def plot_shape_and_roc_stats(da, shp, fpr, tpr, auc, fignum):

    # basic checks
    if not da.any():
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
        
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
            
        # tell user
        print('Generating shapefile point and ROC curve, please hold.')
            
        # generate fig and axes
        width = 0.35
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True, gridspec_kw={'width_ratios': [4, 2]})
        
        # add basic vege map and shapefile points on figure
        class_bins = np.array([0, 1])
        cmap = clrs.ListedColormap(['red', 'blue'])
    
        # plot general veg basemap and shapefile points on top
        cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
        im_veg = da.plot(ax=axes[0], xticks=[], yticks=[], cmap='Greys_r', add_labels=False, add_colorbar=False)
        im_shp = shp.plot(ax=axes[0], column='GDV_ACT', marker='o', markersize=5, cmap=cmap, legend=True, cax=cax, legend_kwds={'orientation': 'horizontal', 'ticks': class_bins})
        axes[0].set_title('Groundtruth Sites (Presence/Absence)', fontsize=12)   
        axes[0].set_aspect('equal')
        axes[0].set_xlabel('')
        axes[0].set_ylabel('')
        
        # plot roc results
        axes[1].set_title('ROC Curve', fontsize=12)
        axes[1].plot(fpr, tpr, 'b', label='AUC (%0.2f)' % auc)        
        axes[1].plot([0, 1], [0, 1], 'r--', label='Random')
        axes[1].legend(loc='lower right')
        axes[1].set_xlim([0, 1])
        axes[1].set_ylim([0, 1])
        axes[1].set_xlabel('False Positive Rate')
        axes[1].set_ylabel('True Positive Rate')
     
        return fig, axes

# plot summary stats for trend raw vegetation data    
def plot_raw_trend_stats(ds, fignum):

    if ds and fignum:
        # tell user
        print('Generating summary statistics and graphs, please hold.')
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            # basic checks
            if not ds:
                print('No cube dataset provided. Please ensure all cells above have been run.')
                return None, None

            # get annual means of veg and mst indices for wet and dry seasons
            ds_means = ds.groupby('time').mean(...)

            # generate fig and axes
            width = 0.35
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12.5, 4), num=fignum, constrained_layout=True)

            # get years and list of positions
            x = ds['time'].dt.strftime('%Y-%m').values
            x_pos = np.arange(len(x))

            # get vegetation data
            means = ds_means['veg_idx'].data

            # vegetation plot
            axes.bar(x_pos, means, width, color='forestgreen', align='center', alpha=0.75)
            axes.set_title('Vegetation Means for Seasons/Months')
            axes.set_xlabel('Date')
            axes.set_ylabel('Vegetation Index')
            axes.set_ylim(np.nanmin(means) - (np.nanmin(means) * 0.05), 
                          np.nanmax(means) + (np.nanmax(means) * 0.05))
            axes.set_xticks(x_pos)
            axes.set_xticklabels(x, rotation=90)
            
            # correct ticks for large amounts of data
            if len(x_pos) > 50:
                every_nth = 4
                for n, label in enumerate(axes.xaxis.get_ticklabels()):
                    if n % every_nth != 0:
                        label.set_visible(False)

            return fig, axes
    else:
        print('Did not provide a valid dataset.')

# plot change point detection results for pelt cbd
def plot_cpd(signal, dates, result, cbd_type, fignum):
     
    if signal.ndim == 1:
        signal = signal.reshape(-1, 1)
    n_samples, n_features = signal.shape
    
    # set options
    figtitle = 'Vegetation Health & Change Point Results ({0})'.format(cbd_type)
    alpha = 0.25                     # transparency of the colored background
    color = "k"                      # color of the lines indicating the computed_chg_pts
    linewidth = 3                    # linewidth of the lines indicating the computed_chg_pts
    linestyle = "--"                 # linestyle of the lines indicating the computed_chg_pts

    # clear existing figure if it exists
    fig = plt.figure(figsize=(12.5, 4), num=fignum)
    plt.figure(fignum).clear()

    # create proper subplot figure
    fig, axarr = plt.subplots(n_features, figsize=(12.5, 4), sharex=True, num=fignum)

    if n_features == 1:
        axarr = [axarr]

    for axe, sig in zip(axarr, signal.T):
        color_cycle = cycle(["#4286f4", "#f44174"])
        
        # plot signal, color each (true) regime
        axe.plot(range(n_samples), sig)
        bkps = [0] + sorted(result)       

        for (start, end), col in zip(rpt.utils.pairwise(bkps), color_cycle):
            axe.axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)

        # vertical lines to mark the result
        if result is not None:
            for bkp in result:
                if bkp != 0 and bkp < n_samples:
                    axe.axvline(x=bkp - 0.5, color=color, linewidth=linewidth, linestyle=linestyle)

    # change axis labels
    if len(dates) > 0:
        date_list = [pd.to_datetime(str(d.data)).strftime('%Y-%m') for d in dates]
        plt.xticks(list(range(0, n_samples)), date_list, rotation='vertical')

    # set title, subtitle, axis labels, etc
    axarr[0].set_title(figtitle, loc='left')
    axarr[0].set_ylabel('Vegetation Health')
    axarr[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axarr[0].axhline(np.nanmedian(signal), label='Mean', linestyle='--', color='orange')

    # final clean up
    fig.tight_layout()
    
    return fig, axarr

# plot summary stats for change raw vegetation data    
def plot_raw_change_stats(ds, fignum):

    if ds and fignum:
        # tell user
        print('Generating summary statistics and graphs, please hold.')
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            # basic checks
            if not ds:
                print('No cube dataset provided. Please ensure all cells above have been run.')
                return None, None

            # get annual means of veg and brt indices
            ds_means = ds.groupby('time').mean(...)

            # generate fig and axes
            width = 0.35
            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)

            # get years and list of positions
            x = ds['time'].dt.strftime('%Y-%m').values
            x_pos = np.arange(len(x))

            # get vegetation data
            veg_means = ds_means['veg_idx'].data

            # vegetation plot
            axes[0].bar(x_pos, veg_means, width, color='forestgreen', align='center', alpha=0.75)
            axes[0].set_title('Greeness Means for Seasons/Months')
            axes[0].set_xlabel('Date')
            axes[0].set_ylabel('Greeness')
            axes[0].set_ylim(np.nanmin(veg_means) - (np.nanmin(veg_means) * 0.05), 
                          np.nanmax(veg_means) + (np.nanmax(veg_means) * 0.05))
            axes[0].set_xticks(x_pos)
            axes[0].set_xticklabels(x, rotation=90)
            
            # correct ticks for large amounts of data
            if len(x_pos) > 50:
                every_nth = 4
                for n, label in enumerate(axes[0].xaxis.get_ticklabels()):
                    if n % every_nth != 0:
                        label.set_visible(False)
            
            # get brightness data
            brt_means = ds_means['brt_idx'].data

            # vegetation plot
            axes[1].bar(x_pos, brt_means, width, color='firebrick', align='center', alpha=0.75)
            axes[1].set_title('Brightness (soil) Means for Seasons/Months')
            axes[1].set_xlabel('Date')
            axes[1].set_ylabel('Brightness')
            axes[1].set_ylim(np.nanmin(brt_means) - (np.nanmin(brt_means) * 0.05), 
                          np.nanmax(brt_means) + (np.nanmax(brt_means) * 0.05))
            axes[1].set_xticks(x_pos)
            axes[1].set_xticklabels(x, rotation=90)
            
            # correct ticks for large amounts of data
            if len(x_pos) > 50:
                every_nth = 4
                for n, label in enumerate(axes[1].xaxis.get_ticklabels()):
                    if n % every_nth != 0:
                        label.set_visible(False)

            return fig, axes
    else:
        print('Did not provide a valid dataset.')

# plot greeness and brightness medians all time for change detection
def plot_median_alltime_idx_change(ds, stretch_min, stretch_max, fignum):
    
    # basic checks
    if not ds:
        print('No cube dataset provided. Please ensure all cells above have been run.')
        return None, None
    
    # generate fig and axes
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 4), num=fignum, constrained_layout=True)

    # plot vegetation median all time
    norm = normalise_values(da=ds['veg_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    im_veg = ds['veg_idx'].plot(ax=axes[0], xticks=[], yticks=[], add_labels=False, robust=False, norm=norm, cmap='RdYlGn', add_colorbar=False)
    axes[0].set_title('Baseline Greeness Median', fontsize=12)   
    axes[0].set_aspect('equal')
    cax = make_axes_locatable(axes[0]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_veg, cax=cax, orientation='horizontal')

    # plot brightness median all time
    norm = normalise_values(da=ds['brt_idx'], stretch_min=stretch_min, stretch_max=stretch_max)
    im_mst = ds['brt_idx'].plot(ax=axes[1], xticks=[], yticks=[], add_labels=False, robust=False, norm=norm, cmap='BrBG_r', add_colorbar=False)
    axes[1].set_title('Baseline Brightness Median', fontsize=12)
    axes[1].set_aspect('equal')
    cax = make_axes_locatable(axes[1]).append_axes('bottom', size='5%', pad=0.1)
    plt.colorbar(im_mst, cax=cax, orientation='horizontal')
    
    return fig, axes

# plot change point detection results for pelt cbd for cva
def plot_cva_cpd(signal_veg, signal_bright, dates, result_veg, result_brt, cbd_type, fignum):
        
    if signal_veg.ndim == 1:
        signal_veg = signal_veg.reshape(-1, 1)
        signal_bright = signal_bright.reshape(-1, 1)
        
    n_samples, n_features = signal_veg.shape
    n_features = 2
            
    # set options
    figtitle = 'Greenness & Brightness with Change Point Results ({0})'.format(cbd_type)
    alpha = 0.25                     # transparency of the colored background
    color = "k"                      # color of the lines indicating the computed_chg_pts
    linewidth = 3                    # linewidth of the lines indicating the computed_chg_pts
    linestyle = "--"                 # linestyle of the lines indicating the computed_chg_pts

    # clear existing figure if it exists
    fig = plt.figure(figsize=(12.5, 5), num=fignum)
    plt.figure(fignum).clear()

    # create proper subplot figure
    fig, axarr = plt.subplots(n_features, figsize=(12.5, 5), sharex=True, num=fignum)

    # set color cycle for break fills
    color_cycle = cycle(["#4286f4", "#f44174"])
    
    # greeness plot
    # plot signal, color each (true) regime
    axarr[0].plot(range(n_samples), signal_veg, color='g')
    bkps = [0] + sorted(result_veg)
    
    # set start and end areas for background color of blue or red
    for (start, end), col in zip(rpt.utils.pairwise(bkps), color_cycle):
        axarr[0].axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
        
    # vertical lines to mark the result
    if result_veg is not None:
        for bkp in result_veg:
            if bkp != 0 and bkp < n_samples:
                axarr[0].axvline(x=bkp - 0.5, color=color, linewidth=linewidth, linestyle=linestyle)
    
    # add mean line
    axarr[0].axhline(np.nanmedian(signal_veg), label='Mean', linestyle='--', color='orange')
    
    # set title, subtitle, axis labels, etc
    axarr[0].set_title(figtitle, loc='left')
    axarr[0].set_ylabel('Greenness')
    axarr[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # set color cycle for break fills
    color_cycle = cycle(["#4286f4", "#f44174"])
    
    # brightness plot
    # plot signal, color each (true) regime
    axarr[1].plot(range(n_samples), signal_bright, color='c')
    bkps = [0] + sorted(result_brt)
    
    # set start and end areas for background color of blue or red
    for (start, end), col in zip(rpt.utils.pairwise(bkps), color_cycle):
        axarr[1].axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
        
    # vertical lines to mark the result
    if result_brt is not None:
        for bkp in result_brt:
            if bkp != 0 and bkp < n_samples:
                axarr[1].axvline(x=bkp - 0.5, color=color, linewidth=linewidth, linestyle=linestyle)
    
    # add mean line
    axarr[1].axhline(np.nanmedian(signal_bright), label='Mean', linestyle='--', color='orange')
    
    # set title, subtitle, axis labels, etc
    axarr[1].set_title(figtitle, loc='left')
    axarr[1].set_ylabel('Brightness')
    axarr[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axarr[1].plot(range(n_samples), signal_bright)
        
    # change axis labels
    if len(dates) > 0:
        date_list = [pd.to_datetime(str(d.data)).strftime('%Y-%m') for d in dates]
        plt.xticks(list(range(0, n_samples)), date_list, rotation='vertical')

    # final clean up
    fig.tight_layout()
    
    return fig, axarr