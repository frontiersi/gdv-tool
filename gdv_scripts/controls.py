# scripts for creation of jupyter controls and ui
import ipywidgets as widgets
import ipyfilechooser as file_widgets
import os
import functools
from gdv_scripts import helpers

## HELPER FUNCTIONS ##
def disable_controls(dash):
    try:
        for v_child in dash.children:
            for h_child in v_child.children:
                h_child.disabled = True
    except Exception as e:
        print('Error occurred during disable_controls of type {0}. Stopping.'.format(e))
        raise e   
        
def change_time_slices(selected_type):
    
    # original options
    if selected_type == 'original':
        options = [('Q1 (DJF)', 'q1'), ('Q2 (MAM)', 'q2'), ('Q3 (JJA)', 'q3'), ('Q4 (SON)', 'q4'),
                   ('JAN', 1), ('FEB', 2), ('MAR', 3), ('APR', 4), ('MAY', 5), ('JUN', 6), ('JUL', 7), 
                   ('AUG', 8), ('SEP', 9), ('OCT', 10), ('NOV', 11), ('DEC', 12)]

    elif selected_type == 'seasonal':
        options = [('Quarter to Quarter', 'qt'), ('Month to Month', 'mt')]
        
    else:
        options = []
        
    return options
    

## EVENT HANDLERS ##

## CONTROL TYPE CONSTRUCTORS ##
# build dropdown widget, handler is event func
def build_dropdown(desc, opts, val, handler, params, tip, d_width='125px', l_width='100%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.Dropdown(description=desc, options=opts, value=val, style=style, layout=layout, 
                                 description_tooltip=tip)
        widge.observe(functools.partial(handler, params=params), names='value')
        
        return widge
    except:
        raise

# build selectionrangeslider widget, handler is event func
def build_select_slider(desc, opts, idx, handler, params, tip, d_width='125px', l_width='100%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.SelectionRangeSlider(description=desc, options=opts, index=idx, style=style, layout=layout, 
                                             description_tooltip=tip)
        widge.observe(functools.partial(handler, params=params), names='value') 
        
        return widge 
    except:
        raise

# build float slider widget, handler is event func  
def build_float_slider(desc, min_val, max_val, step_val, val, form, handler, params, tip, d_width='125px', l_width='100%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.FloatSlider(description=desc, min=min_val, max=max_val, step=step_val, value=val, 
                                    readout_format=form, style=style, layout=layout, description_tooltip=tip)
        widge.observe(functools.partial(handler, params=params), names='value') 
        
        return widge 
    except:
        raise     
        
# build int slider widget, handler is event func  
def build_int_slider(desc, min_val, max_val, step_val, val, form, handler, params, tip, d_width='125px', l_width='100%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.IntSlider(description=desc, min=min_val, max=max_val, step=step_val, value=val, 
                                  readout_format=form, style=style, layout=layout, description_tooltip=tip)
        widge.observe(functools.partial(handler, params=params), names='value') 
        
        return widge 
    except:
        raise  
          
# build radio button widget, handler is event func
def build_radio_button(desc, opts, val, handler, params, tip, d_width='125px', l_width='50%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.RadioButtons(description=desc, options=opts, value=val, style=style, layout=layout, 
                                     description_tooltip=tip)
        widge.observe(functools.partial(handler, params=params), names='value') 
        
        return widge 
    except:
        raise  

# build file dialog widget
def build_file_dialog(directory, filename, title, handler, params, d_width='125px', l_width='100%', l_margin='-5px 0px 0px 37px'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width, margin=l_margin)
        widge = file_widgets.FileChooser(directory, filename=filename, title=title, show_hidden=False, select_default=False)
        widge.layout = layout
        widge.rows=1
        widge.register_callback(functools.partial(handler, params=params))

        return widge
    except:
        raise
        
# build empty cell widget
def build_empty_cell(desc, val, d_width='125px', l_width='100%'):
    try:
        style = {'description_width': d_width}
        layout = widgets.Layout(width=l_width)
        widge = widgets.Label(description=desc, value=val, style=style, layout=layout)
        
        return widge 
    except:
        raise  
               
    
## BUILD SPECIFIC CONTROLS DASHBOARDS ##
# gdv likelihood dashboard
def build_sat_fetch_dashboard(params):
    
    # event handler for platform (ls or s2)
    def platform_handler(event, params):
        if event['type'] == 'change':
            params.platform = event['new']

    # event handler for product (nbar or nbart)
    def product_handler(event, params):
        if event['type'] == 'change':
            params.product = event['new']

    # event handler for query date range (converts year ints to YYYY-MM-DD strings)    
    def query_date_handler(event, params):
        if event['type'] == 'change':
            params.query_dates = helpers.prepare_dates(date_tuple=event['new'])

    # event handler for minimum cloud cover percentage
    def cloud_cover_handler(event, params):
        if event['type'] == 'change':
            params.cloud_cover = helpers.prepare_cloud_percentage(event['new'])

    # event handler for cloud masking boolean
    def cloud_mask_handler(event, params):
        if event['type'] == 'change':
            params.cloud_mask = event['new']
    
    
    # drop down for platform selection param
    tip = 'Select satellite platform. Landsat is recommended due to its larger temporal range.'
    w_platform = build_dropdown(desc='Platform: ', opts=[('Landsat 5, 7, 8', 'ls'), ('Sentinel 2AB', 's2')],
                                val='ls', handler=platform_handler, params=params, tip=tip)

    # drop down for product selection param
    tip = 'Select product. Both are atmospherically corrected. NBAR-T adds terrain illum. correction. NBAR-T recommended.'
    w_product = build_dropdown(desc='Product: ', opts=[('NBAR', 'nbar'), ('NBAR-T', 'nbart')],
                               val='nbart', handler=product_handler, params=params, tip=tip)

    # selection range slider for query years param
    tip = 'Select the date range for analysis. Ten+ years is recommended to reduce seasonal variation.'
    w_query_dates = build_select_slider(desc='Query Dates: ', opts=range(1980, 2030), idx=(19, 39), 
                                        handler=query_date_handler, params=params, tip=tip)

    # float slider for maximum cloud cover param
    tip = 'Select the maximum % of cloud cover. Around 0-2% is recommended for Landsat, and 0-5% for Sentinel.'
    w_cloud_cover = build_float_slider(desc='Max Cloud Cover:', min_val=0.0, max_val=1.0, step_val=0.01, 
                                       val=0.02, form='.0%', handler=cloud_cover_handler, params=params, tip=tip)

    # radiobutton for cloud mask to nan
    tip = 'Choose whether clouds are masked out. If yes, cloud pixels will be set to nan.'
    w_cloud_mask = build_radio_button(desc='Mask Clouds:', opts=[('Yes', True), ('No', False)], val=True, 
                                      handler=cloud_mask_handler, params=params, tip=tip)
    
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0')
    empty_cell = build_empty_cell(desc=' ', val=' ')
    
    # build dashboard
    w_dash_r0 = widgets.HBox([w_platform, w_product, w_query_dates], layout=control_layout)
    w_dash_r1 = widgets.HBox([w_cloud_cover, w_cloud_mask, empty_cell], layout=control_layout)
    w_dash_v0 = widgets.VBox([w_dash_r0, w_dash_r1])

    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_v0], selected_index=0)
    w_dash.set_title(0, 'Satellite platform, query date and cloud masking parameters')
    
    # return dashboard
    return w_dash

# generate indices dashboard
def build_indices_dashboard(params):
    
    # event handler for vegetation index (ndvi, savi, etc.)
    def veg_idx_handler(event, params):
        if event['type'] == 'change':
            params.veg_idx = event['new']

    # event handler for moisture index (ndmi, etc.)
    def mst_idx_handler(event, params):
        if event['type'] == 'change':
            params.mst_idx = event['new']

    # event handler for rescale boolean
    def rescale_handler(event, params):
        if event['type'] == 'change':
            params.rescale = event['new']

    # event handler for soil adjustment (l) value 
    def l_value_handler(event, params):
        if event['type'] == 'change':
            params.l_value = event['new']

    # event handler for show raw index figure boolean
    def plot_raw_idx_handler(event, params):
        if event['type'] == 'change':
            params.plot_raw_idx = event['new']        
    
    
    # drop down for vegetation index selection param
    veg_opts = [('NDVI', 'ndvi'), ('SAVI', 'savi'), ('STVI3', 'stvi3'), ('SLAVI', 'slavi'), ('MSVI1', 'msvi1'), 
                ('MSVI2', 'msvi2'), ('MSVI3', 'msvi3'), ('MAVI', 'mavi'), ('TCAP Green', 'tcap')]
    tip = 'Select vegetation index. NDVI, SAVI, SLAVI and TCAP Green recommended.'
    w_veg_idx = build_dropdown(desc='Vegetation Index: ', opts=veg_opts, val='ndvi', 
                               handler=veg_idx_handler, params=params, tip=tip)

    # drop down for moisture index selection param
    tip = 'Select moisture index. NDMI and GVMI are the only indices currently supported.'
    w_mst_idx = build_dropdown(desc='Moisture Index: ', opts=[('NDMI', 'ndmi'), ('GVMI', 'gvmi')], val='ndmi', 
                               handler=mst_idx_handler, params=params, tip=tip)
    
    # radiobutton for rescale index values param
    tip = 'Choose whether to rescale index values where -1to1 to 0to2 or / by 10000 whenever DEA rescaled. Recommended.'
    w_rescale = build_radio_button(desc='Rescale Values:', opts=[('Yes', True), ('No', False)], val=True, 
                                   handler=rescale_handler, params=params, tip=tip)

    # float slider for soil l value param
    tip = 'Select the soil adjustment (L) value. Applicable only for SAVI.'
    w_l_value = build_float_slider(desc='Soil Adjustment (L):', min_val=0.0, max_val=1.0, step_val=0.1, 
                                   val=0.5, form='.1f', handler=l_value_handler, params=params, tip=tip)  

    # radiobutton for plotting of index values param
    tip = 'Show a panel of figure of all images once data computed computed.'
    w_plot_raw_idx= build_radio_button(desc='Show Figures:', opts=[('Yes', True), ('No', False)], val=False, 
                                   handler=plot_raw_idx_handler, params=params, tip=tip)
    
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0')
    empty_cell = build_empty_cell(desc=' ', val=' ')
    
    # build dashboard
    w_dash_r0 = widgets.HBox([w_veg_idx, w_mst_idx, w_l_value, w_rescale], layout=control_layout)
    w_dash_r1 = widgets.HBox([w_plot_raw_idx, empty_cell, empty_cell, empty_cell], layout=control_layout)

    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r0, w_dash_r1], selected_index=0)
    w_dash.set_title(0, 'Vegetation and moisture index parameters')
    w_dash.set_title(1, 'Show figure parameters')
    
    # return dashboard
    return w_dash

# generate interpolate dashboard
def build_interpolation_dashboard(params):
    
    # event handler for do z-score test boolean
    def zscore_handler(event, params):
        if event['type'] == 'change':
            params.zscore = event['new']

    # event handler for z-score critical value
    def zscore_value_handler(event, params):
        if event['type'] == 'change':
            params.zscore_value = event['new']

    # event handler for missing data handler value
    def miss_data_handle_handler(event, params):
        if event['type'] == 'change':
            params.miss_data_handle = event['new']

    # event handler for filling missing pixels boolean
    def fill_pixels_handler(event, params):
        if event['type'] == 'change':
            params.fill_pixels = event['new'] 

    # event handler for show interp index figure boolean
    def plot_interp_idx_handler(event, params):
        if event['type'] == 'change':
            params.plot_interp_idx = event['new']  
    
    # checkbox for z-score test param
    tip = 'Eliminate outlier scenes using Z-Test. Default set to Yes.'
    w_zscore = build_radio_button(desc='Perform Z-Test:', opts=[('Yes', True), ('No', False)], val=True, 
                                  handler=zscore_handler, params=params, tip=tip)
    
    # drop down for z-score critical value param
    tip = 'Select the Z-score critical value. Options are 1.65, 1.96, 2.58. Equal to pvalue of 0.1, 0.05, 0.01, respectively. '
    tip += 'Recommended 1.96 (p-value 0.05). A Z-score of 1.96 may consider too data to be outliers, 2.58 too few.'
    w_zscore_value = build_dropdown(desc='Z-Score Value: ', opts=[('1.65', 1.65), ('1.96', 1.96), ('2.58', 2.58)], 
                                    val=1.96, handler=zscore_value_handler, params=params, tip=tip)

    # drop down for handling of missing data param
    tip = 'Set how to handle missing data. Either drop it (whole year will be lost), or interpolate it (linear). '
    tip += 'Interpolation is recommended. Dropping data can result in significant data loss.'
    w_miss_data_handle = build_dropdown(desc='Handle Missing Data: ', opts=[('Interpolate', 'interpolate'), ('Drop', 'drop')], 
                                        val='interpolate', handler=miss_data_handle_handler, params=params, tip=tip)

    # radiobutton for show figure values param
    tip = 'Use interpolation to fill in missing pixels if present.'
    w_fill_pixels = build_radio_button(desc='Fill empty pixels:', opts=[('Yes', True), ('No', False)], val=False, 
                                          handler=fill_pixels_handler, params=params, tip=tip)
    
    # radiobutton for show figure values param
    tip = 'Show a panel of figure of all images once data computed computed.'
    w_plot_interp_idx = build_radio_button(desc='Show Figures:', opts=[('Yes', True), ('No', False)], val=False, 
                                          handler=plot_interp_idx_handler, params=params, tip=tip)
   
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0')
    empty_cell = build_empty_cell(desc=' ', val=' ')

    # build dashboard
    w_dash_r0 = widgets.HBox([w_miss_data_handle, w_fill_pixels, empty_cell], layout=control_layout)
    w_dash_r1 = widgets.HBox([w_zscore, w_zscore_value, empty_cell], layout=control_layout)
    w_dash_r2 = widgets.HBox([w_plot_interp_idx, empty_cell, empty_cell], layout=control_layout)

    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r0, w_dash_r1, w_dash_r2], selected_index=0)
    w_dash.set_title(0, 'Interpolation parameters')
    w_dash.set_title(1, 'Z-score analysis parameters')
    w_dash.set_title(2, 'Show figure parameters')
    
    # return dashboard
    return w_dash

# generate standardisation preparation dashboard
def build_prepare_standardisation_dashboard(params):
    
    # event handler for standardisation greennest/moistest percentile float
    def green_moist_percentile_handler(event, params):
        if event['type'] == 'change':
            params.green_moist_percentile = event['new']  

    # event handler for standardisation greennest/moistest percentile float
    def slope_steady_percentile_handler(event, params):
        if event['type'] == 'change':
            params.slope_steady_percentile = event['new']  

    # event handler for standardisation greennest/moistest percentile float
    def max_num_sites_handler(event, params):
        if event['type'] == 'change':
            params.max_num_sites = event['new']          

    # event handler for show standardisation index figure boolean
    def plot_standard_idx_handler(event, params):
        if event['type'] == 'change':
            params.plot_standard_idx = event['new']  
    
    # float slider for greenest/moistestpercentile value param
    tip = 'Select the percentile for obtaining the greennest/moistest pixels. We recommend using the default value (99.5%). '
    tip += 'However, if some features in your study area drive up the vegetation/moisture values, reduce this value 1 or 2%.'
    w_green_moist_percentile = build_float_slider(desc='Greennest/Moistest Percentile:', min_val=80.0, max_val=99.99, step_val=0.1, d_width='175px',
                                                  val=99.5, form='.1f', handler=green_moist_percentile_handler, params=params, tip=tip) 

    # int slider for steadiest slope value param
    tip = 'Select the percentile for obtaining the steadiest pixels. We recommend using the default value (5%). This will obtain areas of the slope '
    tip += 'between with 5% either side of slope = 0 (i.e. no changing pixels). If not enough pixels are returned below, try increasing by 1-2% and re-run.'
    w_slope_steady_percentile = build_int_slider(desc='Steadiest Slope Percentile:', min_val=1, max_val=30, step_val=1, d_width='175px',
                                                   val=5, form='d', handler=slope_steady_percentile_handler, params=params, tip=tip) 

    # float slider for maximum slope value param
    tip = 'Set the maximum number of invariant green/moist sites to be used in standardisation. Default is 50. Fewer sites can provide slightly more accurate '
    tip += 'standardisation, but may suffer if pixels fluctuate across years. Opposite is true for more sites. We recommend 50.'
    w_max_num_sites = build_int_slider(desc='Max Number Sites:', min_val=25, max_val=100, step_val=5, val=50, form='d', d_width='175px',
                                       handler=max_num_sites_handler, params=params, tip=tip) 
    
    # radiobutton for show figure values param
    tip = 'Show a panel of figure of all images once data computed computed.'
    w_plot_standard_idx = build_radio_button(desc='Show Figures:', opts=[('Yes', True), ('No', False)], val=False, d_width='175px', 
                                          handler=plot_standard_idx_handler, params=params, tip=tip)
        
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0')
    empty_cell = build_empty_cell(desc=' ', val=' ')

    # build dashboard
    w_dash_r0 = widgets.HBox([w_green_moist_percentile, w_slope_steady_percentile, w_max_num_sites], layout=control_layout)
    w_dash_r1 = widgets.HBox([w_plot_standard_idx, empty_cell, empty_cell], layout=control_layout)

    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r0, w_dash_r1], selected_index=0)
    w_dash.set_title(0, 'Standardisation invariant site parameters')
    w_dash.set_title(1, 'Show figure parameters')
    
    # return dashboard
    return w_dash
    
# generate likelihood threshold dashboard
def build_likelihood_threshold_dashboard(params, years_list):
    
    # event handler for likelihood threshold type radio
    def like_thresh_type_handler(event, params):
        if event['type'] == 'change':
            params.like_thresh_type = event['new']       

    # event handler for likelihood threshold date drop down
    def like_thresh_date_handler(event, params):
        if event['type'] == 'change':
            params.like_thresh_date = event['new']              

    # event handler for likelihood threshold standard deviation float
    def like_thresh_stdv_handler(event, params):
        if event['type'] == 'change':
            params.like_thresh_stdv = event['new']             

    # event handler for mapping likelihood boolean
    def map_likelihood_handler(event, params):
        if event['type'] == 'change':
            params.map_likelihood = event['new']
    
    # event handler for ground truth file dialog
    def ground_sites_path_handler(params):
        params.ground_sites_path = w_ground_sites_path.selected
        
    # event handler for write data dialog
    def write_data_path_handler(params):
        params.write_data_path = w_write_data_path.selected
    
    # drop down for likelihood threshold type param
    tip = 'Select thresholding type. If groundtruth data for a particular year is available, select Groundtruth Shapefile. If not, select Standard Deviation.'   
    w_like_thresh_type = build_dropdown(desc='Threshold Type:', opts=[('Standard Deviation', 'stdv'), ('Groundtruth Shapefile', 'shape')], 
                                        val='stdv', handler=like_thresh_type_handler, params=params, tip=tip)
    
    # filedialog to select shapefiles
    w_ground_sites_path = build_file_dialog(directory=os.getcwd(), filename='', title='Select Groundtruth Points (.shp):', 
                                            handler=ground_sites_path_handler, params=params)
        
    # drop down for likelihood threshold date
    tip = 'Select date to apply threshold. If groundtruth data is provided, it is best to match the year the data was collected. Ultimately, we recommend using the all option. '
    tip += 'This will find the median of all likelihood maps and threshold this via standard deviation. This will reduce seasonal variation.'
    w_like_thresh_date = build_dropdown(desc='Threshold Date:', opts=years_list, val='all', 
                                        handler=like_thresh_date_handler, params=params, tip=tip)
    
    # float slider for standard deviation value param
    tip = 'Select the standard deviation value for thresholding. We recommend 1.5 to 2. Lower values will return more GDV likelihood pixels. A higher value returns the converse.'
    w_like_thresh_stdv = build_float_slider(desc='Standard Deviation:', min_val=0.5, max_val=3.5, step_val=0.5, 
                                            val=1.5, form='.1f', handler=like_thresh_stdv_handler, params=params, tip=tip)
    
    # radiobutton for show figure values param
    tip = 'Show an interactive map of likelihood (thresholded and non-thresholded) of selected date once data computed.'
    w_map_likelihood = build_radio_button(desc='Show Map:', opts=[('Yes', True), ('No', False)], val=True, 
                                          handler=map_likelihood_handler, params=params, tip=tip)
    
    # write data path
    w_write_data_path = build_file_dialog(directory=os.getcwd(), filename='', title='Export Folder:', 
                                          handler=write_data_path_handler, params=params)
    
    # create throw-away html title for controls 
    lbl = widgets.HTML(value='<h4>Set groundwater dependent vegetation likelihood thresholding parameters and data export options: </h4>')
    empty_cell = build_empty_cell(desc=' ', val=' ')
    
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0px')
    
    # build dashboard
    w_dash_r0 = widgets.HBox([w_ground_sites_path], layout=control_layout)
    w_dash_r1 = widgets.HBox([w_like_thresh_type, w_like_thresh_date, w_like_thresh_stdv, w_map_likelihood], layout=control_layout)
    w_dash_r2 = widgets.HBox([w_write_data_path], layout=control_layout)
    
    # create dash as accordion
    #w_dash = widgets.Accordion(children=[w_dash_r5, w_dash_r0, w_dash_r1, w_dash_r2, w_dash_r3, w_dash_r4], selected_index=1)
    w_dash = widgets.Accordion(children=[w_dash_r0, w_dash_r1, w_dash_r2], selected_index=1)
    w_dash.set_title(0, 'Groundtruthed shapefile import parameters')
    w_dash.set_title(1, 'Likelihood modelling parameters')
    w_dash.set_title(2, 'Export likelihood data parameters')

    # return dashboard
    return w_dash

# generate trend dashboard
def build_trend_dashboard(params):
    
    # event handler for trend platform (ls or s2)
    def trend_platform_handler(event, params):
        if event['type'] == 'change':
            params.platform = event['new']

    # event handler for trend product (nbar or nbart)
    def trend_product_handler(event, params):
        if event['type'] == 'change':
            params.product = event['new']

    # event handler for trend query date range (converts year ints to YYYY-MM-DD strings)    
    def trend_query_date_handler(event, params):
        if event['type'] == 'change':
            params.query_dates = helpers.prepare_dates(date_tuple=event['new'])

    # event handler for trend minimum cloud cover percentage
    def trend_cloud_cover_handler(event, params):
        if event['type'] == 'change':
            params.cloud_cover = helpers.prepare_cloud_percentage(event['new'])      

    # event handler for trend cloud masking boolean
    def trend_cloud_mask_handler(event, params):
        if event['type'] == 'change':
            params.cloud_mask = event['new']

    # event handler for mankendall type handler
    def trend_mk_type_handler(event, params):
        if event['type'] == 'change':
            params.mk_type = event['new']
            w_time_slice.index = None
            w_time_slice.options = change_time_slices(params.mk_type)
    
    # event handler for trend time slice  handler
    def trend_time_slice_handler(event, params):
        if event['type'] == 'change':
            params.time_slice = event['new']

    # event handler for mankendall significant only handler
    def trend_sig_only_handler(event, params):
        if event['type'] == 'change':
            params.mk_sig_only = event['new']
    
    # event handler for mannkendall trend type handler
    def trend_trend_type_handler(event, params):
        if event['type'] == 'change':
            params.trend_type = event['new']
            
    # event handler for trend vegetation index (ndvi, savi, etc.)
    def trend_veg_idx_handler(event, params):
        if event['type'] == 'change':
            params.veg_idx = event['new']
            
    # event handler for trend rescale boolean
    def trend_rescale_handler(event, params):
        if event['type'] == 'change':
            params.rescale = event['new']
            
    # event handler for trend soil adjustment (l) value 
    def trend_l_value_handler(event, params):
        if event['type'] == 'change':
            params.l_value = event['new']

    # event handler for trend do z-score test boolean
    def trend_zscore_handler(event, params):
        if event['type'] == 'change':
            params.zscore = event['new']

    # event handler for trend z-score critical value
    def trend_zscore_value_handler(event, params):
        if event['type'] == 'change':
            params.zscore_value = event['new']

    # event handler for trend missing data handler value
    def trend_miss_data_handle_handler(event, params):
        if event['type'] == 'change':
            params.miss_data_handle = event['new']

    # event handler for trend filling missing pixels boolean
    def trend_fill_pixels_handler(event, params):
        if event['type'] == 'change':
            params.fill_pixels = event['new'] 
            
    # event handler for trend standardisation greennest/moistest percentile float
    def trend_green_moist_percentile_handler(event, params):
        if event['type'] == 'change':
            params.green_percentile = event['new']  

    # event handler for trend standardisation greennest/moistest percentile float
    def trend_max_num_sites_handler(event, params):
        if event['type'] == 'change':
            params.max_num_sites = event['new']          
            
    # event handler for trend likelihood import dialog
    def trend_read_data_path_handler(params):
        params.read_data_path = w_read_data_path.selected
            
    # drop down for platform selection param
    tip = 'Select satellite platform. Landsat is recommended due to its larger temporal range.'
    w_platform = build_dropdown(desc='Platform: ', opts=[('Landsat 5, 7, 8', 'ls'), ('Sentinel 2AB', 's2')],
                                val='ls', handler=trend_platform_handler, params=params, tip=tip)

    # drop down for product selection param
    tip = 'Select product. Both are atmospherically corrected. NBAR-T adds terrain illum. correction. NBAR-T recommended.'
    w_product = build_dropdown(desc='Product: ', opts=[('NBAR', 'nbar'), ('NBAR-T', 'nbart')],
                               val='nbart', handler=trend_product_handler, params=params, tip=tip)

    # selection range slider for query years param
    tip = 'Select the date range for analysis. Five+ years is required for trend analysis.'
    w_query_dates = build_select_slider(desc='Trend Dates: ', opts=range(1980, 2030), idx=(29, 39), 
                                        handler=trend_query_date_handler, params=params, tip=tip)

    # float slider for maximum cloud cover param
    tip = 'Select the maximum % of cloud cover. Around 0-2% is recommended for Landsat, and 0-5% for Sentinel.'
    w_cloud_cover = build_float_slider(desc='Max Cloud Cover:', min_val=0.0, max_val=1.0, step_val=0.01, 
                                       val=0.02, form='.0%', handler=trend_cloud_cover_handler, params=params, tip=tip)

    # radiobutton for cloud mask to nan
    tip = 'Choose whether clouds are masked out. If yes, cloud pixels will be set to nan.'
    w_cloud_mask = build_radio_button(desc='Mask Clouds:', opts=[('Yes', True), ('No', False)], val=True, 
                                      handler=trend_cloud_mask_handler, params=params, tip=tip)
    
    # radiobutton for trend type values param
    tip = 'Perform a Mann-Kendall trend analysis using seasonal or original processor. Seasonal is recommended.'
    w_mk_type = build_radio_button(desc='Mann-Kendall Type:', opts=[('Seasonal', 'seasonal'), ('Original', 'original')], val='seasonal', 
                                          handler=trend_mk_type_handler, params=params, tip=tip)
    
    # drop down for time slice selection param
    tip = 'Select time slice to which to perform trend alaysis across. For example: selecting MAR will look at trends across years for every March scene.'
    w_time_slice = build_dropdown(desc='Time Slice: ', opts=[('Quarter to Quarter', 'qt'), ('Month to Month', 'mt')], val='qt', 
                                  handler=trend_time_slice_handler, params=params, tip=tip)
    
    # radiobutton for trend sigificant only param
    tip = 'Significant (p-value < 0.05) trend pixels returned only (select Yes), or all trends returned regardless of significance (select No).'
    w_sig_only = build_radio_button(desc='Sig. Trend Only:', opts=[('Yes', True), ('No', False)], val=False, 
                                          handler=trend_sig_only_handler, params=params, tip=tip)
    
    # drop down for trend type selection param
    tip = 'Show trend that is increasing only, decreasing only, or any (increase, decreasing, no trend).'
    w_trend_type = build_dropdown(desc='Show Trend: ', opts=[('Increasing', 'increasing'), ('Decreasing', 'decreasing'), ('All', 'all')], val='all', 
                                  handler=trend_trend_type_handler, params=params, tip=tip)
    
    # drop down for vegetation index selection param
    veg_opts = [('NDVI', 'ndvi'), ('SAVI', 'savi'), ('STVI3', 'stvi3'), ('SLAVI', 'slavi'), ('MSVI1', 'msvi1'), 
                ('MSVI2', 'msvi2'), ('MSVI3', 'msvi3'), ('MAVI', 'mavi'), ('TCAP Green', 'tcap')]
    tip = 'Select vegetation index. NDVI, SAVI, SLAVI and TCAP Green recommended.'
    w_veg_idx = build_dropdown(desc='Vegetation Index: ', opts=veg_opts, val='ndvi', 
                               handler=trend_veg_idx_handler, params=params, tip=tip)
    
    # radiobutton for rescale index values param
    tip = 'Choose whether to rescale index values where -1to1 to 0to2 or / by 10000 whenever DEA rescaled. Recommended.'
    w_rescale = build_radio_button(desc='Rescale Values:', opts=[('Yes', True), ('No', False)], val=True, 
                                   handler=trend_rescale_handler, params=params, tip=tip)

    # float slider for soil l value param
    tip = 'Select the soil adjustment (L) value. Applicable only for SAVI.'
    w_l_value = build_float_slider(desc='Soil Adjustment (L):', min_val=0.0, max_val=1.0, step_val=0.1, 
                                   val=0.5, form='.1f', handler=trend_l_value_handler, params=params, tip=tip)  
    
    # checkbox for z-score test param
    tip = 'Eliminate outlier scenes using Z-score test. Default set to Yes.'
    w_zscore = build_radio_button(desc='Perform Z-Test:', opts=[('Yes', True), ('No', False)], val=True, 
                                  handler=trend_zscore_handler, params=params, tip=tip)
    
    # drop down for z-score critical value param
    tip = 'Select the Z-score value. Options are 1.65, 1.96, 2.58. Equal to pvalue of 0.1, 0.05, 0.01, respectively. '
    tip += 'Recommended 1.96 (p-value 0.05). A Z-score of 1.96 may consider too data to be outliers, 2.58 too few.'
    w_zscore_value = build_dropdown(desc='Z-Score Value: ', opts=[('1.65', 1.65), ('1.96', 1.96), ('2.58', 2.58)], 
                                    val=1.96, handler=trend_zscore_value_handler, params=params, tip=tip)

    # drop down for handling of missing data param
    tip = 'Set how to handle missing data. Either drop it (whole year will be lost), or interpolate it (linear). '
    tip += 'Interpolation is recommended. Dropping data can result in significant data loss.'
    w_miss_data_handle = build_dropdown(desc='Handle Missing Data: ', opts=[('Interpolate', 'interpolate'), ('Drop', 'drop')], 
                                        val='interpolate', handler=trend_miss_data_handle_handler, params=params, tip=tip)

    # radiobutton for show figure values param
    tip = 'Use interpolation to fill in missing pixels if present.'
    w_fill_pixels = build_radio_button(desc='Fill empty pixels:', opts=[('Yes', True), ('No', False)], val=False, 
                                       handler=trend_fill_pixels_handler, params=params, tip=tip)
    
    # float slider for greenest percentile value param
    tip = 'Select the percentile for obtaining the greennest/moistest pixels. We recommend using the default value (99.5%). '
    tip += 'However, if some features in your study area drive up the vegetation values, reduce this value 1 or 2%.'
    w_green_moist_percentile = build_float_slider(desc='Greenest Percentile:', min_val=80.0, max_val=99.99, step_val=0.1, d_width='175px',
                                                  val=99.5, form='.1f', handler=trend_green_moist_percentile_handler, params=params, tip=tip) 

    # float slider for maximum slope value param
    tip = 'Set the maximum number of invariant green/moist sites to be used in standardisation. Default is 50. Fewer sites can provide slightly more accurate '
    tip += 'standardisation, but may suffer if pixels fluctuate across years. Opposite is true for more sites. We recommend 50.'
    w_max_num_sites = build_int_slider(desc='Max Number Sites:', min_val=25, max_val=100, step_val=5, val=50, form='d', d_width='175px',
                                       handler=trend_max_num_sites_handler, params=params, tip=tip) 
    
    # file dialog for read netcdf data path
    w_read_data_path = build_file_dialog(directory=os.getcwd(), filename='', title='Select Existing Likelihood Data:', 
                                         handler=trend_read_data_path_handler, params=params)
        
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0px')
    empty_cell = build_empty_cell(desc=' ', val=' ')
        
    # build dashboard
    w_dash_r0 = widgets.HBox([w_platform, w_product, w_cloud_cover, w_cloud_mask], layout=control_layout)
    w_dash_r1 = widgets.VBox([widgets.HBox([w_mk_type, w_time_slice, w_query_dates], layout=control_layout), 
                              widgets.HBox([w_sig_only, w_trend_type, empty_cell], layout=control_layout)])
    w_dash_r2 = widgets.HBox([w_veg_idx, w_rescale, w_l_value])
    w_dash_r3 = widgets.HBox([w_zscore, w_zscore_value, w_miss_data_handle, w_fill_pixels])
    w_dash_r4 = widgets.HBox([w_green_moist_percentile, w_max_num_sites])
    w_dash_r5 = widgets.HBox([w_read_data_path])
    
    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r5, w_dash_r0, w_dash_r1, w_dash_r2, w_dash_r3, w_dash_r4], selected_index=1)
    w_dash.set_title(0, 'Likelihood import parameters')
    w_dash.set_title(1, 'Satellite data parameters')
    w_dash.set_title(2, 'Trend analysis parameters')
    w_dash.set_title(3, 'Vegetation index parameters')
    w_dash.set_title(4, 'Interpolation and outlier correction parameters')
    w_dash.set_title(5, 'Basic standardisation parameters')

    # return dashboard
    return w_dash

# generate trend change point detection dashboard
def build_trend_cpd_dashboard(params, handle_button):
        
    # event handler for cpd model selection
    def cpd_model_handler(event, params):
        if event['type'] == 'change':
            params.model = event['new']

    # event handler for cpd min size selection
    def cpd_min_size_handler(event, params):
        if event['type'] == 'change':
            params.min_size = event['new']
            
    # event handler for jump selection
    def cpd_jump_handler(event, params):
        if event['type'] == 'change':
            params.jump = event['new']
            
    # event handler for pelt breakpoint penelty selection
    def cpd_pelt_penelty_handler(event, params):
        if event['type'] == 'change':
            params.pelt_penelty = event['new']
            
    # event handler for dynp breakpoint number selection
    def cpd_dynp_n_bkps_handler(event, params):
        if event['type'] == 'change':
            params.dynp_n_bkps = event['new']
            

    # drop down for model selection param
    tip = 'Select change point detection model. L1: detects changes in the median of a signal. L2: Cost function detects mean-shifts in a signal. '
    tip += 'RBF: detects changes in the mean of the embedded signal. Recommend RBF.'
    w_model = build_dropdown(desc='Model: ', opts=[('Least absolute deviation (l1)', 'l1'), ('Least squared deviation (l2)', 'l2'), ('Kernelized mean change (rbf)', 'rbf')],
                                val='rbf', handler=cpd_model_handler, params=params, tip=tip)
    
    # int slider for min size between rbeaks value param
    tip = 'Minimum distance between change points (e.g. minimum number of dates between detected breaks). Default is 3.'
    w_min_size = build_int_slider(desc='Min Break Distance:', min_val=1, max_val=25, step_val=1, val=2, 
                                  handler=cpd_min_size_handler, form='d', params=params, tip=tip) 
    
    # int slider for jump value param
    tip = 'Consider only change point locations that are multiple of this particular value. Default is 1.'
    w_jump = build_int_slider(desc='Jump Value:', min_val=1, max_val=25, step_val=1, val=1, 
                                  handler=cpd_jump_handler, form='d', params=params, tip=tip)

    # int slider for pelt penelty param
    tip = 'Set the penelty of break points for pelt change point detection method. Default is 3 for penelty.'
    w_pelt_penelty = build_int_slider(desc='Break Penelty (Pelt):', min_val=1, max_val=25, step_val=1, val=3, 
                                      handler=cpd_pelt_penelty_handler, form='d', params=params, tip=tip)
    
    # int slider for dynamic programming num breaks param
    tip = 'Set the number of break points for dynamic programming change point detection method. Default is 4 breakpoints.'
    w_dynp_n_bkps = build_int_slider(desc='Num Breaks (DynPro):', min_val=2, max_val=24, step_val=2, val=2, 
                                  handler=cpd_dynp_n_bkps_handler, form='d', params=params, tip=tip)
    
    # button for cpd processing
    w_trend_cpd_button = widgets.Button(description='Show')
    w_trend_cpd_button.on_click(handle_button)

    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0px')
    empty_cell = build_empty_cell(desc=' ', val=' ')
        
    # build dashboard
    w_dash_r0 = widgets.VBox([widgets.HBox([w_model, w_min_size, w_jump], layout=control_layout), 
                              widgets.HBox([w_dynp_n_bkps, w_pelt_penelty], layout=control_layout),
                              widgets.HBox([empty_cell, empty_cell, w_trend_cpd_button], layout=control_layout)])
    
    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r0], selected_index=0)
    w_dash.set_title(0, 'Change point detection parameters')

    # return dashboard
    return w_dash

# generate cva dashboard
def build_change_dashboard(params):
    
    # event handler for change platform (ls or s2)
    def change_platform_handler(event, params):
        if event['type'] == 'change':
            params.platform = event['new']

    # event handler for change product (nbar or nbart)
    def change_product_handler(event, params):
        if event['type'] == 'change':
            params.product = event['new']

    # event handler for change baseline date range (converts year ints to YYYY-MM-DD strings)    
    def change_base_date_handler(event, params):
        if event['type'] == 'change':
            params.base_dates = helpers.prepare_dates(date_tuple=event['new'])

    # event handler for change comparison date range (converts year ints to YYYY-MM-DD strings)    
    def change_comp_date_handler(event, params):
        if event['type'] == 'change':
            params.comp_dates = helpers.prepare_dates(date_tuple=event['new'])
             
    # event handler for change minimum cloud cover percentage
    def change_cloud_cover_handler(event, params):
        if event['type'] == 'change':
            params.cloud_cover = helpers.prepare_cloud_percentage(event['new'])      

    # event handler for change cloud masking boolean
    def change_cloud_mask_handler(event, params):
        if event['type'] == 'change':
            params.cloud_mask = event['new']
    
    # event handler for change time slice  handler
    def change_time_slice_handler(event, params):
        if event['type'] == 'change':
            params.time_slice = event['new']
            
    # event handler for change angle range
    def change_angles_range_handler(event, params):
        if event['type'] == 'change':
            params.change_angles_range = event['new']
            
    # event handler for change threshold value
    def change_threshold_handler(event, params):
        if event['type'] == 'change':
            params.change_threshold = event['new']            
            
    # event handler for change vegetation index (tcap only)
    def change_veg_idx_handler(event, params):
        if event['type'] == 'change':
            params.veg_idx = event['new']
            
    # event handler for change rescale boolean
    def change_rescale_handler(event, params):
        if event['type'] == 'change':
            params.rescale = event['new']
            
    # event handler for change do z-score test boolean
    def change_zscore_handler(event, params):
        if event['type'] == 'change':
            params.zscore = event['new']

    # event handler for change z-score critical value
    def change_zscore_value_handler(event, params):
        if event['type'] == 'change':
            params.zscore_value = event['new']
            
    # event handler for change missing data handler value
    def change_miss_data_handle_handler(event, params):
        if event['type'] == 'change':
            params.miss_data_handle = event['new']

    # event handler for change filling missing pixels boolean
    def change_fill_pixels_handler(event, params):
        if event['type'] == 'change':
            params.fill_pixels = event['new'] 
            
    # event handler for change standardisation greennest/moistest percentile float
    def change_green_moist_percentile_handler(event, params):
        if event['type'] == 'change':
            params.green_percentile = event['new']  

    # event handler for trend standardisation greennest/moistest percentile float
    def change_max_num_sites_handler(event, params):
        if event['type'] == 'change':
            params.max_num_sites = event['new']          
            
    # event handler for change likelihood import dialog
    def change_read_data_path_handler(params):
        params.read_data_path = w_read_data_path.selected
            
    # drop down for platform selection param
    tip = 'Select satellite platform. Landsat is recommended due to its larger temporal range.'
    w_platform = build_dropdown(desc='Platform: ', opts=[('Landsat 5, 7, 8', 'ls'), ('Sentinel 2AB', 's2')],
                                val='ls', handler=change_platform_handler, params=params, tip=tip)

    # drop down for product selection param
    tip = 'Select product. Both are atmospherically corrected. NBAR-T adds terrain illum. correction. NBAR-T recommended.'
    w_product = build_dropdown(desc='Product: ', opts=[('NBAR', 'nbar'), ('NBAR-T', 'nbart')],
                               val='nbart', handler=change_product_handler, params=params, tip=tip)

    # selection range slider for baseline query years param
    tip = 'Select the baseline date range for analysis. 20+ years are recommended.'
    w_base_dates = build_select_slider(desc='Baseline Dates: ', opts=range(1980, 2030), idx=(9, 39), 
                                       handler=change_base_date_handler, params=params, tip=tip)

    # selection range slider for comparison query years param
    tip = 'Select the comparison date range for analysis. 5+ years are recommended.'
    w_comp_dates = build_select_slider(desc='Comparison Dates: ', opts=range(1980, 2030), idx=(29, 39), 
                                       handler=change_comp_date_handler, params=params, tip=tip)
    
    # float slider for maximum cloud cover param
    tip = 'Select the maximum % of cloud cover. Around 0-2% is recommended for Landsat, and 0-5% for Sentinel.'
    w_cloud_cover = build_float_slider(desc='Max Cloud Cover:', min_val=0.0, max_val=1.0, step_val=0.01, 
                                       val=0.02, form='.0%', handler=change_cloud_cover_handler, params=params, tip=tip)

    # radiobutton for cloud mask to nan
    tip = 'Choose whether clouds are masked out. If yes, cloud pixels will be set to nan.'
    w_cloud_mask = build_radio_button(desc='Mask Clouds:', opts=[('Yes', True), ('No', False)], val=True, 
                                      handler=change_cloud_mask_handler, params=params, tip=tip)
    
    # drop down for time slice selection param
    tip = 'Select time slice to which to perform change vector analysis across. For example: selecting MAR will look at trends across years for every March scene.'
    opts = [('Q1 (DJF)', 'q1'), ('Q2 (MAM)', 'q2'), ('Q3 (JJA)', 'q3'), ('Q4 (SON)', 'q4'), ('JAN', 1), ('FEB', 2), ('MAR', 3), 
            ('APR', 4), ('MAY', 5), ('JUN', 6), ('JUL', 7), ('AUG', 8), ('SEP', 9), ('OCT', 10), ('NOV', 11), ('DEC', 12)]
    w_time_slice = build_dropdown(desc='Time Slice: ', opts=opts, val='q1', handler=change_time_slice_handler, params=params, tip=tip)
    
    # selection range slider for change angle range param
    tip = 'Select the range of angles to be returned from the change vector analysis. 0-90 is moisture reduction. 90-180 is vegetation decline. '
    tip += '180-270 is moisture increase. 270-359 is vegetation increase.'
    w_change_angles_range = build_select_slider(desc='Angle Range: ', opts=range(0, 360), idx=(90, 180), 
                                                handler=change_angles_range_handler, params=params, tip=tip)
    
    # drop down for magnitude threshold drop down param
    tip = 'Select threshold median factor. Used to calculate a threshold magnitude for which pixels are considered stable, i.e. no change. ' 
    tip += 'Defaults to 1.5 times the median non-zero magnitude.'
    w_change_threshold = build_dropdown(desc='Magnitude Threshold: ', opts=[('0.5', 0.5), ('1', 1), ('1.5', 1.5), ('2', 2), ('2.5', 2.5), ('3', 3)],
                               val=1.5, handler=change_threshold_handler, params=params, tip=tip)
    
    # drop down for vegetation index selection param
    tip = 'Select vegetation index. Only TCAP available for change vector analysis.'
    w_veg_idx = build_dropdown(desc='Vegetation Index: ', opts=[('TCAP', 'tcap')], val='tcap', 
                               handler=change_veg_idx_handler, params=params, tip=tip)
    
    # radiobutton for rescale index values param
    tip = 'Choose whether to rescale index values where -1to1 to 0to2 or / by 10000 whenever DEA rescaled. Recommended.'
    w_rescale = build_radio_button(desc='Rescale Values:', opts=[('Yes', True), ('No', False)], val=True, 
                                   handler=change_rescale_handler, params=params, tip=tip)
    
    # checkbox for z-score test param
    tip = 'Eliminate outlier scenes using Z-score test. Default set to Yes.'
    w_zscore = build_radio_button(desc='Perform Z-Test:', opts=[('Yes', True), ('No', False)], val=True, 
                                  handler=change_zscore_handler, params=params, tip=tip)
    
    # drop down for z-score critical value param
    tip = 'Select the Z-score value. Options are 1.65, 1.96, 2.58. Equal to pvalue of 0.1, 0.05, 0.01, respectively. '
    tip += 'Recommended 1.96 (p-value 0.05). A Z-score of 1.96 may consider too data to be outliers, 2.58 too few.'
    w_zscore_value = build_dropdown(desc='Z-Score Value: ', opts=[('1.65', 1.65), ('1.96', 1.96), ('2.58', 2.58)], 
                                    val=1.96, handler=change_zscore_value_handler, params=params, tip=tip)

    # drop down for handling of missing data param
    tip = 'Set how to handle missing data. Either drop it (whole year will be lost), or interpolate it (linear). '
    tip += 'Interpolation is recommended. Dropping data can result in significant data loss.'
    w_miss_data_handle = build_dropdown(desc='Handle Missing Data: ', opts=[('Interpolate', 'interpolate'), ('Drop', 'drop')], 
                                        val='interpolate', handler=change_miss_data_handle_handler, params=params, tip=tip)

    # radiobutton for show figure values param
    tip = 'Use interpolation to fill in missing pixels if present.'
    w_fill_pixels = build_radio_button(desc='Fill empty pixels:', opts=[('Yes', True), ('No', False)], val=False, 
                                       handler=change_fill_pixels_handler, params=params, tip=tip)
        
    # float slider for greenest percentile value param
    tip = 'Select the percentile for obtaining the greenest pixels. We recommend using the default value (99.5%). '
    tip += 'However, if some features in your study area drive up the vegetation values, reduce this value 1 or 2%.'
    w_green_moist_percentile = build_float_slider(desc='Greenest/Brightest Percentile:', min_val=80.0, max_val=99.99, step_val=0.1, d_width='175px',
                                                  val=99.5, form='.1f', handler=change_green_moist_percentile_handler, params=params, tip=tip) 

    # float slider for maximum slope value param
    tip = 'Set the maximum number of invariant green/moist sites to be used in standardisation. Default is 50. Fewer sites can provide slightly more accurate '
    tip += 'standardisation, but may suffer if pixels fluctuate across years. Opposite is true for more sites. We recommend 50.'
    w_max_num_sites = build_int_slider(desc='Max Number Sites:', min_val=25, max_val=100, step_val=5, val=50, form='d', d_width='175px',
                                       handler=change_max_num_sites_handler, params=params, tip=tip) 
    
    # file dialog for read netcdf data path
    w_read_data_path = build_file_dialog(directory=os.getcwd(), filename='', title='Select Existing Likelihood Data:', 
                                         handler=change_read_data_path_handler, params=params)
        
    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0px')
    empty_cell = build_empty_cell(desc=' ', val=' ')
        
    # build dashboard
    w_dash_r0 = widgets.HBox([w_platform, w_product, w_cloud_cover, w_cloud_mask], layout=control_layout)
    w_dash_r1 = widgets.VBox([widgets.HBox([w_base_dates, w_comp_dates, w_time_slice], layout=control_layout), 
                              widgets.HBox([w_change_angles_range, w_change_threshold, empty_cell], layout=control_layout)])
    w_dash_r2 = widgets.HBox([w_veg_idx, w_rescale, empty_cell])
    w_dash_r3 = widgets.HBox([w_zscore, w_zscore_value, w_miss_data_handle, w_fill_pixels])
    w_dash_r4 = widgets.HBox([w_green_moist_percentile, w_max_num_sites])
    w_dash_r5 = widgets.HBox([w_read_data_path])
    
    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r5, w_dash_r0, w_dash_r1, w_dash_r2, w_dash_r3, w_dash_r4], selected_index=1)
    w_dash.set_title(0, 'Likelihood import parameters')
    w_dash.set_title(1, 'Satellite data parameters')
    w_dash.set_title(2, 'Change vector analysis parameters')
    w_dash.set_title(3, 'Vegetation index parameters')
    w_dash.set_title(4, 'Interpolation and outlier correction parameters')
    w_dash.set_title(5, 'Basic standardisation parameters')

    # return dashboard
    return w_dash

# generate cva change point detection dashboard
def build_change_cpd_dashboard(params, handle_button):
        
    # event handler for cpd model selection
    def cpd_model_handler(event, params):
        if event['type'] == 'change':
            params.model = event['new']

    # event handler for cpd min size selection
    def cpd_min_size_handler(event, params):
        if event['type'] == 'change':
            params.min_size = event['new']
            
    # event handler for jump selection
    def cpd_jump_handler(event, params):
        if event['type'] == 'change':
            params.jump = event['new']
            
    # event handler for pelt breakpoint penelty selection
    def cpd_pelt_penelty_handler(event, params):
        if event['type'] == 'change':
            params.pelt_penelty = event['new']
            
    # event handler for dynp breakpoint number selection
    def cpd_dynp_n_bkps_handler(event, params):
        if event['type'] == 'change':
            params.dynp_n_bkps = event['new']
            

    # drop down for model selection param
    tip = 'Select change point detection model. L1: detects changes in the median of a signal. L2: Cost function detects mean-shifts in a signal. '
    tip += 'RBF: detects changes in the mean of the embedded signal. Recommend RBF.'
    w_model = build_dropdown(desc='Model: ', opts=[('Least absolute deviation (l1)', 'l1'), ('Least squared deviation (l2)', 'l2'), ('Kernelized mean change (rbf)', 'rbf')],
                                val='rbf', handler=cpd_model_handler, params=params, tip=tip)
    
    # int slider for min size between rbeaks value param
    tip = 'Minimum distance between change points (e.g. minimum number of dates between detected breaks). Default is 3.'
    w_min_size = build_int_slider(desc='Min Break Distance:', min_val=1, max_val=25, step_val=1, val=2, 
                                  handler=cpd_min_size_handler, form='d', params=params, tip=tip) 
    
    # int slider for jump value param
    tip = 'Consider only change point locations that are multiple of this particular value. Default is 1.'
    w_jump = build_int_slider(desc='Jump Value:', min_val=1, max_val=25, step_val=1, val=1, 
                                  handler=cpd_jump_handler, form='d', params=params, tip=tip)

    # int slider for pelt penelty param
    tip = 'Set the penelty of break points for pelt change point detection method. Default is 3 for penelty.'
    w_pelt_penelty = build_int_slider(desc='Break Penelty (Pelt):', min_val=1, max_val=25, step_val=1, val=3, 
                                      handler=cpd_pelt_penelty_handler, form='d', params=params, tip=tip)
    
    # int slider for dynamic programming num breaks param
    tip = 'Set the number of break points for dynamic programming change point detection method. Default is 4 breakpoints.'
    w_dynp_n_bkps = build_int_slider(desc='Num Breaks (DynPro):', min_val=2, max_val=24, step_val=2, val=2, 
                                  handler=cpd_dynp_n_bkps_handler, form='d', params=params, tip=tip)
    
    # button for cpd processing
    w_trend_cpd_button = widgets.Button(description='Show')
    w_trend_cpd_button.on_click(handle_button)

    # stylise dashboard
    control_layout = widgets.Layout(margin='10px 15px 10px 0px')
    empty_cell = build_empty_cell(desc=' ', val=' ')
        
    # build dashboard
    w_dash_r0 = widgets.VBox([widgets.HBox([w_model, w_min_size, w_jump], layout=control_layout), 
                              widgets.HBox([w_dynp_n_bkps, w_pelt_penelty], layout=control_layout),
                              widgets.HBox([empty_cell, empty_cell, w_trend_cpd_button], layout=control_layout)])
    
    # create dash as accordion
    w_dash = widgets.Accordion(children=[w_dash_r0], selected_index=0)
    w_dash.set_title(0, 'Change point detection parameters')

    # return dashboard
    return w_dash