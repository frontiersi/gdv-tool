# scripts for creation of jupyter interactive maps and functionality
import os
import functools
import numpy as np
import pandas as pd
import ipyleaflet as leaf
import ipywidgets as widgets
from gdv_scripts import helpers
from gdv_scripts import data

## EVENT HANDLERS ##
# event handler for study area selection
def clicked(params, study_area_map, **kwargs):
    try:
        if kwargs['event'] == 'click':
            study_area_name = kwargs['properties']['name']

            for lyr in study_area_map.layers:
                if isinstance(lyr, leaf.leaflet.GeoJSON):
                    if lyr.data['properties']['name'] == study_area_name:

                        # change style to show selected
                        lyr.style = {'color': 'blue'}

                        # get top left and bottom right coord lon, lat, 
                        tl = lyr.data['geometry']['coordinates'][0][1]
                        br = lyr.data['geometry']['coordinates'][0][3]

                        # add to storage
                        params.gcs_nw = (tl[0], tl[1])
                        params.gcs_se = (br[0], br[1])
                        
                        # transform to albers also
                        params.alb_nw = (helpers.transform_point(tl[0], tl[1], 'epsg:4326', 'epsg:3577'))
                        params.alb_se = (helpers.transform_point(br[0], br[1], 'epsg:4326', 'epsg:3577'))
                    else:
                        lyr.style = {'color': 'red'}
    except Exception as e:
        print('Error occurred during of map click event of type {0}. Stopping.'.format(e))
        raise e

# event handler for study area drawing
def drawing_handler(self, action, geo_json, params, study_area_map):   
    try:
        # get geom and coordinates
        geom = geo_json['geometry']
        coords = geom['coordinates']
        
        # get tl and br coords
        tl = coords[0][1]
        br = coords[0][3] 
        
        # add to storage param
        params.gcs_nw = (tl[0], tl[1])
        params.gcs_se = (br[0], br[1])
        
        # transform to albers also
        params.alb_nw = (helpers.transform_point(tl[0], tl[1], 'epsg:4326', 'epsg:3577'))
        params.alb_se = (helpers.transform_point(br[0], br[1], 'epsg:4326', 'epsg:3577'))
        
        
        # clear existing selections
        for lyr in study_area_map.layers:
            if isinstance(lyr, leaf.leaflet.GeoJSON):
                lyr.style = {'color': 'red'}
                
        # clear drawing
        self.clear()
    except Exception as e:
        print('Error occurred during of map drawing event of type {0}. Stopping.'.format(e))
        raise e


## MAP AND OBJECT TYPE CONSTRUCTORS ##
# build basic template map for further data additions
def generate_basic_map(study_center, w_full_screen=True, w_measure=True, l_width='100%', l_height='500px'):
    
    # generate and rename esri and otm basemaps
    def create_basemaps():
        try:
            esri_basemap = leaf.basemap_to_tiles(leaf.basemaps.Esri.WorldImagery)
            otm_basemap = leaf.basemap_to_tiles(leaf.basemaps.OpenTopoMap)
            esri_basemap.name, otm_basemap.name = 'ESRI Satellite', 'OTM'   # rename layers

            return esri_basemap, otm_basemap
        except:
            print('Problem during creation of leaflet basemaps. Please check.')
            raise

    # basic check
    if not all(study_center):
        print('No study center. Did you run the above cells?')
        return None
        
    # begin creation of basic map
    try:
        # get esri and otm basemaps
        bm_esri, bm_otm = create_basemaps()

        # build map
        leaf_map = leaf.Map(center=study_center, zoom=9, box_zoom=True, zoom_control=False,
                            layers=(bm_otm, bm_esri), layout={'width': l_width, 'height': l_height})

        # add layers control panel
        layers_control = leaf.LayersControl(position='topleft')
        leaf_map.add_control(layers_control)

        # add full screen widget if requested
        if w_full_screen:
            leaf_map.add_control(leaf.FullScreenControl())

        # add measure widget if requested
        if w_measure:
            leaf_map.add_control(leaf.MeasureControl(primary_length_unit='meters',
                                                     secondary_length_unit='kilometers'))

        # for cleaner ui - show zoom control after layers control panel
        leaf_map.zoom_control = True

        return leaf_map
    except Exception as e:
        print('Error occurred during create_basemaps of type {0}. Stopping.'.format(e))
        raise e      

# build geojson study area boundaries for demo
def build_study_areas_geojson_list():
    
    # set default style
    style = {'stroke': True, 'color': 'red', 'weight': 4, 'opacity': 0.5, 'fill': True, 'fillColor': None, 'fillOpacity': 0.2, 'clickable': True}

    # add raw geojson into list
    json = [{'type': 'Feature', 
             'properties': {'name': 'yandi'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[118.959804, -22.810513], [118.959804, -22.670243], [119.180606, -22.670243], [119.180606, -22.810513], [118.959804, -22.810513]]]
             }
            }, 
            {'type': 'Feature', 
             'properties': {'name': 'royhill'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.870154, -22.689841], [119.870154, -22.419704], [120.151737, -22.419704], [120.151737, -22.689841], [119.870154, -22.689841]]]
             }
            }, 
            {'type': 'Feature', 
             'properties': {'name': 'ophthalmia'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.800822, -23.398774], [119.800822, -23.160984], [119.917576, -23.160984], [119.917576, -23.398774], [119.800822, -23.398774]]]
                         }
            },

            {'type': 'Feature', 
             'properties': {'name': 'royhill_sub_1'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[120.018101, -22.611297], [120.018101, -22.547896], [120.133426, -22.547896], [120.133426, -22.611297], [120.018101, -22.611297]]]
                         }
            },
            {'type': 'Feature', 
             'properties': {'name': 'royhill_sub_2'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[120.121012, -22.521259], [120.121012, -22.428623], [120.2027, -22.428623], [120.2027, -22.521259], [120.121012, -22.521259]]]
                         }
            },
            {'type': 'Feature', 
             'properties': {'name': 'yandi_sub_1'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.097814, -22.798882], [119.097814, -22.715934], [119.204901, -22.715934], [119.204901, -22.798882], [119.097814, -22.798882]]]
                         }
            },
            {'type': 'Feature', 
             'properties': {'name': 'yandi_sub_2'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.248854, -22.819137], [119.248854, -22.736201], [119.342898, -22.736201], [119.342898, -22.819137], [119.248854, -22.819137]]]
                         }
            },
            {'type': 'Feature', 
             'properties': {'name': 'ophthalmia_sub_1'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.789221, -23.326405], [119.789221, -23.272485], [119.929945, -23.272485], [119.929945, -23.326405], [119.789221, -23.326405]]]
                         }
            },            
            {'type': 'Feature', 
             'properties': {'name': 'ophthalmia_sub_2'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.884639, -23.232109], [119.884639, -23.132694], [119.951568, -23.132694], [119.951568, -23.232109], [119.884639, -23.232109]]]
                         }
            },
            {'type': 'Feature', 
             'properties': {'name': 'extra_sub_1'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.177669, -22.939019], [119.177669, -22.858688], [119.261073, -22.858688], [119.261073, -22.939019], [119.177669, -22.939019]]]
                         }
            },            
            {'type': 'Feature', 
             'properties': {'name': 'extra_sub_2'}, 
             'geometry': {'type': 'Polygon', 
                          'coordinates': [[[119.499759, -22.729779], [119.499759, -22.657562], [119.592431, -22.657562], [119.592431, -22.729779], [119.499759, -22.729779]]]
                         }
            }           
           ] 
    
    # generate leaf geojson layers, add on_click event
    gjson = [leaf.GeoJSON(data=j, style=style) for j in json]

    return gjson
  

## BUILD SPECIFIC INTERACTIVE MAPS ##
# initialise study area interactive map
def load_sat_fetch_map(params, study_center):
    
    # basic checks
    if not params:
        print('Parameter variable is empty. Did you miss any previous cells?')
        return None
    elif not all(study_center):
        print('No study center provided. Did you miss any previous cells?')
        return None
        
    try:
        # load basic map
        study_area_map = generate_basic_map(study_center=study_center)
                  
        # generate study area locations for workshops
        study_area_gjson = build_study_areas_geojson_list()
        
        # add handler to each gjson and add to map
        if not study_area_map:
            print('No study area could be generated. Aborting map creation.')
            return None
            
        for lyr in study_area_gjson:
            lyr.on_click(functools.partial(clicked, params=params, study_area_map=study_area_map))
            study_area_map.add_layer(lyr)
               
        # add interactive draw control
        draw_control = leaf.DrawControl(rectangle={'shapeOptions': {'color': '#FFF700'}}, polyline={}, 
                                        polygon={}, circlemarker={})
        
        # add draw handler
        draw_control.on_draw(functools.partial(drawing_handler, params=params, study_area_map=study_area_map))
        study_area_map.add_control(draw_control)
        
        return study_area_map
    except Exception as e:
        print('Error occurred during load_sat_fetch_map of type {0}. Stopping.'.format(e))
        raise e      

# initialise likelihood interactive map
def load_likelihood_map(da_raw, da_thresh, params):
    
    # basic checks
    if not da_raw.any() or not da_thresh.any():
        print('Dataset variable is empty. Did you miss any previous cells?')
        return None
    elif not params:
        print('Parameter variable is empty. Did you miss any previous cells?')
        return None
        
    try:
        # setup output
        out_folder = os.getcwd()
        mask_value = 0

        # get study center, load basic map, update its zoom
        study_center = (params.gcs_nw[1] + params.gcs_se[1]) / 2, (params.gcs_nw[0] + params.gcs_se[0]) / 2
        likelihood_map = generate_basic_map(study_center=study_center)
        likelihood_map.zoom = 12

        # fix zero values, reproject tif and save as png
        da_raw = da_raw.where(da_raw != 0, 0.0001)
        like_img, like_bounds = data.create_and_reproject_image_from_array(da=da_raw, out_folder=out_folder, filename='like_raw', filetype='.tif', mask_value=mask_value)
        thresh_img, thresh_bounds = data.create_and_reproject_image_from_array(da=da_thresh, out_folder=out_folder, filename='like_thresh', filetype='.tif', mask_value=mask_value)
        
        # convert png to url
        like_url = data.generate_png_url(in_img=like_img, color_map='Spectral')
        thresh_url = data.generate_png_url(in_img=thresh_img, color_map='Spectral')
        
        # add images to map they exist
        if like_url and like_bounds:
                n, l, b, r = like_bounds.top, like_bounds.left, like_bounds.bottom, like_bounds.right
                io_like = leaf.ImageOverlay(name='Likeihood Full', url=like_url, bounds=((n, l), (b, r)))
                likelihood_map.add_layer(io_like)      
        if thresh_url and thresh_bounds:
                n, l, b, r = thresh_bounds.top, thresh_bounds.left, thresh_bounds.bottom, thresh_bounds.right
                io_thresh = leaf.ImageOverlay(name='Likeihood Thresholded', url=thresh_url, bounds=((n, l), (b, r)))
                likelihood_map.add_layer(io_thresh)  
        
        return likelihood_map
    
    except Exception as e:
        print('Error occurred during load_likelihood_map of type {0}. Stopping.'.format(e))
        raise e     

# initialise trend interactive map   
def load_trend_map(ds, mk_result, params):
        
    # basic checks
    if not ds or not mk_result:
        print('Dataset is empty. Did you miss any previous cells?')
        return None
    
    if not params:
        print('Parameter variable is empty. Did you miss any previous cells?')
        return None
    
    try:
        # setup output
        out_folder = os.getcwd()
        mask_value = 0

        # get study center, load basic map, update its zoom
        study_center = (params.gcs_nw[1] + params.gcs_se[1]) / 2, (params.gcs_nw[0] + params.gcs_se[0]) / 2
        trend_map = generate_basic_map(study_center=study_center)
        trend_map.zoom = 12

        # fix zero values, reproject tif and save as png
        mk_result['p'] = mk_result['p'].where(mk_result['p'] != 1.0, np.nan)
        
        # change tau color map depending on trend type
        if params.trend_type == 'increasing':
            tau_cmap = 'YlGn'
        elif params.trend_type == 'decreasing':
            tau_cmap = 'YlOrRd_r'
        else:
            tau_cmap = 'RdYlGn'
        
        # create reprojected tifs, bounds
        p_img, p_bounds = data.create_and_reproject_image_from_array(da=mk_result['p'], out_folder=out_folder, filename='p', filetype='.tif', mask_value=mask_value)
        tau_img, tau_bounds = data.create_and_reproject_image_from_array(da=mk_result['tau'], out_folder=out_folder, filename='tau', filetype='.tif', mask_value=mask_value)
        z_img, z_bounds = data.create_and_reproject_image_from_array(da=mk_result['z'], out_folder=out_folder, filename='z', filetype='.tif', mask_value=mask_value)
        slope_img, slope_bounds = data.create_and_reproject_image_from_array(da=mk_result['slope'], out_folder=out_folder, filename='slope', filetype='.tif', mask_value=mask_value)
        s_img, s_bounds = data.create_and_reproject_image_from_array(da=mk_result['s'], out_folder=out_folder, filename='s', filetype='.tif', mask_value=mask_value) 

        # convert png to url
        p_url = data.generate_png_url(in_img=p_img, color_map='hot_r')
        tau_url = data.generate_png_url(in_img=tau_img, color_map=tau_cmap)
        z_url = data.generate_png_url(in_img=z_img, color_map='Spectral')
        slope_url = data.generate_png_url(in_img=slope_img, color_map='PiYG')
        s_url = data.generate_png_url(in_img=s_img, color_map='Spectral')
        
        # add images to map they exist
  
        if z_url and z_bounds:
            n, l, b, r = z_bounds.top, z_bounds.left, z_bounds.bottom, z_bounds.right
            io_z = leaf.ImageOverlay(name='Normalised Test Statistics (Z)', url=z_url, bounds=((n, l), (b, r)))
            #trend_map.add_layer(io_z)      
        if slope_url and slope_bounds:
            n, l, b, r = slope_bounds.top, slope_bounds.left, slope_bounds.bottom, slope_bounds.right
            io_slope = leaf.ImageOverlay(name='Sens Slope', url=slope_url, bounds=((n, l), (b, r)))
            trend_map.add_layer(io_slope)
        if s_url and s_bounds:
            n, l, b, r = s_bounds.top, s_bounds.left, s_bounds.bottom, s_bounds.right
            io_s = leaf.ImageOverlay(name='Mann-Kendall Score', url=s_url, bounds=((n, l), (b, r)))
            #trend_map.add_layer(io_s)
        if p_url and p_bounds:
            n, l, b, r = p_bounds.top, p_bounds.left, p_bounds.bottom, p_bounds.right
            io_p = leaf.ImageOverlay(name='P-Value', url=p_url, bounds=((n, l), (b, r)))
            trend_map.add_layer(io_p)    
        if tau_url and tau_bounds:
            n, l, b, r = tau_bounds.top, tau_bounds.left, tau_bounds.bottom, tau_bounds.right
            io_tau = leaf.ImageOverlay(name='Tau', url=tau_url, bounds=((n, l), (b, r)))
            trend_map.add_layer(io_tau)  
        
        return trend_map
    
    except Exception as e:
        print('Error occurred during load_trend_map of type {0}. Stopping.'.format(e))
        raise e   

# initialise change interactive map   
def load_change_map(ds, cva_ds, params):
        
    # basic checks
    if not ds or not cva_ds:
        print('Dataset is empty. Did you miss any previous cells?')
        return None
    
    if not params:
        print('Parameter variable is empty. Did you miss any previous cells?')
        return None
    
    # slider img change handler
    def change_handler(control):
        if control['type'] == 'change':
            try:
                # get old and new slider info, and data for new from dict
                old_cva, new_cva = control['old'], control['new']
                
                ang_url, ang_bounds = cva_dict[new_cva][0], cva_dict[new_cva][2]
                mag_url, mag_bounds = cva_dict[new_cva][1], cva_dict[new_cva][3]

                # remove old layer
                [change_map.remove_layer(lyr) for lyr in change_map.layers if lyr.name == 'Angles' or lyr.name == 'Magnitudes']
                
                # get new ang image data
                n, l, b, r = ang_bounds.top, ang_bounds.left, ang_bounds.bottom, ang_bounds.right
                io_ang = leaf.ImageOverlay(name='Angles', url=ang_url, bounds=((n, l), (b, r)))
                change_map.add_layer(io_ang)
                
                # get new mag image data
                n, l, b, r = mag_bounds.top, mag_bounds.left, mag_bounds.bottom, mag_bounds.right
                io_mag = leaf.ImageOverlay(name='Magnitudes', url=mag_url, bounds=((n, l), (b, r)))
                change_map.add_layer(io_mag)
              
            except:
                print('Problem during change slider functionality. Skipping.')
    
    try:
        # setup output
        out_folder = os.getcwd()
        mask_value = 0

        # get study center, load basic map, update its zoom
        study_center = (params.gcs_nw[1] + params.gcs_se[1]) / 2, (params.gcs_nw[0] + params.gcs_se[0]) / 2
        change_map = generate_basic_map(study_center=study_center)
        change_map.zoom = 12
        
        # get color map depending on whether angle is inc or dec
        if params.change_angles_range[0] >= 0 and params.change_angles_range[0] <= 180:
            cmap = 'YlOrRd'
        else:
            cmap = 'YlGnBu'
        
        # build dict of mag img and bound info
        cva_dict = {}
        for idx in cva_ds['time']:
            scene = cva_ds.sel(time=idx)
            
            # get date for dict key
            dt = pd.to_datetime(str(scene['time'].values))
            dt = str(dt.year) + '-' + str(dt.month)  
            
            # create reprojected tifs, bounds
            ang_img, ang_bounds = data.create_and_reproject_image_from_array(da=scene['angout'], out_folder=out_folder, filename='angout', 
                                                                             filetype='.tif', mask_value=mask_value)
                
            mag_img, mag_bounds = data.create_and_reproject_image_from_array(da=scene['magout'], out_folder=out_folder, filename='magout', 
                                                                             filetype='.tif', mask_value=mask_value)
            
            # convert png to url
            ang_url = data.generate_png_url(in_img=ang_img, color_map='gist_rainbow')
            mag_url = data.generate_png_url(in_img=mag_img, color_map=cmap)
            
            # add to dict
            cva_dict[dt] = [ang_url, mag_url, ang_bounds, mag_bounds]
            
        # get lsit of dict keys (dates)
        cva_dates = list(cva_dict)
      
        # add slider and observer handler
        cva_slider = widgets.SelectionSlider(description='Date: ', options=cva_dates, value=cva_dates[0]) 
        cva_slider.observe(change_handler, names='value')

        # creat map widget and add
        w_cva_date_widget = leaf.WidgetControl(widget=cva_slider, position='topright')
        change_map.add_control(w_cva_date_widget)

        return change_map

    except Exception as e:
        print('Error occurred during load_change_map of type {0}. Stopping.'.format(e))
        raise e   
    