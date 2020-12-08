#!/glade/u/home/hawbecke/local/envs/mmc/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from os import path
import xarray as xr
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.colors as colors
import scipy.stats as stats
from netCDF4 import Dataset
from scipy.stats import pearsonr
import wrf
import sys
from datetime import date

from mmctools.plotting import TaylorDiagram
from mmctools.helper_functions import theta_to_T, lowess_mean, calc_uv, calc_wind, w_s, T_d
from mmctools.wrf.utils import Tower, tsout_seriesReader, write_tslist_file
import matplotlib.gridspec as gridspec
from string import ascii_lowercase
import matplotlib.patches as patches
import skimage.morphology
import matplotlib.colors
#sys.path.append('/Users/hawbecke/Code/Python/')
sys.path.append('/glade/u/home/hawbecke/Code/Python/')
from pydicts.baybreezedict import DetectBayBreeze, spatial_breeze_check
from pydicts.obsdict import get_FINO_obs

wrf_dir    = '/glade/scratch/hawbecke/WRF/ATEC/Chesapeake/20190716to20190801/'
#wrf_dir    = '/Users/hawbecke/Research/Chesapeake/Data/WRF/20190716to20190801/'

restarts   = ['CBB_2019071518', 'CBB_2019071718', 'CBB_2019071918', 'CBB_2019072118', 
              'CBB_2019072318', 'CBB_2019072518', 'CBB_2019072718', 'CBB_2019072918']

wrf_start  = ['2019-07-15 18:00:00','2019-07-17 18:00:00','2019-07-19 18:00:00','2019-07-21 18:00:00',
              '2019-07-23 18:00:00','2019-07-25 18:00:00','2019-07-27 18:00:00','2019-07-29 18:00:00',]


from CBB_case_dict import case_dict
cases = list(case_dict.keys())

ncases   = len(cases)
case_dom = [3]*ncases



wrfinput = xr.open_dataset('{}{}/{}/wrfinput_d0{}'.format(wrf_dir,cases[0],restarts[0],case_dom[0]))
land_mask  = np.squeeze(wrfinput.LANDMASK)
water_mask = land_mask.copy().where(land_mask==0.0) + 1.0
lat  = wrfinput.XLAT
lon  = wrfinput.XLONG

nx = len(wrfinput.west_east)
ny = len(wrfinput.south_north)
dx = wrfinput.DX/1000.0
dy = wrfinput.DY/1000.0
x = np.arange(0,nx)*dx
y = np.arange(0,nx)*dy
xy,yx = np.meshgrid(x,y)


onshore_min_max_path = '{}onshore_min_max.nc'.format(wrf_dir)


# Maxiumum radius to look for water:
max_water_dist = 80.0
low_pct_0 = 95.0
upr_pct_0 = 5.0
max_deg_range = 180.0
show_plots = False

window_len = int(max_water_dist/dx)*2
half_window_len = int(window_len/2)
window_center = int((window_len)/2)
window_dist = ((xy[:window_len+1,:window_len+1] - xy[window_center,window_center])**2 + 
               (yx[:window_len+1,:window_len+1] - yx[window_center,window_center])**2)**0.5
window_dist[np.where(window_dist > max_water_dist)] = np.nan
window_filter = window_dist / window_dist

window_x,window_y = np.meshgrid(np.arange(0,np.shape(window_dist)[1]+1)*dy - max_water_dist - 1.5,
                                np.arange(0,np.shape(window_dist)[0]+1)*dx - max_water_dist - 1.5)

window_deg = -1*(180.0*np.arctan(((yx[:window_len+1,:window_len+1] - yx[window_center,window_center])/
                     (xy[:window_len+1,:window_len+1] - xy[window_center,window_center])))/(np.pi) - 90.0)
window_deg[:,:half_window_len] = window_deg[:,:half_window_len] + 180.0
window_deg[np.where(np.isnan(window_dist))] = np.nan

math_deg = 180.0*np.arctan(((yx[:window_len+1,:window_len+1] - yx[window_center,window_center])/
                     (xy[:window_len+1,:window_len+1] - xy[window_center,window_center])))/(np.pi)
math_deg[:,:half_window_len] = math_deg[:,:half_window_len] + 180.0
math_deg[:half_window_len,half_window_len:] = 360 + math_deg[:half_window_len,half_window_len:]
math_deg[np.where(np.isnan(window_dist))] = np.nan


if path.exists(onshore_min_max_path):
    print('loading in onshore min/max dataset!')
    onshore_min_max_ds = xr.open_dataset(onshore_min_max_path)
    onshore_min = onshore_min_max_ds.onshore_min
    onshore_max = onshore_min_max_ds.onshore_max
else:
    onshore_min = np.zeros((ny,nx))*np.nan
    onshore_max = np.zeros((ny,nx))*np.nan
    for ii in np.arange(half_window_len,nx-half_window_len):
        for jj in np.arange(half_window_len,ny-half_window_len): 
            if land_mask[jj,ii] == 1.0:
                loc_water_mask = water_mask[jj-half_window_len:jj+half_window_len+1, ii-half_window_len:ii+half_window_len+1]
                dist_water = loc_water_mask * window_dist
                deg_water  = loc_water_mask * window_deg

                # Break down the water bodies into groups:
                water_bodies = skimage.morphology.label(~np.isnan(deg_water)).astype(np.float32)
                water_bodies[water_bodies==0.0] = np.nan

                water_body_size = {}
                water_body_dist = {}
                min_water_distance = 999.9
                closest_water_body = 0.0

                # If we have water bodies to check, enter loop:
                if ~np.all(np.isnan(water_bodies)): 
                    for i in np.arange(1.0,np.nanmax(water_bodies)+1.0): # Loop over all identified water bodies
                        water_size = len(water_bodies[water_bodies==i])
                        if water_size < 8: # Only check for large bodies
                            water_bodies[water_bodies==i] = np.nan
                        else:
                            water_body_size[i] = water_size
                            water_body = water_bodies.copy()
                            water_body[water_bodies!=i] = np.nan # Check only this water body
                            water_body[~np.isnan(water_body)] = 1.0 # Set values to 1 to get mask for distance calculation
                            water_body_min_dist = np.nanpercentile(water_body*dist_water,50)
                            water_body_dist[i] = water_body_min_dist

                    # Small water bodies were removed, check to see if there are any large ones:
                    if ~np.all(np.isnan(water_bodies)):
                        # Find the largest and closest water bodies:
                        largest_water_body = max(water_body_size,key=water_body_size.get)
                        water_body_id = largest_water_body
                        # Loop over all other water bodies
                        for i in water_body_size.keys():
                            if i != water_body_id:
                                # Check to see if this water body is still relatively large:
                                if water_body_size[i] >= 0.5*water_body_size[water_body_id]:
                                    # Assign this the same water body ID
                                    water_bodies[water_bodies==i] = water_body_id

                        # Set the selected water body (bodies) to 1.0 for masking
                        water_bodies[water_bodies==water_body_id] = 1.0

                        # Multiply water body mask by direction to water:
                        deg_water *= water_bodies

                        # Check to see if there are negative and positive values in the same water body:
                        deg_range = float(np.nanmax(deg_water)) - float(np.nanmin(deg_water))
                        if deg_range > 300:
                            deg_water[np.where(deg_water>300)] -= 360.0

                        # Set limits for upper and lower bounds.
                        # If range is too big (> max_deg_range) then we iterate by 5 degrees
                        # ... on the upper and lower limits until the range is sufficient.
                        if np.nanmax(water_bodies) > 0:
                            good_lims = False
                            low_pct = low_pct_0
                            upr_pct = upr_pct_0
                            while good_lims == False:
                                lowr_lim = np.nanpercentile(deg_water,low_pct)
                                uppr_lim = np.nanpercentile(deg_water,upr_pct)
                                if lowr_lim - uppr_lim < max_deg_range:
                                    good_lims = True
                                else:
                                    low_pct -= 5.0
                                    upr_pct += 5.0                                
                        else: # Set limits to nan when water_bodies is all nan
                            lowr_lim = np.nan
                            uppr_lim = np.nan

                        onshore_min[jj,ii] = uppr_lim 
                        onshore_max[jj,ii] = lowr_lim

    onshore_min_da = land_mask.copy()
    onshore_min_da.data = onshore_min
    onshore_max_da = land_mask.copy()
    onshore_max_da.data = onshore_max
    onshore_min_max_ds = xr.Dataset({'onshore_min':onshore_min_da,
                                     'onshore_max':onshore_max_da})
    onshore_min_max_ds.to_netcdf(onshore_min_max_path)


for cc,case in enumerate(cases[8:9]):
    print(case)
    for rr,rst in enumerate(restarts):
        print(rst)
        breeze_ds_path = '{}{}/{}/spatial_breeze_detection_d0{}.nc'.format(wrf_dir,case,rst,case_dom[cc])
        if path.exists(breeze_ds_path):
            print('Numerical detection dataset already created!')
            print(breeze_ds_path)
            breeze_ds_f = xr.open_dataset(breeze_ds_path)
        else:
            wrfout_files = sorted(glob.glob('{}{}/{}/wrfout_d0{}*'.format(wrf_dir,case,rst,case_dom[cc])))
            init_vars = True
            for ww,wrfout_file in enumerate(wrfout_files[12:]):
                out_time = pd.to_datetime(wrfout_file.split('/')[-1][11:],format='%Y-%m-%d_%H:%M:%S')
                if (out_time.hour >=13) & (out_time.hour <=22):
                    wrfout = xr.open_dataset(wrfout_file)
                    wrfout = np.squeeze(wrfout)
                    time_str = date.strftime(pd.to_datetime(wrfout.XTIME.data),'%Y-%m-%d %H:%M')
                    print(time_str)
                    if init_vars:
                        wrfinput = xr.open_dataset('{}{}/{}/wrfinput_d0{}'.format(wrf_dir,cases[6],restarts[0],case_dom[6]))
                        land_mask  = np.squeeze(wrfinput.LANDMASK)
                        hgt        = np.squeeze(wrfinput.HGT)
                        water_mask = land_mask.copy().where(land_mask==0.0) + 1.0
                        lat  = wrfout.XLAT
                        lon  = wrfout.XLONG
                        z_f = (np.squeeze(wrfout.PH) + np.squeeze(wrfout.PHB))/9.8 - np.squeeze(wrfout.HGT)
                        zs_f = 0.5*(z_f[1:,:,:]+z_f[:-1,:,:])

                    breeze_ds = spatial_breeze_check(onshore_min,
                                                     onshore_max,
                                                     wrfout,
                                                     land_mask=land_mask,
                                                     wdir_check='vertical',
                                                     dT_cutoff_pct=50.0
                                                    )

                    if init_vars:
                        breeze_ds_f = breeze_ds.copy()
                        init_vars = False
                    else:
                        breeze_ds_f = breeze_ds.combine_first(breeze_ds_f) 

            breeze_ds_f['onshore_min'] = onshore_min
            breeze_ds_f['onshore_max'] = onshore_max
            breeze_ds_f['land_mask'] = land_mask

            breeze_ds_f.attrs = {'WRF_DIR':'/'.join(wrfout_files[0].split('/')[:-1]),
                                  'domain':wrfout_files[0].split('/')[-1].split('_')[1]}

            breeze_ds_f.to_netcdf(breeze_ds_path)
    wefwef
