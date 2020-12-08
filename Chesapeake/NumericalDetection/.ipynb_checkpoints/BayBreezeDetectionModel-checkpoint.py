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
from pydicts.baybreezedict import DetectBayBreeze, spatial_breeze_check, find_onshore_minmax
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


bad_cases  = []
cases_of_interest = cases.copy()
for cc,case in enumerate(cases_of_interest):
    numerical_detection_path = '{}/spatial_breeze_detection_CLDnRAIN_d0{}.nc'.format('{}{}/'.format(wrf_dir, case),case_dom[cc])
    if path.exists(numerical_detection_path):
        print('already did ',case)
    else:
    

        init_vars = True
        wrfout_files = sorted(glob.glob('{}{}/CBB_2019*/wrfout_d0{}*'.format(wrf_dir,case,case_dom[0])))
        if len(wrfout_files) <= 1:
            cases_of_interest.remove(case)

            bad_cases.append(case)
            
cases = ['ERA5_YSU_CHRN_ER5_NOSK_3DOM','ERAI_YSU_CHRN_OST_NOSK_4DOM','GFSR_YSU_CHRN_GFS_NOSK_3DOM']
cases_of_interest = cases.copy()


wrfout_files = sorted(glob.glob('{}{}/{}/wrfout_d0{}*'.format(wrf_dir,cases_of_interest[0],restarts[0],case_dom[0])))

time_of_interest = '2019-07-16 19:00:00'
#time_of_interest = '2019-07-16 13:00:00'

for wrfout_file in wrfout_files:
    if time_of_interest.replace(' ','_') in wrfout_file:
        file_of_interest = wrfout_file

wrfout = xr.open_dataset(file_of_interest)
wrfout = np.squeeze(wrfout)
wrfinput = xr.open_dataset('{}{}/{}/wrfinput_d0{}'.format(wrf_dir,cases_of_interest[0],restarts[0],case_dom[0]))

land_mask  = np.squeeze(wrfinput.LANDMASK)
hgt        = np.squeeze(wrfinput.HGT)
water_mask = land_mask.copy().where(land_mask==0.0) + 1.0
lat  = wrfout.XLAT
lon  = wrfout.XLONG

nx = len(wrfout.west_east)
ny = len(wrfout.south_north)
dx = wrfout.DX/1000.0
dy = wrfout.DY/1000.0
x = np.arange(0,nx)*dx
y = np.arange(0,nx)*dy
xy,yx = np.meshgrid(x,y)

t2   = wrfout.T2.where(land_mask == 1.0)
u10  = wrfout.U10.where(land_mask == 1.0)
v10  = wrfout.V10.where(land_mask == 1.0)
sfcP = wrfout.PSFC.where(land_mask == 1.0)
temp = np.squeeze(wrfout.T)
z_f = (np.squeeze(wrfout.PH) + np.squeeze(wrfout.PHB))/9.8 - np.squeeze(wrfout.HGT)
zs_f = 0.5*(z_f[1:,:,:]+z_f[:-1,:,:])

vel10,dir10 = calc_wind(wrfout,u='U10',v='V10')
vel10 = vel10.where(land_mask == 1.0)
dir10 = dir10.where(land_mask == 1.0)

top_ind = 18
bot_ind = 0

u = wrfout.U[top_ind,:,:].data
v = wrfout.V[top_ind,:,:].data

u = 0.5*(u[:,1:] + u[:,:-1])
v = 0.5*(v[1:,:] + v[:-1,:])
wdir1km = 180. + np.degrees(np.arctan2(u, v))

obs_dir = '/glade/u/home/hawbecke/Research/Chesapeake/Data/Obs/'
awos_ds = xr.open_dataset('{}AWOS_2019to2020.nc'.format(obs_dir))
asos_ds = xr.open_dataset('{}ASOS_2019to2020.nc'.format(obs_dir))
apg_ds  = xr.open_dataset('{}APG_data_2019.nc'.format(obs_dir))

#obs_dir = '/Users/hawbecke/Research/Chesapeake/Data/Obs/'
#awos_ds = xr.open_dataset('{}AWOS/AWOS_2019to2020.nc'.format(obs_dir))
#asos_ds = xr.open_dataset('{}ASOS/ASOS_2019to2020.nc'.format(obs_dir))
## Add APG data into AWOS
#apg_ds  = xr.open_dataset('{}APG/APG_data_2019.nc'.format(obs_dir))

apg_stn_list = list(apg_ds.station.data)
for ss,apg_stn in enumerate(apg_stn_list):
    if apg_stn == 'PAA':
        apg_stn_list[ss] = 'APG'
apg_ds = apg_ds.assign_coords({'station': apg_stn_list})
non_apg_list = apg_stn_list.copy()
non_apg_list.remove('APG')
apg_ds = apg_ds.drop_sel(station=non_apg_list)
apg_ds = apg_ds.drop(['rh','gust','radt','rain','alt'])

temp = awos_ds.sel(station='APG').combine_first(apg_ds)
temp = temp.squeeze()
temp = temp.expand_dims('station')

lat_da = xr.DataArray([float(apg_ds.lat)], dims=('station'))
lon_da = xr.DataArray([float(apg_ds.lon)], dims=('station'))

temp = temp.assign_coords({'lat':lat_da})
temp = temp.assign_coords({'lon':lon_da})

awos_ds = awos_ds.drop_sel(station=['APG'])

awos_ds = xr.merge([awos_ds,temp])

near_shore_stations = ['APG', 'BWI', 'MTN', 'NAK', 'NHK', 'NUI']

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
new_cmap = truncate_colormap(plt.cm.terrain, 0.3, 0.8)
new_cmap.set_bad(color='darkblue')

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return colors.LinearSegmentedColormap('colormap',cdict,1024)

light_cmap = cmap_map(lambda x: x/2 + 0.5, new_cmap)
light_cmap.set_bad(color='lightsteelblue')

#met_deg = 270.0
#expected_deg = 180.0
#print('Met:\t\t{}\nExpected:\t{}\ngot:\t\t{}'.format(met_deg,expected_deg,convert_met_to_math(met_deg)))



onshore_min_max_path = '{}onshore_min_max_new.nc'.format(wrf_dir)
if path.exists(onshore_min_max_path):
    print('loading in onshore min/max dataset!')
    onshore_min_max_ds = xr.open_dataset(onshore_min_max_path)
    onshore_min = onshore_min_max_ds.onshore_min
    onshore_max = onshore_min_max_ds.onshore_max
else:
    onshore_min_max_ds = find_onshore_minmax(land_mask,
                                             max_water_dist = 110.0,
                                             low_pct_0 = 90.0,
                                             upr_pct_0 = 10.0,
                                             max_deg_range = 350.0,
                                             min_water_size=20)
    
    onshore_min = onshore_min_max_ds.onshore_min
    onshore_max = onshore_min_max_ds.onshore_max
    onshore_min_max_ds.to_netcdf(onshore_min_max_path)
    
    
    
def convert_breeze_dict_to_xr(breeze_dict):
    for kk,key in enumerate(breeze_dict.keys()):
        if key == 'breeze':
            var = breeze_dict['good_wdir'].copy()
            var.data = breeze_dict[key]
        else:
            var = breeze_dict[key]
            var.name = key
        if kk == 0:
            ds = xr.Dataset({key: var})
        else:
            ds = xr.merge([ds,var])
            
    ds = ds.expand_dims('datetime')
    dtime = ds.XTIME.expand_dims('datetime')
    ds = ds.drop('XTIME')
    ds['datetime'] = dtime
    return(ds)



breeze_cmap = plt.cm.RdYlGn(np.linspace(0.1,0.9,3))
#breeze_cmap = [[0.0, 0.0, 0.0, 1.0]]
#for ll,lvl in enumerate(np.arange(0,2)):
#    breeze_cmap = np.append(breeze_cmap,[list(breeze_cmap_f[ll])],axis=0)
#breeze_cmap = np.append(breeze_cmap,[[1.0, 0.0, 0.0, 1.0]],axis=0)
#for ll,lvl in enumerate(np.arange(0,3)):
#    plt.plot(np.arange(0,10),np.arange(0,10)*ll,c=breeze_cmap[ll])
    
    
stn_colors = plt.cm.nipy_spectral(np.linspace(0.05,0.95,len(near_shore_stations)))
#for ss,stn in enumerate(near_shore_stations):    
#    plt.plot(np.arange(0,10),np.arange(0,10)*ss,c=stn_colors[ss])

stn_dict = {}
for ss,stn in enumerate(near_shore_stations):
    stn_dict[stn] = {}
    #if stn == 'APG':
    #    stn_ds = apg_ds.sel(station=stn)
    if stn in awos_ds.station:
        stn_ds = awos_ds.sel(station=stn)
    elif stn in asos_ds.station:
        stn_ds = asos_ds.sel(station=stn)
    stn_lat = stn_ds.lat.data
    stn_lon = stn_ds.lon.data
    dist_from_stn = ((lon-stn_lon)**2 + (lat-stn_lat)**2)**0.5
    stn_j,stn_i = np.where(dist_from_stn == np.nanmin(dist_from_stn))
    stn_j = int(stn_j)
    stn_i = int(stn_i)
    stn_dict[stn]['lat'] = float(stn_ds.lat.data)
    stn_dict[stn]['lon'] = float(stn_ds.lon.data)
    stn_dict[stn]['i'] = stn_i
    stn_dict[stn]['j'] = stn_j

show_plot = False
bad_cases = []


check_clouds = True
check_rain = True
print(cases_of_interest,'PSH')

for cc,case in enumerate(cases_of_interest):
    numerical_detection_path = '{}/spatial_breeze_detection_CLDnRAIN_d0{}.nc'.format('{}{}/'.format(wrf_dir, case),case_dom[0])
    if path.exists(numerical_detection_path):
        print('already did ',case)
    else:
    

        init_vars = True
        wrfout_files = sorted(glob.glob('{}{}/CBB_2019*/wrfout_d0{}*'.format(wrf_dir,case,case_dom[0])))
        if len(wrfout_files) > 1:
            for ww,wrfout_file in enumerate(wrfout_files):
                wrfout = xr.open_dataset(wrfout_file)
                wrfout = np.squeeze(wrfout)
                out_time = pd.to_datetime(wrfout.XTIME.data)
                if (out_time.hour >=13) & (out_time.hour <=22):
                    time_str = date.strftime(pd.to_datetime(wrfout.XTIME.data),'%Y-%m-%d %H:%M')
                    print(time_str)
                    if init_vars:
                        wrfinput = xr.open_dataset('{}{}/{}/wrfinput_d0{}'.format(wrf_dir,cases[1],restarts[0],case_dom[0]))
                        land_mask  = np.squeeze(wrfinput.LANDMASK)
                        hgt        = np.squeeze(wrfinput.HGT)
                        water_mask = land_mask.copy().where(land_mask==0.0) + 1.0
                        lat  = wrfout.XLAT
                        lon  = wrfout.XLONG
                        z_f = (np.squeeze(wrfout.PH) + np.squeeze(wrfout.PHB))/9.8 - np.squeeze(wrfout.HGT)
                        zs_f = 0.5*(z_f[1:,:,:]+z_f[:-1,:,:])


                        rain = wrfout.RAINC + wrfout.RAINNC + wrfout.RAINSH
                        rain_prev = rain.copy()*0.0
                    else:
                        rain = wrfout.RAINC + wrfout.RAINNC + wrfout.RAINSH - rain_prev

                    rain_prev += rain 
                    print(np.nanmin(rain))

                    breeze_ds = spatial_breeze_check(onshore_min,
                                                     onshore_max,
                                                     wrfout,
                                                     land_mask=land_mask,
                                                     wdir_check='vertical',
                                                     dT_cutoff_pct=50.0,
                                                     check_rain=check_rain,
                                                     rain_da=rain,
                                                     check_clouds=check_clouds,
                                                     cloud_cutoff=0.7
                                                    )

                    
                    if init_vars:
                        breeze_ds_f = breeze_ds.copy()
                        init_vars = False
                    else:
                        breeze_ds_f = breeze_ds.combine_first(breeze_ds_f) 



                    if show_plot:
                        u10 = wrfout.U10
                        v10 = wrfout.V10
                        vel10,dir10 = calc_wind(wrfout,u='U10',v='V10')
                        vel10 = vel10.where(land_mask == 1.0)
                        dir10 = dir10.where(land_mask == 1.0)

                        fig = plt.figure(figsize=(16,20))
                        nrow,ncol = 4,4
                        tmp_plt = plt.subplot2grid((nrow,ncol),(0,0),colspan=2,aspect='equal')
                        dT_lvl = np.nanmax(abs(breeze_ds.dT_full.sel(datetime=time_str)))
                        tmp_plt_pc = tmp_plt.pcolormesh(lon,lat,breeze_ds.dT_full.sel(datetime=time_str),
                                                        cmap=plt.cm.RdYlBu,norm=Normalize(-dT_lvl,dT_lvl))
                        tmp_cbar = plt.colorbar(tmp_plt_pc,ax=tmp_plt)
                        tmp_cbar.ax.tick_params(labelsize=13)
                        tmp_plt.contour(lon,lat,land_mask,levels=[0.5],colors='k',alpha=0.5)
                        tmp_plt.contour(lon,lat,breeze_ds.dT.sel(datetime=time_str),levels=[breeze_ds.dT_cutoff],colors='limegreen')
                        tmp_plt.tick_params(labelsize=15,labelbottom=False)
                        tmp_plt.set_title('$\Delta_T$',size=18)


                        vel_plt = plt.subplot2grid((nrow,ncol),(1,0),colspan=2,aspect='equal',sharex=tmp_plt,sharey=tmp_plt)
                        vel_plt_pc = vel_plt.pcolormesh(lon,lat,vel10,cmap=plt.cm.YlOrRd)
                        vel_cbar = plt.colorbar(vel_plt_pc,ax=vel_plt)
                        vel_cbar.ax.tick_params(labelsize=13)
                        vel_plt.contour(lon,lat,land_mask,levels=[0.5],colors='k',alpha=0.5)
                        vel_plt.contour(lon,lat,breeze_ds.dU.sel(datetime=time_str),levels=[0.5],colors='b')
                        vel_plt.tick_params(labelsize=15)
                        vel_plt.set_title('Wind Speed',size=18)

                        ons_plt = plt.subplot2grid((nrow,ncol),(0,2),colspan=2,aspect='equal',sharex=tmp_plt,sharey=tmp_plt)
                        ons_plt_pc = ons_plt.pcolormesh(lon,lat,onshore_min,cmap=plt.cm.viridis,norm=Normalize(0,360))
                        ons_cbar = plt.colorbar(ons_plt_pc,ax=ons_plt)
                        ons_cbar.ax.tick_params(labelsize=13)
                        ons_plt.contour(lon,lat,land_mask,levels=[0.5],colors='k',alpha=0.5)
                        ons_plt.tick_params(labelsize=15,labelbottom=False)
                        ons_plt.set_title('Onshore Min',size=18)

                        dir_plt = plt.subplot2grid((nrow,ncol),(1,2),colspan=2,aspect='equal',sharex=tmp_plt,sharey=tmp_plt)
                        dir_plt_pc = dir_plt.pcolormesh(lon,lat,dir10,cmap=plt.cm.viridis,norm=Normalize(0,360))
                        dir_cbar = plt.colorbar(dir_plt_pc,ax=dir_plt)
                        dir_cbar.ax.tick_params(labelsize=13)
                        dir_plt.contour(lon,lat,land_mask,levels=[0.5],colors='k',alpha=0.5)
                        #dir_plt.contour(lon,lat,onshore_winds,levels=[1.0],colors='r')
                        qint = 7
                        dir_plt.quiver(lon[::qint,::qint],lat[::qint,::qint],u10[::qint,::qint],v10[::qint,::qint])
                        dir_plt.contour(lon,lat,breeze_ds.good_wdir.sel(datetime=time_str),levels=[1.0],colors='r')
                        dir_plt.tick_params(labelsize=15)
                        dir_plt.set_title('Wind Direction',size=18)


                        crt_plt = plt.subplot2grid((nrow,ncol),(2,0),rowspan=2,colspan=4,aspect='equal',sharex=tmp_plt,sharey=tmp_plt)
                        crt_plt_pc = crt_plt.pcolormesh(lon,lat,land_mask,cmap=plt.cm.binary,norm=Normalize(0,5.0))
                        #plt.colorbar(crt_plt_pc,ax=crt_plt)
                        crt_plt.contour(lon,lat,breeze_ds.good_wdir.sel(datetime=time_str),levels=[1.0],colors='r',alpha=0.5)
                        crt_plt.contour(lon,lat,breeze_ds.dT.sel(datetime=time_str),levels=[breeze_ds.dT_cutoff],colors='limegreen',alpha=0.5)
                        crt_plt.contour(lon,lat,breeze_ds.dU.sel(datetime=time_str),levels=[0.5],colors='b',alpha=0.5)
                        crt_plt.plot(0.0,0.0,c='r',label='Wind Direction')
                        crt_plt.plot(0.0,0.0,c='limegreen',label='∆T > 0')
                        crt_plt.plot(0.0,0.0,c='b',label='Wind Ramp')
                        crt_plt.set_xlim(np.nanmin(lon),np.nanmax(lon))
                        crt_plt.set_ylim(np.nanmin(lat),np.nanmax(lat))
                        good_areas = crt_plt.contourf(lon,lat,breeze_ds.breeze.sel(datetime=time_str),levels=[1.9,4.1],colors='k',
                                                      alpha=0.25,hatches=['+++', '/'],label='Test')
                        artists, labels = good_areas.legend_elements()
                        #crt_plt.legend(artists, ['Good Points'], handleheight=2)
                        crt_plt.legend(frameon=False,fontsize=18,loc=4)
                        crt_plt.tick_params(labelsize=15)
                        crt_plt.set_title(time_str,size=18)


                        for ss,stn in enumerate(near_shore_stations):
                            stn_color = int(breeze_ds.breeze.sel(datetime=time_str,
                                                               south_north=stn_dict[stn]['j'],
                                                               west_east=stn_dict[stn]['i']).data)
                            crt_plt.scatter(stn_dict[stn]['lon'],stn_dict[stn]['lat'],
                                        facecolor=breeze_cmap[stn_color],s=100,
                                        edgecolor='k',zorder=4,alpha=0.65)
                            crt_plt.text(stn_dict[stn]['lon']-0.1,stn_dict[stn]['lat'],stn,size=13,ha='right',va='center')


                        plt.show()


            breeze_ds_f['onshore_min'] = onshore_min
            breeze_ds_f['onshore_max'] = onshore_max
            breeze_ds_f['land_mask'] = land_mask

            breeze_ds_f.attrs = {'WRF_DIR':'/'.join(wrfout_files[0].split('/')[:-2]),
                                  'domain':wrfout_files[0].split('/')[-1].split('_')[1]}

            breeze_ds_f.to_netcdf(numerical_detection_path)
        else:
            print('{} has no wrfout files...'.format(case))
            bad_cases.append(case)
