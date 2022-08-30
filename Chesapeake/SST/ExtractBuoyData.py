#!/glade/u/home/hawbecke/local/envs/mmc/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from os import path
import xarray as xr
#from matplotlib.colors import Normalize, LinearSegmentedColormap
import sys
#sys.path.append('../../')
#from pydicts.obsdict import get_FINO_obs
from mmctools.helper_functions import theta_to_T, lowess_mean, calc_uv,w_s,T_d
from mmctools.wrf.utils import Tower, tsout_seriesReader, write_tslist_file
from mmctools.helper_functions import calc_wind
#from mmctools.plotting import TaylorDiagram
#import wrf
#from netCDF4 import Dataset
#import matplotlib.colors as colors
#import scipy.stats as stats
#import matplotlib.gridspec as gridspec
#from scipy.stats import pearsonr
#import matplotlib.patches as patches
#from string import ascii_lowercase


if len(sys.argv) > 1:
    if len(sys.argv) == 3:
        case_ind_s = int(sys.argv[1])
        case_ind_e = int(sys.argv[2])
    else:
        raise ValueError('ExtractBuoyData.py [case index start] [case index end]')
else:
    case_ind_s = 0
    case_ind_e = -1


sys.path.append('/glade/u/home/hawbecke/Code/Python/')
from publications.Chesapeake.CBB_case_dict import case_dict


cases = list(case_dict.keys())
cases = cases[3:-4]

cases = cases[case_ind_s:case_ind_e]

ncases   = len(cases)
case_dom = [3]*ncases

obs_f_dir = '/glade/work/hawbecke/ATEC/Chesapeake/Data/Obs/BUOY/'
obs_type = 'combined' # combined or CBIBS or NDBC
#obs_dir = '{}{}/'.format(obs_f_dir,obs_type)
buoy_f = '{}{}_buoy_data_res_2019.nc'.format(obs_f_dir,obs_type)
buoy_ds = xr.open_dataset(buoy_f)


buoy_list = list(buoy_ds.station.data)

bad_buoys = ['PMC','GRF','SRP','WDC','CAM','SLM','JTN']

for bad_buoy in bad_buoys:
    if bad_buoy in buoy_list:
        buoy_list.remove(bad_buoy)
buoy_loc_dict = {}
for stn in buoy_list:
    buoy_loc_dict[stn] = {'lat':float(buoy_ds.sel(station=stn).lat.data),
                          'lon':float(buoy_ds.sel(station=stn).lon.data)}

#buoy_loc_dict['JTN']['lat'] = 37.21137
#buoy_loc_dict['JTN']['lon'] = -76.78677


shallow_buoys = buoy_list.copy()

shallow_buoys.remove('DEB')
shallow_buoys.remove('VAB')

scratch_dir = '/glade/scratch/hawbecke/WRF/ATEC/Chesapeake/20190716to20190801/SST_SENSITIVITY/'
work_dir = '/glade/work/hawbecke/ATEC/Chesapeake/SENSITIVITY_STUDY/20190716to20190801/'
wrf_dir = scratch_dir
#work_dir   = '/glade/work/hawbecke/ATEC/Chesapeake/SENSITIVITY_STUDY/20190716to20190801/'

restarts = ['ERAI_20190715', 'ERAI_20190717', 'ERAI_20190719', 'ERAI_20190721', 
            'ERAI_20190723', 'ERAI_20190725', 'ERAI_20190727', 'ERAI_20190729']

wrf_start  = ['2019-07-15 18:00:00','2019-07-17 18:00:00','2019-07-19 18:00:00','2019-07-21 18:00:00',
              '2019-07-23 18:00:00','2019-07-25 18:00:00','2019-07-27 18:00:00','2019-07-29 18:00:00',]



case = 'DEFLT_NOFL_NOSK'
dom = 3
case_str = '{}'.format(case)
if dom == 3:
    dt = 9.0
elif dom == 4:
    dt = 3.0
else:
    print('Add logic for domain {}'.format(dom))
print('Starting {} d0{}: time step = {}'.format(case_str,dom,dt))


case_dir = '{}{}/'.format(wrf_dir,case)
twr_path = '{}{}_d0{}_towers.nc'.format(case_dir,case,dom)    
if path.exists(twr_path):
    print('loading in tower data!')
    towers = xr.open_dataset(twr_path)
else:
    towers = tsout_seriesReader(case_dir,restarts,wrf_start,'d0{}'.format(dom),structure='unordered',
                                        select_tower=np.append(near_shore_stations,ref_stn),time_step=dt,
                                        #select_tower=near_shore_stations[:2],time_step=dt,
                                        heights=[10.0],height_var='ph')
    towers['temp'] = theta_to_T(wrf_twrs[case_str].theta,wrf_twrs[case_str].pr/100.0)-273.15
    towers['wspd10'],towers['wdir10'] = calc_wind(wrf_twrs[case_str],u='u10',v='v10')
    towers['t2'] += -273.15 





wrf_buoys = {}
for cc,case in enumerate(cases[:]):
    new_f_name = '{0}{1}/{1}_extracted_{2}_buoy_data_d0{3}.nc'.format(wrf_dir,case,obs_type,case_dom[cc])
    if path.exists('{}'.format(new_f_name)):
        print('Data for {} already created!'.format(new_f_name.split('/')[-1]))
        wrf_buoys['{}_{}'.format(case,'d0{}'.format(case_dom[cc]))] = xr.open_dataset('{}'.format(new_f_name))
    else:
        print(case, new_f_name)
        get_wrf_locs = True
        rst_dict = {}
        for rr,rst in enumerate(restarts):
            f_list = sorted(glob.glob('{}{}/{}/wrfout_d0{}*'.format(scratch_dir,case,rst,case_dom[cc])))

            # Need to account for spinup...
            f_list = f_list[12:]

            if get_wrf_locs:
                wrf_loc_dict = {}
                wrf_i = xr.open_dataset('{}/wrfinput_d0{}'.format('/'.join(f_list[0].split('/')[:-1]),case_dom[cc]))
                is_water = [True]*len(buoy_loc_dict)
                hgt      = np.zeros(len(buoy_loc_dict))
                for ss,stn in enumerate(buoy_loc_dict):
                    #if stn in list(towers.station):
                    wrf_loc_dict[stn] = {}
                    loc_x = int(towers.sel(station=stn).i.data)
                    loc_y = int(towers.sel(station=stn).j.data)
                    # Some buoy locations are on land... grab the water cell next to it
                    stn_landmask = wrf_i.LANDMASK.sel(south_north=loc_y,west_east=loc_x).data
                    if stn_landmask == 1.0:
                        if stn == 'JTN':
                            loc_y -= 1
                            loc_x += 1
                        if stn == 'CAM':
                            loc_x -= 1
                            loc_y -= 1

                    wrf_loc_dict[stn]['x'] = loc_x
                    wrf_loc_dict[stn]['y'] = loc_y
                    wrf_stn = wrf_i.sel(south_north=loc_y,west_east=loc_x)
                    wrf_lat = wrf_stn.XLAT
                    wrf_lon = wrf_stn.XLONG
                    landmask = wrf_stn.LANDMASK
                    hgt[ss]  = wrf_stn.HGT
                    if landmask == 1:
                        is_water[ss] = False
                get_wrf_locs = False
                wrf_i.close()
            wrf_times = []
            wspd = np.zeros((len(f_list),len(buoy_loc_dict)))
            wdir = np.zeros((len(f_list),len(buoy_loc_dict)))
            t2   = np.zeros((len(f_list),len(buoy_loc_dict)))
            tsk  = np.zeros((len(f_list),len(buoy_loc_dict)))
            pres = np.zeros((len(f_list),len(buoy_loc_dict)))
            dwpt = np.zeros((len(f_list),len(buoy_loc_dict)))
            sst  = np.zeros((len(f_list),len(buoy_loc_dict)))
            

            for ff,fname in enumerate(f_list):
                print(fname)
                wrf_f = xr.open_dataset(fname,decode_times=False).sel(Time=0)
                wrf_time = pd.to_datetime(str(wrf_f.Times.data).replace("b'",'').replace("'",'').replace('_',' '))
                wrf_times.append(wrf_time)
                for ss,stn in enumerate(buoy_loc_dict):
                    u10 = wrf_f.U10.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])
                    v10 = wrf_f.V10.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])

                    wspd[ff,ss] = (u10**2 + v10**2)**0.5
                    wdir[ff,ss] = 180. + np.degrees(np.arctan2(u10, v10))
                    t2[ff,ss]   = wrf_f.T2.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])
                    tsk[ff,ss]  = wrf_f.TSK.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])
                    pres[ff,ss] = wrf_f.PSFC.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])/100.0
                    mixingratio = wrf_f.Q2.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])
                    w_sat       = w_s(t2[ff,ss],pres[ff,ss])
                    rh          = (mixingratio/w_sat)*100.0
                    dwpt[ff,ss] = T_d(t2[ff,ss],rh,celsius=True)
                    sst[ff,ss]  = wrf_f.SST.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])

                    #try:
                    #    sst[ff,ss]  = wrf_f.SST.sel(south_north=wrf_loc_dict[stn]['y'],west_east=wrf_loc_dict[stn]['x'])
                    #except KeyError:
                    #    if ff == 0 and ss == 0: print('no SST data found... setting to zero')
                    #    sst[ff,ss] = tsk[ff,ss]*0.0
                    

                wrf_f.close()
            rst_dict[rst] = xr.Dataset({'wspd': (['datetime','station'],wspd),
                                        'wdir': (['datetime','station'],wdir),
                                        't2'  : (['datetime','station'],t2),
                                        'tsk' : (['datetime','station'],tsk),
                                        'sst' : (['datetime','station'],sst),
                                        'pres': (['datetime','station'],pres),
                                        'dwpt': (['datetime','station'],dwpt),
                                    'is_water': (['station'],is_water),
                                         'hgt': (['station'],hgt)},
                                        coords={'datetime':wrf_times,
                                                'station': list(buoy_loc_dict.keys())})             
            if rr == 0:
                full_ds = rst_dict[rst]
            else:
                full_ds = full_ds.combine_first(rst_dict[rst])

        full_ds.to_netcdf(new_f_name)
        wrf_buoys['{}_{}'.format(case,'d0{}'.format(case_dom[cc]))]  = full_ds


        #fig = plt.figure(figsize=(18,6))
        #full_ds.sel(station='ANN').dwpt.plot(lw=2.0,c='k')
        #for rr,rst in enumerate(rst_dict):
        #    rst_dict[rst].sel(station='ANN').dwpt.plot(lw=1.0)
        #plt.show()
        
