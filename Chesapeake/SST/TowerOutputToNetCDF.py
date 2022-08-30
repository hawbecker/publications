#!/glade/u/home/hawbecke/local/envs/mmc/bin/python
from os import path
import xarray as xr
from mmctools.wrf.utils import tsout_seriesReader, write_tslist_file
from mmctools.helper_functions import theta_to_T, lowess_mean, calc_uv,w_s,T_d
import numpy as np

import sys
sys.path.append('/glade/u/home/hawbecke/Code/Python/publications/')
from Chesapeake.CBB_case_dict import case_dict
cases = list(case_dict.keys())
cases = cases[3:7]
print(cases)

ncases   = len(cases)
case_dom = [3]*ncases

wrf_dir    = '/glade/scratch/hawbecke/WRF/ATEC/Chesapeake/20190716to20190801/SST_SENSITIVITY/'
#work_dir   = '/glade/work/hawbecke/ATEC/Chesapeake/SENSITIVITY_STUDY/20190716to20190801/'

restarts = ['ERAI_20190715', 'ERAI_20190717', 'ERAI_20190719', 'ERAI_20190721', 
            'ERAI_20190723', 'ERAI_20190725', 'ERAI_20190727', 'ERAI_20190729']

wrf_start  = ['2019-07-15 18:00:00','2019-07-17 18:00:00','2019-07-19 18:00:00','2019-07-21 18:00:00',
              '2019-07-23 18:00:00','2019-07-25 18:00:00','2019-07-27 18:00:00','2019-07-29 18:00:00',]

ref_stn = 'IAD'
wrf_twrs = {}
for cc,case in enumerate(cases):
    dom = case_dom[cc]
    case_str = '{}_{}'.format(case,'d0{}'.format(dom))
    if dom == 3:
        dt = 9.0
    elif dom == 4:
        dt = 3.0
    else:
        print('Add logic for domain {}'.format(dom))
    print('Starting {} d0{}: time step = {}'.format(case_str,dom,dt))
    

    case_dir = '{}{}/'.format(wrf_dir,case)
    twr_path = '{}{}_d0{}_towers2.nc'.format(case_dir,case,dom)    
    if path.exists(twr_path):
        print('loading in full dataset!')
        towers = xr.open_dataset(twr_path)
    else:
        case_dir = '{}{}/'.format(wrf_dir,case)

        towers = tsout_seriesReader(case_dir,restarts,wrf_start,'d0{}'.format(dom),structure='unordered',
                                    time_step=dt,
                                    #select_tower=np.append(near_shore_stations,ref_stn),time_step=dt,
                                    heights=[10.0],height_var='ph')
        towers['temp'] = theta_to_T(towers.theta,towers.pr/100.0)-273.15
        towers['wspd10'],towers['wdir10'] = calc_wind(towers,u='u10',v='v10')
        towers['t2'] += -273.15 
        towers.to_netcdf(twr_path)
