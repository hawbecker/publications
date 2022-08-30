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

#cases = cases[case_ind_s:case_ind_e]
#for case in cases: print(case)
#wefwef
ncases   = len(cases)
case_dom = [3]*ncases

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


wrf_buoys = {}

init_vars = True
for cc,case in enumerate(cases[:]):
    new_f_name = '{0}{1}/{1}_extracted_flux_data_d0{2}.nc'.format(wrf_dir,case,case_dom[cc])
    if path.exists('{}'.format(new_f_name)):
        print('Data for {} already created!'.format(new_f_name.split('/')[-1]))
    else:
        print(case, new_f_name)
        get_wrf_locs = True
        rst_dict = {}
        for rr,rst in enumerate(restarts):
            f_list = sorted(glob.glob('{}{}/{}/wrfout_d0{}*'.format(scratch_dir,case,rst,case_dom[cc])))

            # Need to account for spinup...
            f_list = f_list[1:]
            
            wrf_times = []
            hfx_l = np.zeros(len(f_list))
            hfx_w = np.zeros(len(f_list))
            qfx_l = np.zeros(len(f_list))
            qfx_w = np.zeros(len(f_list))
            lh_l  = np.zeros(len(f_list))
            lh_w  = np.zeros(len(f_list))
            
            for ff,fname in enumerate(f_list):
                print(fname)
                wrf_f = xr.open_dataset(fname,decode_times=False).sel(Time=0)
                wrf_time = pd.to_datetime(str(wrf_f.Times.data).replace("b'",'').replace("'",'').replace('_',' '))
                wrf_times.append(wrf_time)
                
                wrf_f = wrf_f[['LANDMASK','HFX','QFX','LH']]
                wrf_f = wrf_f.sel(south_north=slice(40,135),
                                  west_east=slice(75,128))
                
                if init_vars:
                    wrf_landmask = wrf_f.LANDMASK
                    init_vars = False
                    
                wrf_hfx = wrf_f.HFX
                wrf_qfx = wrf_f.QFX
                wrf_lh  = wrf_f.LH
                wrf_f.close()
                
                wrf_hfx_land  = wrf_hfx.where(wrf_landmask==1.0)
                wrf_hfx_water = wrf_hfx.where(wrf_landmask==0.0)
                wrf_qfx_land  = wrf_qfx.where(wrf_landmask==1.0)
                wrf_qfx_water = wrf_qfx.where(wrf_landmask==0.0)
                wrf_lh_land  = wrf_lh.where(wrf_landmask==1.0)
                wrf_lh_water = wrf_lh.where(wrf_landmask==0.0)
                
                hfx_l[ff] = wrf_hfx_land.mean().data
                hfx_w[ff] = wrf_hfx_water.mean().data
                qfx_l[ff] = wrf_qfx_land.mean().data
                qfx_w[ff] = wrf_qfx_water.mean().data
                lh_l[ff] = wrf_lh_land.mean().data
                lh_w[ff] = wrf_lh_water.mean().data
                
            rst_dict[rst] = xr.Dataset({'hfx_land': (['datetime'],hfx_l),
                                        'qfx_land': (['datetime'],qfx_l),
                                        'lh_land' : (['datetime'],lh_l),
                                        'hfx_water': (['datetime'],hfx_w),
                                        'qfx_water': (['datetime'],qfx_w),
                                        'lh_water' : (['datetime'],lh_w)},
                                        coords={'datetime':wrf_times}
                                        )

            if rr == 0:
                full_ds = rst_dict[rst]
            else:
                full_ds = full_ds.combine_first(rst_dict[rst])

        full_ds.to_netcdf(new_f_name)
