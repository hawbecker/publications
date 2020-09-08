import netCDF4 as nc
from netCDF4 import Dataset as ncdf
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from datetime import datetime,timedelta
import glob
import shutil
from os import path
import pandas as pd


# -------------------------------------------------------------------------------- #
#fdir = '/glade/scratch/hawbecke/WRF/ATEC/Chesapeake/SENSITIVITY_STUDY/'
fdir = '/glade/work/hawbecke/MMC/FINO/'
#fdir = '/glade/work/hawbecke/ATEC/Chesapeake/SENSITIVITY_STUDY/'
reanalysis = 'ERA5'
#met_dir = '{}20190716to20190801/met_em/{}/orig/'.format(fdir,reanalysis)
met_dir = '{}met_em/{}/orig/'.format(fdir,reanalysis)

smoothed = False
#new_data = 'GHRSST' 
#new_data = 'MODIS'
new_data = 'OVERWRITE'
#if reanalysis == 'ERA5': new_data = 'OVERWRITE'

#new_dir = '{}20190716to20190801/met_em/{}/{}/'.format(fdir,reanalysis,new_data)
new_dir = '{}met_em/{}/{}/'.format(fdir,reanalysis,new_data)
sst_dir = '{}{}/'.format(fdir,new_data)
# -------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------- #

ghr_valid = 9  # GHRSST is valid at 09 Z
mod_valid = 12 # Assuming it is valid at 12 Z for interpolation purposes...


met_list = sorted(glob.glob('{0}met_em.d0*'.format(met_dir)))

print(met_list)
wefwef
def fill_new_sst(orig_sst,orig_lat,orig_lon,orig_mask,
                 new_sst,new_lat,new_lon,window=None):
    '''
    example call: fill_new_sst(met_sst,met_lat,met_lon,met_mask,
                               newdata_sst,newdata_lat,newdata_lon,
                               window=2)
    
    '''
    if np.shape(np.shape(orig_lat))[0] == 1:
        ny = np.shape(orig_lat)[0]
    elif np.shape(np.shape(orig_lat))[0] == 2:
        ny = np.shape(orig_lat)[0]
    else:
        print('Original latitude must be 1 or 2 dimensions.')
        return
    
    if np.shape(np.shape(orig_lon))[0] == 1:
        nx = np.shape(orig_lon)[0]
    elif np.shape(np.shape(orig_lon))[0] == 2:
        nx = np.shape(orig_lon)[1]
    else:
        print('Original latitude must be 1 or 2 dimensions.')
        return

    sst = orig_sst.copy()
    for jj in np.arange(0,ny):
        for ii in np.arange(0,nx):
            if orig_mask[jj,ii] == 0.0:
                dist_lat = (new_lat - orig_lat[jj,ii])**2
                lat_ind = np.where(dist_lat==np.min(dist_lat))[0][0]
                dist_lon = (new_lon - orig_lon[jj,ii])**2
                lon_ind = np.where(dist_lon==np.min(dist_lon))[0][0]
                
                if window == None:
                    new_val = float(new_sst[lat_ind,lon_ind])
                    if np.isnan(new_val):
                        new_val = orig_sst[jj,ii]
                    sst[jj,ii] = new_val
                else:
                    ywindow_s = max([0,lat_ind-window])
                    ywindow_e = min([np.shape(new_sst)[0],lat_ind+window+1])

                    xwindow_s = max([0,lon_ind-window])
                    xwindow_e = min([np.shape(new_sst)[1],lon_ind+window+1])
                    new_val = np.asarray(new_sst[ywindow_s:ywindow_e, lon_ind-window:lon_ind+window+1],dtype=float)
                    new_val[new_val<230] = np.nan
                    new_val[new_val>1e20] = np.nan
                    
                    new_val = np.nanmean(new_val)
                    if np.size(new_val) != 1 or np.isnan(new_val):
                        sst[jj,ii] = orig_sst[jj,ii]
                    else:
                        sst[jj,ii] = new_val
    sst[sst == 0] = np.nan
    return(sst.data)



for met_file in met_list:
    if smoothed:
        smooth_str = 'smooth'
    else:
        smooth_str = 'raw'

    if reanalysis == 'ERA5':
        new_file = '{}{}'.format(new_dir,met_file.split('/')[-1])
    else:
        new_file = '{}{}/{}'.format(new_dir,smooth_str,met_file.split('/')[-1])
    print(new_file)
    if path.exists(new_file):
        print('skipping...')
    else:

        met = ncdf(met_file,'r')

        m_lat = met.variables['XLAT_M'][0,:,:]
        m_lon = met.variables['XLONG_M'][0,:,:]
        if reanalysis == 'GFS':
            m_sst = met.variables['SKINTEMP'][0,:,:]
        else:
            m_sst = met.variables['SST'][0,:,:]
        m_mask = met.variables['LANDMASK'][0,:,:]
        met_timeb = met.variables['Times'][:].data[0]
        if reanalysis == 'ERA5': 
            new_sst = met.variables['SKINTEMP'][0,:,:]
            new_sst[np.where(m_mask == 1.0)] = new_sst[np.where(m_mask == 1.0)]*0.0
        met_dx = met.DX
        met_dy = met.DY
        met.close()

#        if reanalysis != 'ERA5':
        if new_data != 'OVERWRITE':
            min_m_lat = np.min(m_lat)
            max_m_lat = np.max(m_lat)
            min_m_lon = np.min(m_lon)
            max_m_lon = np.max(m_lon)

            m_sst[m_sst==0.0] = np.nan


            met_time = ''
            for tt in met_timeb:
                met_time = met_time+tt.decode("utf-8")
            met_time = datetime.strptime(met_time,'%Y-%m-%d_%H:%M:%S')


            if new_data == 'GHRSST':
                valid = ghr_valid
            elif new_data == 'MODIS':
                valid = mod_valid
            else:
                print('Must be GHRSST or MODIS')

            if met_time.hour < valid:
                time_s = met_time - pd.to_timedelta(1,'D')
                time_e = met_time
            elif met_time.hour == valid:
                time_s = met_time
                time_e = met_time
            else:
                time_s = met_time
                time_e = met_time + pd.to_timedelta(1,'D')

            print('Getting {}'.format(new_data))
            if new_data == 'GHRSST':
                print('{0}{1:04d}{2:02d}{3:02d}*'.format(sst_dir,time_s.year,time_s.month,time_s.day))
                new_first = glob.glob('{0}{1:04d}{2:02d}{3:02d}*'.format(sst_dir,time_s.year,time_s.month,time_s.day))[0]
                new_last  = glob.glob('{0}{1:04d}{2:02d}{3:02d}*'.format(sst_dir,time_e.year,time_e.month,time_e.day))[0]
                valid = ghr_valid
            elif new_data == 'MODIS':
                new_first = glob.glob('{0}composite.{1:04d}{2:02d}{3:02d}*'.format(
                                      sst_dir,time_s.year,time_s.month,time_s.day))[0]
                new_last  = glob.glob('{0}composite.{1:04d}{2:02d}{3:02d}*'.format(
                                      sst_dir,time_e.year,time_e.month,time_e.day))[0]
                valid = mod_valid
            else:
                print('Must be GHRSST or MODIS')

            if new_data == 'GHRSST':
                new_1 = ncdf('{}'.format(new_first),'r')
                new_2 = ncdf('{}'.format(new_last),'r')
                n_lat = new_1.variables['lat']
                n_lon = new_1.variables['lon']

                n_lat_s = np.where(n_lat>=min_m_lat)[0][0]
                n_lat_e = np.where(n_lat<=max_m_lat)[0][-1]
                n_lon_s = np.where(n_lon>=min_m_lon)[0][0]
                n_lon_e = np.where(n_lon<=max_m_lon)[0][-1]


                n_sst_1 = new_1.variables['analysed_sst'][0,n_lat_s:n_lat_e,n_lon_s:n_lon_e]
                n_sst_2 = new_2.variables['analysed_sst'][0,n_lat_s:n_lat_e,n_lon_s:n_lon_e]
                n_mask = new_1.variables['mask'][0,n_lat_s:n_lat_e,n_lon_s:n_lon_e]

                nlon,nlat = np.meshgrid(n_lon[n_lon_s:n_lon_e],n_lat[n_lat_s:n_lat_e])

                new_dx = 1000.0
                new_dy = 1000.0

            elif new_data == 'MODIS':
                new_1 = ncdf('{}'.format(new_first),'r')
                new_2 = ncdf('{}'.format(new_last),'r')

                n_lat = new_1.variables['latitude']
                n_lon = new_1.variables['longitude']

                n_lat_e = np.where(n_lat>=min_m_lat)[0][-1]
                n_lat_s = np.where(n_lat<=max_m_lat)[0][0]
                n_lon_s = np.where(n_lon>=min_m_lon)[0][0]
                n_lon_e = np.where(n_lon<=max_m_lon)[0][-1]


                n_sst_1 = new_1.variables['sst_data'][0,n_lat_s:n_lat_e,n_lon_s:n_lon_e]+273.15
                n_sst_2 = new_2.variables['sst_data'][0,n_lat_s:n_lat_e,n_lon_s:n_lon_e]+273.15

                nlon,nlat = np.meshgrid(n_lon[n_lon_s:n_lon_e],n_lat[n_lat_s:n_lat_e])

                new_dx = 4625.0
                new_dy = 4625.0
            else:
                print('Must be GHRSST or MODIS')
           
            start_time = datetime.strptime('{0:04d}-{1:02d}-{2:02d} {3:02d}'.format(
                                           time_s.year,time_s.month,time_s.day,valid),'%Y-%m-%d %H')
            end_time   = datetime.strptime('{0:04d}-{1:02d}-{2:02d} {3:02d}'.format(
                                           time_e.year,time_e.month,time_e.day,valid),'%Y-%m-%d %H')

            d1 = (met_time - start_time).seconds
            d2 = (end_time - met_time).seconds
            if d1 == 0 & d2 == 0:
                n_sst = n_sst_1
            else:
                w1 = d2/(d1+d2)
                w2 = d1/(d1+d2)

                n_sst = n_sst_1*w1 + n_sst_2*w2

            n_lat_trim = nlat[:,0]
            n_lon_trim = nlon[0,:]

            if not smoothed:
                new_sst = fill_new_sst(m_sst,m_lat,m_lon,m_mask,n_sst,n_lat_trim,n_lon_trim)
            else:
                new_sst = fill_new_sst(m_sst,m_lat,m_lon,m_mask,n_sst,n_lat_trim,n_lon_trim,
                                    window=int(min([met_dx/new_dx,met_dy/new_dy])/2.0))

        new_sst[np.isnan(new_sst)] = 0.0

        shutil.copy2(met_file,new_file)
        new = ncdf(new_file,'a')
        if reanalysis == 'GFS':
            sst = new.variables['SKINTEMP']
        else:
            sst = new.variables['SST']
        sst[0,:,:] = new_sst
        sst.source = '{}{}'.format(new_data,smooth_str)
        new.close()

