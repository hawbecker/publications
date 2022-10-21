#!/Users/hawbecke/.conda/envs/pyhawbeck/bin/python
import pandas as pd
from mmctools.wrf.preprocessing import ERA5
import sys
import xarray as xr
import pygrib
import glob
import requests
import os
import matplotlib.pyplot as plt

hrrr_type = 'analysis'
hrrr_dir = '/Users/hawbecke/Research/Chesapeake/Data/HRRR/YEAR/{}/'.format(hrrr_type)
xs = 1400
ys = 500
span = 201
h_s, h_e = 13,22


if hrrr_type == 'analysis':
    # For paper:
    #date_range = pd.date_range('2019-07-16 {0:02d}:00:00'.format(h_s),
    #                           '2019-07-31 {0:02d}:00:00'.format(h_e),freq='1h')
    date_range = pd.date_range('2020-01-01 {0:02d}:00:00'.format(h_s),
                               '2020-12-31 {0:02d}:00:00'.format(h_e),freq='1h')
elif hrrr_type == 'forecast':
    forecast_s = 10
    forecast_hrs = range(0,19)
    date_range = pd.date_range('2019-07-16 {0:02d}:00:00'.format(forecast_s),
                               '2019-07-31 {0:02d}:00:00'.format(forecast_s),freq='1d')


def download_HRRR(datetime,save_dir,forecast_hr=None):
    if forecast_hr is None:
        url = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/{0:04d}{1:02d}{2:02d}/hrrr.t{3:02d}z.wrfsfcf00.grib2'.format(
                                                                datetime.year,datetime.month,datetime.day,datetime.hour)
        new_f_name = 'hrrr.{0:04d}{1:02d}{2:02d}.t{3:02d}z.wrfsfcf00.grib2'.format(
                        datetime.year,datetime.month,datetime.day,datetime.hour)
    else:
        url = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/{0:04d}{1:02d}{2:02d}/hrrr.t{3:02d}z.wrfsfcf{4:02d}.grib2'.format(
                                                                datetime.year,datetime.month,datetime.day,datetime.hour,forecast_hr)
        new_f_name = 'hrrr.{0:04d}{1:02d}{2:02d}.t{3:02d}z.wrfsfcf{4:02d}.grib2'.format(
                        datetime.year,datetime.month,datetime.day,datetime.hour,forecast_hr)
        

    print(url)
    if os.path.exists('{}{}'.format(save_dir,new_f_name)) or (os.path.exists('{}{}'.format(save_dir,new_f_name.replace('grib2','nc')))):
        print('{}{} already downloaded.'.format(save_dir,new_f_name))
    else:
        r = requests.get(url, allow_redirects=True)
        #new_f_name = 'hrrr.{0:04d}{1:02d}{2:02d}.t{3:02d}z.wrfsfcf00.grib2'.format(
        #                datetime.year,datetime.month,datetime.day,datetime.hour)
        open('{}{}'.format(save_dir,new_f_name), 'wb').write(r.content)
    return(new_f_name)
    
def convert_GRIB2_NetCDF(grib_dir,grib_name,
                         ncl_loc='/Users/hawbecke/.conda/envs/ncl_stable/bin/ncl',
                         run_script=True):
    
    aux_cmd = 'fname="{}"'.format(grib_name.replace('.grib2',''))
    scrpt_n = 'grib2netcdf.ncl'
    arg_str = 'cd {} && '.format(grib_dir) + ncl_loc + " '"+aux_cmd+"' " + scrpt_n #+ ' && cd -'
    if run_script:
        os.system(arg_str)
    else:
        print(arg_str)
    
    new_f_name = grib_name.replace('grib2','nc')
    os.remove('{}{}'.format(grib_dir,grib_name))
    return(new_f_name)
    
def crop_NetCDF(fname,xs,ys,span):
    ncd = xr.open_dataset(fname)
    print(ncd)
    ncd_small = ncd.sel(xgrid_0=slice(xs,xs+span),ygrid_0=slice(ys,ys+span))

    # Moved the variable extraction into the NCL code...
    '''
    good_vars = ['LAND_P0_L1_GLC0',
                 'HGT_P0_L1_GLC0',
                 'gridlat_0',
                 'gridlon_0',
                 'TMP_P0_L103_GLC0',
                 'POT_P0_L103_GLC0',
                 'UGRD_P0_L103_GLC0',
                 'VGRD_P0_L103_GLC0',
                 'PRES_P0_L1_GLC0',
                 'HGT_P0_L100_GLC0',
                 'TMP_P0_L100_GLC0',
                 'POT_P0_L103_GLC0',
                 'UGRD_P0_L100_GLC0',
                 'VGRD_P0_L100_GLC0',
                 'lv_HTGL2',
                 'lv_ISBL0',
                 'lv_ISBL1',
                 'lv_ISBL5',
                ]
    for varn in ncd_small.variables:
        if varn not in good_vars:
            ncd_small = ncd_small.drop(varn)
    '''
    
    os.system('rm {}'.format(fname))
    ncd_small.to_netcdf('{}'.format(fname))
    
    return(ncd_small)

def write_grib2ncdf_code(grib_dir):
    if 'forecast' in grib_dir:
        precip_var = 'APCP_P8_L1_GLC0_acc1h'
    else:
        precip_var = 'APCP_P8_L1_GLC0_acc","REFC_P0_L10_GLC0'
        
    ncl_file = open('{}grib2netcdf.ncl'.format(grib_dir),'w')
    ncl_file.write('begin\n')
    ncl_file.write(';***********************************************\n')
    ncl_file.write('; get variable names from grib file\n')
    ncl_file.write(';***********************************************\n')
    ncl_file.write(';    fname    = "hrrr.t15z.wrfsfcf00"\n')
    ncl_file.write('    grib_in  = addfile("./"+fname+".grib2","r")\n')
    ncl_file.write(';    names    = getfilevarnames(grib_in); extract all variable names\n')
    ncl_file.write('    names = (/"LAND_P0_L1_GLC0","HGT_P0_L1_GLC0","gridlat_0","gridlon_0","TMP_P0_L103_GLC0","POT_P0_L103_GLC0","UGRD_P0_L103_GLC0","VGRD_P0_L103_GLC0","PRES_P0_L1_GLC0","HGT_P0_L100_GLC0","TMP_P0_L100_GLC0","POT_P0_L103_GLC0","UGRD_P0_L100_GLC0","VGRD_P0_L100_GLC0","lv_HTGL2","lv_ISBL0","lv_ISBL1","lv_ISBL5","LCDC_P0_L214_GLC0","TCDC_P0_L10_GLC0","MCDC_P0_L224_GLC0","HCDC_P0_L234_GLC0","{}"/)\n'.format(precip_var))
    ncl_file.write('    ;***********************************************\n')
    ncl_file.write('    ; create output netcdf file\n')
    ncl_file.write('    ;***********************************************\n')
    ncl_file.write('    system("rm " + fname + ".nc") ; remove any pre-existing file\n')
    ncl_file.write('    ncdf_out = addfile(fname + ".nc" ,"c")       ; create output netCDF file\n')
    ncl_file.write('    ;***********************************************\n')
    ncl_file.write('    ; loop through variables and output each to netcdf\n')
    ncl_file.write('    ;***********************************************\n')
    ncl_file.write('    do i = 0, dimsizes(names)-1\n')
    ncl_file.write('        ncdf_out->$names(i)$ = grib_in->$names(i)$\n')
    ncl_file.write('    end do\n')
    ncl_file.write('end\n')
    ncl_file.close()


# - - - - - - - - - - - - - - - #


write_grib2ncdf_code(hrrr_dir)
for dd,date in enumerate(date_range):
    if hrrr_type == 'analysis':
        if (date.hour >= h_s) & (date.hour <= h_e):
            print(date)
            grib_n = download_HRRR(date,hrrr_dir)
            if not os.path.exists('{}{}'.format(hrrr_dir,grib_n.replace('grib2','nc'))):
                hrrr_fname = convert_GRIB2_NetCDF(hrrr_dir,grib_n)
                new_nc = crop_NetCDF('{}{}'.format(hrrr_dir,hrrr_fname),xs,ys,span)
    elif hrrr_type == 'forecast':
        for forecast_hr in forecast_hrs:
            hr = (date + pd.to_timedelta(forecast_hr,'h')).hour
            if (hr >= h_s) & (hr <= h_e):
                print('{} - forecast hour: {}\t hour: {}'.format(date,forecast_hr,hr))

                grib_n = download_HRRR(date,hrrr_dir,forecast_hr=forecast_hr)
                hrrr_fname = convert_GRIB2_NetCDF(hrrr_dir,grib_n)
                new_nc = crop_NetCDF('{}{}'.format(hrrr_dir,hrrr_fname),xs,ys,span)
                
