#!/glade/u/home/hawbecke/local/envs/mmc/bin/python
from mmctools.wrf.preprocessing import OverwriteSST
import os
import glob
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

smoothed = False

reanalysis = 'ERAI'
# Base directory where everything is stored:
fdir = '/glade/work/hawbecke/ATEC/Chesapeake/SENSITIVITY_STUDY/20190716to20190801/'

# Location of met_em files:
met_dir = '{}met_em/{}/orig/'.format(fdir,reanalysis)


for new_data in ['NCEI', 'G1SST', 'MUR', 'MODIS', 'OSPO', 'OSTIA', 'NAVO', 'CMC', 'FILL'][4:8]:
    print('starting {}'.format(new_data))
    # Location of where you want to save new met_em files:
    out_dir = '{}met_em/{}/{}/'.format(fdir,reanalysis,new_data)
    num_new_met = len(glob.glob('{}raw/met_em*'.format(out_dir)))
    num_org_met = len(glob.glob('{}met_em*'.format(met_dir)))

    if num_new_met < num_org_met:
        # Location of SST data to be included in new met_em files:
        sst_dir = '{}SST/{}/'.format(fdir,new_data)

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            
        test = OverwriteSST(met_type=reanalysis,
                            overwrite_type=new_data,
                            met_directory=met_dir,
                            sst_directory=sst_dir,
                            out_directory=out_dir,
                            smooth_opt=False,
                            fill_missing=True,
                            skip_finished=True)
        
