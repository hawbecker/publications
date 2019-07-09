#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sowfadict as sowfa
import tools
import pandas as pd
from matplotlib.colors import Normalize
import gzip
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from sys import stdout
from time import sleep
import pickle


simstr  = ['PGx2.2En4_zi550','PGx2.2En4_zi1100','PGx4.4En4_zi550']#,'PGx4.4En4_zi1100']
#cases   = ['M2_SR_PertBC','M2_D2_PertBC']
nsims   = np.shape(simstr)[0]

fdir     = '/projects/mmc/NWTCRegion/'
savedir  = '/home/phawbeck/SOWFA/Terrain/data/'
saveData = True
saveFigs = True

timedir = 30600
zi_cut = 302.5
resamp_norm_z = np.arange(0.0,1.21,0.01)
nx = 896; ny = 208; nz = 75

for ss in np.arange(0,nsims):
    for cs in np.arange(0,2):
        if cs == 0: 
            terr_str = 'SR'
        if cs == 1:
            terr_str = 'D2'
        cases = ['M2_{}_TurbBC'.format(terr_str),'M2_{}_NoTurbBC'.format(terr_str),'M2_{}_PertBC'.format(terr_str)]
        ncases  = np.shape(cases)[0]
        print('Starting {}'.format(simstr[ss]))
        print('Getting x,y,z')

        xxf = gzip.open('{}{}/NWTC.run.{}/constant/ccx.gz'.format(fdir,simstr[ss],cases[0]),'r')
        xx = xxf.readlines()[22:-15]
        xxf.close()

        yyf = gzip.open('{}{}/NWTC.run.{}/constant/ccy.gz'.format(fdir,simstr[ss],cases[0]),'r')
        yy = yyf.readlines()[22:-15]
        yyf.close()

        zzf = gzip.open('{}{}/NWTC.run.{}/constant/ccz.gz'.format(fdir,simstr[ss],cases[0]),'r')
        zz = zzf.readlines()[22:-15]
        zzf.close()

        print('Getting Temperature for zi calculation')

        Thf = gzip.open('{}{}/NWTC.run.{}/{}/TAvg.gz'.format(fdir,simstr[ss],cases[0],timedir),'r')
        Th1 = Thf.readlines()[22:-15]
        Thf.close()
        Thf = gzip.open('{}{}/NWTC.run.{}/{}/TAvg.gz'.format(fdir,simstr[ss],cases[1],timedir),'r')
        Th2 = Thf.readlines()[22:-15]
        Thf.close()
        Thf = gzip.open('{}{}/NWTC.run.{}/{}/TAvg.gz'.format(fdir,simstr[ss],cases[2],timedir),'r')
        Th3 = Thf.readlines()[22:-15]
        Thf.close()

        

        x = np.zeros((nx,ny,nz))
        y = np.zeros((nx,ny,nz))
        z = np.zeros((nx,ny,nz))
        T = np.zeros((3,nx,ny,nz))
        tke = np.zeros((3,nx,ny,nz))


        print('shaping arrays')
        ct = 0
        for kk in range(0,nz):
            for jj in range(0,ny):
                for ii in range(0,nx):
                    x[ii,jj,kk] = np.float(xx[ct])
                    y[ii,jj,kk] = np.float(yy[ct])
                    z[ii,jj,kk] = np.float(zz[ct])

                    Tline   = Th1[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[0,ii,jj,kk]   = np.float(Tline[0])

                    Tline   = Th2[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[1,ii,jj,kk]   = np.float(Tline[0])

                    Tline   = Th3[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[2,ii,jj,kk]   = np.float(Tline[0])
                    ct += 1


        print('Calculating zi')
        zi     = np.zeros((ncases,nx,ny))
        zi_agl = np.zeros((ncases,nx,ny))
        norm_z = np.zeros((ncases,nx,ny,nz))
        for jj in range(0,ny):
            for ii in range(0,nx):
                zloc = z[ii,jj,:]*1.0
                for case in np.arange(0,ncases):
                    zi_agl[case,ii,jj] = zloc[T[case,ii,jj,:]>=zi_cut][0]
                zloc -= zloc[0]
                for case in np.arange(0,ncases):
                    zi[case,ii,jj] = zloc[T[case,ii,jj,:]>=zi_cut][0]
                    norm_z[case,ii,jj,:] = zloc/zi[case,ii,jj]
        pickle.dump( zi, open( "{}{}_{}_zi.p".format(savedir,simstr[ss].replace('.','_'),terr_str), "wb" ))

