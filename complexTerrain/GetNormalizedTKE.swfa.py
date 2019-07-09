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


simstr  = ['PGx4.4En4_zi550']#'PGx2.2En4_zi550','PGx2.2En4_zi1100']#,,'PGx4.4En4_zi550']#,'PGx4.4En4_zi1100']
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

        
        print('Getting TKE')

        kf = gzip.open('{}{}/NWTC.run.{}/{}/kResolved.gz'.format(fdir,simstr[ss],cases[0],timedir),'r')
        tkef1 = kf.readlines()[22:-15]
        kf.close()
        kf = gzip.open('{}{}/NWTC.run.{}/{}/kResolved.gz'.format(fdir,simstr[ss],cases[1],timedir),'r')
        tkef2 = kf.readlines()[22:-15]
        kf.close()
        kf = gzip.open('{}{}/NWTC.run.{}/{}/kResolved.gz'.format(fdir,simstr[ss],cases[2],timedir),'r')
        tkef3 = kf.readlines()[22:-15]
        kf.close()


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
                    tkeline = tkef1[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[0,ii,jj,kk]   = np.float(Tline[0])
                    tke[0,ii,jj,kk] = np.float(tkeline[0])

                    Tline   = Th2[ct].replace(b'(',b'').replace(b')',b'').split()
                    tkeline = tkef2[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[1,ii,jj,kk]   = np.float(Tline[0])
                    tke[1,ii,jj,kk] = np.float(tkeline[0])

                    Tline   = Th3[ct].replace(b'(',b'').replace(b')',b'').split()
                    tkeline = tkef3[ct].replace(b'(',b'').replace(b')',b'').split()
                    T[2,ii,jj,kk]   = np.float(Tline[0])
                    tke[2,ii,jj,kk] = np.float(tkeline[0])
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

        if saveFigs:
            print('Saving Figures')
            yloc = int(ny/2)
            fig  = plt.figure(figsize=(16,4))

            for case in np.arange(0,ncases):
                pltnum = 231 + case
                plt1 = plt.subplot(pltnum)#,aspect='equal')
                cont = plt1.pcolormesh(x[:,yloc,:],z[:,yloc,:],tke[int(case),:,yloc,:],norm=colors.LogNorm(vmin=1e-1,vmax=1e1))#,cmap=plt.cm.RdBu)
                pltnum = 234 + case
                plt2 = plt.subplot(pltnum,sharex=plt1)#,aspect='equal')
                plt2.pcolormesh(x[:,yloc,:],z[:,yloc,:],T[case,:,yloc,:])#,norm=Normalize(0,3))#,cmap=plt.cm.RdBu)
                plt2.scatter(x[:,yloc,0],zi_agl[case,:,yloc],s=4,c='k',alpha=0.4)

            cax1 = fig.add_axes([0.92, 0.15, 0.01, 0.7])
            plt.colorbar(cont,cax1)
            figname1 = '{}_{}_PBLHeight'.format(simstr[ss],terr_str)
            print('{}{}.png'.format(savedir,figname1))
            plt.savefig('{}{}.png'.format(savedir,figname1))
            plt.clf()

#            plt.figure(figsize=(9,6))
#            plt.pcolormesh(x[:,yloc,:],norm_z[:,yloc,:],T[:,yloc,:])#,norm=Normalize(0,1.1))
#            plt.colorbar()
#            plt.ylim(-0.01,1.2)
#            plt.axhline(1.0,c='r',ls=':')
#            figname2 = '{}_{}_NormHeight'.format(simstr[ss],cases[cs])
#            print('{}{}.png'.format(savedir,figname2))
#            plt.savefig('{}{}.png'.format(savedir,figname2))
#            plt.clf()

        print('Interpolating to common normalized levels')
        int_tke_f     = np.zeros((ncases+2,nx,ny,resamp_norm_z.size))
        int_T_f       = np.zeros((ncases+2,nx,ny,resamp_norm_z.size))
        avg_int_tke   = np.zeros((ncases+2,nx,resamp_norm_z.size))
        avg_int_T     = np.zeros((ncases+2,nx,resamp_norm_z.size))
        for ii in range(0,nx):
            for jj in range(0,ny):
                for case in np.arange(0,ncases):
                    z_spline           = interp1d(norm_z[case,ii,jj,:],tke[case,ii,jj,:],kind='cubic')
                    T_spline           = interp1d(norm_z[case,ii,jj,:],T[case,ii,jj,:],kind='cubic')
                    int_tke_f[case,ii,jj,:] = z_spline(resamp_norm_z)
                    int_T_f[case,ii,jj,:]   = T_spline(resamp_norm_z)
                int_tke_f[3,ii,jj,:] = int_tke_f[1,ii,jj,:]/int_tke_f[0,ii,jj,:]
                int_tke_f[4,ii,jj,:] = int_tke_f[2,ii,jj,:]/int_tke_f[0,ii,jj,:]
            for case in np.arange(0,ncases):
                avg_int_tke[case,ii,:] = np.mean(int_tke_f[case,ii,:,:],axis=0) 
                avg_int_T[case,ii,:]   = np.mean(int_T_f[case,ii,:,:],axis=0) 
            avg_int_tke[3,ii,:] = np.mean(int_tke_f[3,ii,:,:],axis=0) 
            avg_int_tke[4,ii,:] = np.mean(int_tke_f[4,ii,:,:],axis=0) 
            stdout.write("\r%2.2f %% " % (100.0*np.float(ii+1)/np.float(nx)))
            stdout.flush()
            sleep(0.0001)
        stdout.write("\n")

        zx_resamp,xz_resamp = np.meshgrid(resamp_norm_z,x[:,0,0])

        fname_znorm = '{}_{}_normz_tke'.format(simstr[ss],terr_str)
        nz_norm = resamp_norm_z.size
        if saveData:
            fz  = open('{}{}'.format(savedir,fname_znorm),'w')
            print('{}{}'.format(savedir,fname_znorm))
            fz.write('x {}\t\tz {}\ttke_{}\t\ttke_{}\t\ttke_{}\t\ttkefrac_{}\t\ttkefrac_{}\n'.format(nx,nz_norm,cases[0],cases[1],cases[2],cases[1],cases[2]))
            for ii in range(0,nx):
                for kk in range(0,nz_norm):
                    fz.write('{0:.7f}\t{1:.2f}\t{2:.13f}\t{3:.13f}\t{4:.13f}\t{5:.13f}\t{6:.13f}\n'.format(
                        xz_resamp[ii,kk],zx_resamp[ii,kk],avg_int_tke[0,ii,kk],avg_int_tke[1,ii,kk],avg_int_tke[2,ii,kk],avg_int_tke[3,ii,kk],avg_int_tke[4,ii,kk]))    
            fz.close()

