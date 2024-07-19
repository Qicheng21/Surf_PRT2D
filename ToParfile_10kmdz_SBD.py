#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 15:14:24 2022

@author: u1318104
"""


#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
from cartopy.io.img_tiles import GoogleTiles, Stamen#, StamenTerrain
import matplotlib.cm as cm
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

import os
import sys
dir_script=os.path.dirname(os.path.abspath(__file__))
if not dir_script in sys.path:
    sys.path.append(dir_script)

xmax=-100
ymax=50 #36.2
xmin=-125 #-121.0
ymin=30 #32.2


#%%
def ReadEMC():
    dir_EMC=dir_script #'/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/'
    if 'VsEMC' in locals():
        # print('EMC2015 exists')
        pass
    else:
        LatsEMC, LonsEMC, DepsEMC, VsEMC=np.loadtxt(dir_EMC+'US_CrustVs_SLK_GRL_2015.txt',unpack=True)
        # print('EMC2015 NOT exists')
#%
# from PlotAll import tightgaussiansmooth

    ''' # for 2D Vs plot at 60 km
    Lats60km=LatsEMC[DepsEMC==60]
    Lons60km=LonsEMC[DepsEMC==60]
    Vs60km=VsEMC[DepsEMC==60]
    
    Lon1=np.unique(Lons60km);Lon1.sort()
    Lat1=np.unique(Lats60km);Lat1.sort()
    DATA=np.zeros((len(Lat1),len(Lon1)))
    # if scale==[]:
    #     scale=[min(data),max(data)]
        
    for k in range(len(Vs60km)):
        iLon=np.where(Lon1==Lons60km[k])
        iLat=np.where(Lat1==Lats60km[k])
        DATA[iLat,iLon]=Vs60km[k]
        
    DATA=np.ma.array(DATA,mask=~np.array(DATA,dtype=bool))
    '''
    # plt.figure()
    # plt.contourf(Lon1,Lat1,DATA,np.linspace(4.1,4.9,100),cmap=cmnew,extend='both')
    lonSRPtmp0=-112.826942
    latSRPtmp0=43.384350
    lonSRPtmp1=-114.461650
    latSRPtmp1=45.427743
    
    i0=np.argmin((LonsEMC-lonSRPtmp0)**2+(LatsEMC-latSRPtmp0)**2)
    i1=np.argmin((LonsEMC-lonSRPtmp1)**2+(LatsEMC-latSRPtmp1)**2)
    lonSRP0=LonsEMC[i0]; latSRP0=LatsEMC[i0]
    lonSRP1=LonsEMC[i1]; latSRP1=LatsEMC[i1]
    
    tmpDepSRP0=DepsEMC[np.logical_and(LonsEMC==lonSRP0,LatsEMC==latSRP0)]
    tmpVsSRP0=VsEMC[np.logical_and(LonsEMC==lonSRP0,LatsEMC==latSRP0)]
    tmpVsSRP1=VsEMC[np.logical_and(LonsEMC==lonSRP1,LatsEMC==latSRP1)]
    
    DepSRP0=np.arange(0,tmpDepSRP0[-1]+1,10)
    VsSRP0=np.zeros(DepSRP0.shape)
    VsSRP1=np.zeros(DepSRP0.shape)
    for ii in range(len(DepSRP0)-1): # averaged every 10 km
        idx1=np.where(tmpDepSRP0<=DepSRP0[ii])[0][-1]
        idx2=np.where(tmpDepSRP0>=DepSRP0[ii+1])[0][0]
        if idx1==idx2:
            idx2=idx2+1
    
        VsSRP0[ii]=np.mean(tmpVsSRP0[idx1:idx2])
        VsSRP1[ii]=np.mean(tmpVsSRP1[idx1:idx2])
    VsSRP0[-1]=tmpVsSRP0[-1]
    VsSRP1[-1]=tmpVsSRP1[-1]    
        
    # tmpDep=np.zeros(len(DepSRP0)*2)
    # tmpVs0=np.zeros(len(DepSRP0)*2)
    # tmpVs1=np.zeros(len(DepSRP0)*2)
    # for ii in range(len(DepSRP0)):
    #     if ii==0:
    #         tmpDep[ii]=DepSRP0[ii]
    #     else:
    #         tmpDep[2*ii-1]=DepSRP0[ii]
    #         tmpDep[2*ii]=DepSRP0[ii]
    #     tmpVs0[2*ii]=VsSRP0[ii]
    #     tmpVs0[2*ii+1]=VsSRP0[ii]
    #     tmpVs1[2*ii]=VsSRP1[ii]
    #     tmpVs1[2*ii+1]=VsSRP1[ii]
    # tmpDep[-1]=200
    # tmpDepSRP0=np.arange(0,600,10)
    # tmpDepSRP0=np.arange(0,4000,10)
    tmpDepSRP0=np.arange(0,1000,10)

    tmpVsSRP0=np.zeros(len(tmpDepSRP0))
    tmpVsSRP1=np.zeros(len(tmpDepSRP0))
    tmpVsSRP0[:len(VsSRP0)]=VsSRP0; tmpVsSRP0[len(VsSRP0):]=VsSRP0[-1]
    tmpVsSRP1[:len(VsSRP1)]=VsSRP1; tmpVsSRP1[len(VsSRP1):]=VsSRP1[-1]
    # return DepSRP0, VsSRP0, VsSRP1
    return tmpDepSRP0, tmpVsSRP0, tmpVsSRP1 #, DepSRP0, VsSRP0, VsSRP1


# plt.figure()
# # plt.plot(VsSRP0,DepSRP0)
# # plt.plot(VsSRP1,DepSRP0)
# plt.plot(tmpVs0,tmpDep,'.-',color='orange',label='SRP')
# plt.plot(tmpVs1,tmpDep,'.--',color='cornflowerblue',label='out')
# plt.xlabel('Vs (km/s)')
# plt.ylabel('Depth (km)')
# plt.legend()
# plt.ylim([200,0])
# plt.yticks(np.arange(0,201,20))
# plt.grid()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/03252022_SRP_10kmdxdz/Vs10kmdz.png')
#%% 1km grid for now
def Brocher_rho_vp(vstmp):
    vptmp=0.9409+(2.0947*vstmp)+(-0.8206*vstmp**2)+(0.2683*vstmp**3)+(-0.0251*vstmp**4) #Brocher 2005
    # rhotmp = 1.22679 + (1.53201*vstmp) -(0.83668*vstmp**2) + (0.20673*vstmp**3) -(0.01656*vstmp**4) # EMB source?
    rho1=1.6612*vptmp-0.4721*vptmp**2+0.0671*vptmp**3-0.0043*vptmp**4+0.000106*vptmp**5 #Brocher 2005
    return vptmp, rho1

def ToParfile(dVs=0,dVp=0,dRho=0,PLayer=120,Opt_Side='L',Opt_AbsPercent='Abs',SBD=0):
    '''
    #!!!!! Was used to modify Par_file in specfem2d for 2 laterally homogeneous model in numerical wavefield simulation, 
    #!!!!! Not needed after modification    
    
    Parameters
    ----------
    dVs : TYPE, optional
        DESCRIPTION. The default is 0.
    dVp : TYPE, optional
        DESCRIPTION. The default is 0.
    dRho : TYPE, optional
        DESCRIPTION. The default is 0.
    PLayer : TYPE, optional
        DESCRIPTION. The default is 120.
    Opt_Side : TYPE, optional
        DESCRIPTION. The default is 'L'.
    Opt_AbsPercent : TYPE, optional
        DESCRIPTION. The default is 'Abs'.
    SBD : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    '''
    # Original velocity Structure (Vs, Vp, Rho) from Schmandt et al. (2015)
    # Left side (SPECFEM2D) structure from a point inside SRP, right side (SPECFEM2D) from outside SRP ()
    # 10 km dx, dz in simulation, averaged
    # dVs - km/s if Abs, -% if Percent
    # dVp - km/s if Abs, -% if Percent
    # Rho - g/cm^3 if Abs, -% if Percent  
    # PLayer perturb Layer id 1-60 1-bottom 60-top
    # SBD: Smooth Boundary, grid points, how many blocks, not in km
    tmpVs=0
    tmpDep0=0
    tmpk=0
    Material0=[]
    
    Nz=120

    # Nz=400
    # Nz=100
    # Nz=120
    #Smooth Boundary
    BD0=600 #Change if needed, dependent on the nx, number of grids in x direction
    # BD0=400 #Change if needed, dependent on the nx, number of grids in x direction
    # BD0=300 #Change if needed, dependent on the nx, number of grids in x direction

    BDL=BD0-int(SBD/2)+1#276
    BDR=BD0+int(SBD/2)#325

    DepSRP0, VsSRP0, VsSRP1=ReadEMC()
    Region=[]
    for ii in range(len(DepSRP0)):
        # if tmpVs==VsSRP0[ii]:
        #     continue
        tmpk+=1
        tmpVs=VsSRP0[ii]
        tmpVsSK=tmpVs.copy()
        tmpVpSK, tmpRhoSK=Brocher_rho_vp(tmpVsSK)
        # Perturbation for Numerical Sensitivity Kernels
        if Opt_Side=='L' and ii==Nz-PLayer:
            if Opt_AbsPercent=='Abs':
                tmpVsSK=tmpVsSK+dVs
                tmpVpSK=tmpVpSK+dVp
                tmpRhoSK=tmpRhoSK+dRho
            else:
                tmpVsSK=tmpVsSK*(1+dVs/100.)
                tmpVpSK=tmpVpSK*(1+dVp/100.)
                tmpRhoSK=tmpRhoSK*(1+dRho/100.)
            
        Material0.append('%d 1 %d.d0 %d.d0 %d.d0 0 0 9999 9999 0 0 0 0 0 0'%(tmpk,round(tmpRhoSK*1e3),round(tmpVpSK*1e3),round(tmpVsSK*1e3)) )
        
        if ii==0:
            continue
        tmpDep1=DepSRP0[ii]
        # Region.append('%d %d %d %d %d'%(1,300,60-int(tmpDep1/10)+1,60-int(tmpDep0/10),tmpk-1))
        # Region.append('%d %d %d %d %d'%(1,100,60-int(tmpDep1/10)+1,60-int(tmpDep0/10),tmpk-1))
        Region.append('%d %d %d %d %d'%(1,BDL-1,Nz-int(tmpDep1/10)+1,Nz-int(tmpDep0/10),tmpk-1))

        tmpDep0=tmpDep1
    # Region.append('%d %d %d %d %d'%(1,300,1,60-int(tmpDep0/10),tmpk))
    # Region.append('%d %d %d %d %d'%(1,100,1,60-int(tmpDep0/10),tmpk))
    Region.append('%d %d %d %d %d'%(1,BDL-1,1,Nz-int(tmpDep0/10),tmpk))


    
        
    tmpVs=0
    tmpDep0=0
    Material1=[]
    for ii in range(len(DepSRP0)):
        # if tmpVs==VsSRP1[ii]:
        #     continue
        tmpk+=1
        tmpVs=VsSRP1[ii]
        tmpVsSK=tmpVs.copy()
        tmpVpSK, tmpRhoSK=Brocher_rho_vp(tmpVsSK)
        # Perturbation for Numerical Sensitivity Kernels
        if Opt_Side=='R' and ii==Nz-PLayer:
            if Opt_AbsPercent=='Abs':
                tmpVsSK=tmpVsSK+dVs
                tmpVpSK=tmpVpSK+dVp
                tmpRhoSK=tmpRhoSK+dRho
            else:
                tmpVsSK=tmpVsSK*(1+dVs/100.)
                tmpVpSK=tmpVpSK*(1+dVp/100.)
                tmpRhoSK=tmpRhoSK*(1+dRho/100.)
            
        Material1.append('%d 1 %d.d0 %d.d0 %d.d0 0 0 9999 9999 0 0 0 0 0 0'%(tmpk,round(tmpRhoSK*1e3),round(tmpVpSK*1e3),round(tmpVsSK*1e3)) )
        
        if ii==0:
            continue
        tmpDep1=DepSRP0[ii]
        # Region.append('%d %d %d %d %d'%(301,600,60-int(tmpDep1/10)+1,60-int(tmpDep0/10),tmpk-1))
        # Region.append('%d %d %d %d %d'%(101,200,60-int(tmpDep1/10)+1,60-int(tmpDep0/10),tmpk-1))
        Region.append('%d %d %d %d %d'%(BDR+1,BD0*2,Nz-int(tmpDep1/10)+1,Nz-int(tmpDep0/10),tmpk-1))

        tmpDep0=tmpDep1
    # Region.append('%d %d %d %d %d'%(301,600,1,60-int(tmpDep0/10),tmpk))
    # Region.append('%d %d %d %d %d'%(101,200,1,60-int(tmpDep0/10),tmpk))
    Region.append('%d %d %d %d %d'%(BDR+1,BD0*2,1,Nz-int(tmpDep0/10),tmpk))

    #Smooth Boundary transition
    if SBD==0:
        pass
    else:
        flag_end=0
        for ii in range(len(DepSRP0)):
            vsL=VsSRP0[ii]
            vsR=VsSRP1[ii]
            if vsL==VsSRP0[-1]:
                flag_end=1
            for xj in range(BDL,BDR+1):
                tmpk+=1
                tmpVs=vsL+(vsR-vsL)*(1-np.cos((xj-BDL+0.5)/(BDR-BDL+1)*np.pi))/2
                tmpVp, tmpRho=Brocher_rho_vp(tmpVs)
                Material1.append('%d 1 %d.d0 %d.d0 %d.d0 0 0 9999 9999 0 0 0 0 0 0'%(tmpk,round(tmpRho*1e3),round(tmpVp*1e3),round(tmpVs*1e3)) )
                if flag_end:
                    Region.append('%d %d %d %d %d'%(xj,xj,1,Nz-ii,tmpk))
                else:
                    Region.append('%d %d %d %d %d'%(xj,xj,Nz-ii,Nz-ii,tmpk))
            
            if flag_end:
                break
        
    #%
    dir_SRP='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_LocWid/DATA/'
    Par=[]
    Start_MatReg=0
    End_MatReg=0
    MatReg=0 # 0 Material, 1 Region
    
    # with open(dir_SRP+'Ori_Par_file_10kmdxdz') as fp:
    with open(dir_SRP+'Ori_Par_file_X12000kmZ1200kmS10km') as fp:
        for tmp_line in fp:
            line=tmp_line.split('\n')[0]
            if 'nbmodels' in line and '=' in line:
                Par.append('nbmodels\t\t= %d'%(len(Material0)+len(Material1)))
                continue
            if 'nbregions' in line and '=' in line:
                Par.append('nbregions\t\t= %d'%(len(Region)))
                continue
            
            if line=='# Start_Material':
                Start_MatReg=1
                MatReg=0
            if line=='# End_Material':
                End_MatReg=1
            if line=='# Start_Region':
                Start_MatReg=1
                MatReg=1
            if line=='# End_Region':
                End_MatReg=1
                
            if not Start_MatReg: # General case, just append
                Par.append(line)
            elif not End_MatReg: # Start Material or Region, haven't ended yet
                continue
            else: #End Material or Region
                if not MatReg: 
                    Par=Par+['# Start_Material']+Material0+Material1+['# End_Material']
                else:
                    Par=Par+['# Start_Region']+Region+['# End_Region']
                Start_MatReg=0
                End_MatReg=0
    
    np.savetxt(dir_SRP+'/Par_file',Par,fmt='%s')
    
#%%
if __name__ == '__main__':
    ToParfile(SBD=0)
    # DepSRP0, VsSRP0, VsSRP1=ReadEMC()
    pass

