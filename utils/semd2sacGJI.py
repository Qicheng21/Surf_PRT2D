#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:01:47 2023

@author: u1318104
"""



from obspy import read, Trace
import obspy
import numpy as np
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth
from tqdm import tqdm
from geopy import Point
from geopy.distance import distance, geodesic
import os
from multiprocessing import Pool

#%%
def semd2sac(dir_rec,net,sta,xx,comp,SOURCE='L'):
    # component X (radial) or Z
    if comp=='X':
        COMP='EHR'
    elif comp=='Z':
        COMP='EHZ'   
    else:
        print('component NOT recognized')
        return
    
        
    ttXX,aaXX=np.loadtxt(dir_rec+net+'.'+sta+'.BX'+comp+'.semd',unpack=True)
    # ttXX,aaXZ=np.loadtxt(dir_rec+net+'.'+sta+'.BXZ.semd',unpack=True)
    st=read()
    st.clear()
    st.append(Trace())
    
    if comp=='X' and SOURCE=='R':# Radial component direction
        st[0].data=-aaXX
    else:
        st[0].data=aaXX
    st[0].stats.delta=0.1
    st[0].stats.sac = obspy.core.AttribDict()
    st[0].stats.sac.delta=0.1 #ttXX[1]-ttXX[0]
    # st[0].stats.delta=0.01
    # st[0].stats.sac.delta=0.01 #ttXX[1]-ttXX[0]
    # st[0].resample(10.0)
    st[0].write(dir_rec+'/SAC/'+sta+'.'+COMP+'.sac',format='SAC')
    stX=read(dir_rec+'/SAC/'+sta+'.'+COMP+'.sac')
    
    stX[0].stats.sac.b=ttXX[0] #!!!NOT working, not saved to SAC !!!
    stX[0].stats.sac.evla=0
    stX[0].stats.sac.evlo=0
    stX[0].stats.sac.evdp=1.


    # tmpLatLonBr=distance(meters=xx).destination(Point(0, 0), 0) #lat1, lon1, Bearing
    # stX[0].stats.sac.stla=tmpLatLonBr[0]
    # stX[0].stats.sac.stlo=tmpLatLonBr[1]
    # stX[0].stats.sac.lcalda=True
    stX[0].stats.sac.dist=xx/1000.0 # m (semd) to km (sac)
    stX[0].stats.sac.lcalda=False
    
    stX[0].stats.channel=COMP

    stX[0].stats.sac.kcmpnm=COMP
    stX[0].stats.sac.kstnm=sta
    stX[0].stats.sac.knetwk=net
    stX[0].write(dir_rec+'SAC/'+sta+'.'+COMP+'.sac',format='SAC')
    # stX[0].write(dir_rec+'../tmp/'+sta+'.'+COMP+'.sac',format='SAC')

    return ttXX, aaXX, stX
#%%
# LocSL=1e6 # Left Source Location 
# LocSR=5e6 # Right Source Location

# LocSL=0.5e6 # Left Source Location 
# LocSR=5.5e6 # Right Source Location
LocSL=0.8e6 # Left Source Location 
LocSR=5.2e6 # Right Source Location

if __name__ == '__main__':
    Para_Pool=[]
    dir_SK='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/RayleighWaveGJI/'
    Stas,Network = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=str,usecols=[0,1],unpack=True)
    # XXs,ZZs = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=float,usecols=[2,3],unpack=True) #ZZs should be vertical axis, not source dis?! 
    # SOURCE='L'
    
    for SOURCE in ['L','R']:
        XXs,ZZs = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=float,usecols=[2,3],unpack=True) #ZZs should be vertical axis, not source dis?! 
        # if SOURCE=='L':
        #     XXs=abs(XXs-2000e3) #1000e3 SOURCE L, 5000e3 SOURCE R, Receivers within SOURCE L and SOURCE R
        # elif SOURCE=='R':
        #     XXs=abs(XXs-10000e3)
        if SOURCE=='L':
            XXs=abs(XXs-LocSL) #1000e3 SOURCE L, 5000e3 SOURCE R, Receivers within SOURCE L and SOURCE R
            LocS=LocSL
        elif SOURCE=='R':
            XXs=abs(XXs-LocSR)
            LocS=LocSR
            
        # for Mech in ['Mw']:#['Mw','Mw0km']:  #Focal Mechanism
        # for Mech in ['Mw0km','Mw1km','Mw5km','Mw10km','Mw20km']:
        # for Mech in ['X8000kmZ600km','X8000kmZ1000km']:
        # for Mech in ['txt']:
                # dir_rec=dir_SK+'OUTPUT_SRP_SK_SOURCEL_Layer'+str(tmpLayer)+'_dVs0.00_dVp0.00_dRho10.00_PerturbSideR/'
                # dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_'+Mech+'/'
        dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+str(int(LocS/1e3))+'km/'

        if not os.path.isdir(dir_rec+'/SAC'): #os.path.exists
            try:
                os.system('mkdir '+dir_rec+'/SAC')
            except:
                pass

        for i in range(len(Stas)):
            sta=Stas[i]
            net=Network[i]
            xx=XXs[i]
            zz=ZZs[i]
            
            Para_Pool.append((dir_rec,net,sta,xx,'X',SOURCE))
            Para_Pool.append((dir_rec,net,sta,xx,'Z',SOURCE))
        
    pool=Pool(10)
    results=pool.starmap(semd2sac,Para_Pool)
    pool.close()
    pool.join()
    print('process finished')
