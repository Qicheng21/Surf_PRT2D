#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 12:02:29 2022

@author: u1318104
"""

# sys.path.append('')
import sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi')
# from SWRT_Zeng_ex8_LockedMode import EigenRead, single_freq
# from SWRT_Zeng_ex13_LockedDeep import EigenRead, single_freq, Csingle_freq, Csingle_freq2
from SWRT_Zeng_ex14_LockedModes import EigenRead, single_freq, Csingle_freq, Csingle_freq2, Csingle_freqA
from obspy import read, Trace, Stream
import numpy as np
import matplotlib.pyplot as plt

from numpy import cos, pi, ceil, exp
from math import log
from scipy.fft import fft, ifft
from scipy.signal import hilbert
from scipy import interpolate
from scipy.integrate import simps

from multiprocessing import Pool
import os
from tqdm import tqdm
#%% Parameters and Functions
dir_CPS='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/test_CPS330/'
# dirM=dir_CPS+'ex8_ModeSumWave'
# dirM=dir_CPS+'ex9_ModeSumLock'
# dirM=dir_CPS+'ex11_ModeSumUB0'
# dirM=dir_CPS+'ex7_SRP'
# dirM=dir_CPS+'ex8_ModeSumWave'
# dirM=dir_CPS+'ex12B_Vs6p9_Dep900'

dirM=dir_CPS+'ex14A_Vs6p9_Dep900_M10/'


eps=1e-5

# else:
X_SL=1000 #Left Souce Location
X_SR=5000 #Right Source Location
X_C=3000 #Center or structural location

X_RL=1000 #Leftmost receiver location
X_RR=5000 #Rightmost receiver location
# X_SL=2000 #Left Souce Location
# X_SR=6000 #Right Source Location
# X_C=4000 #Center or structural location

# X_RL=2000 #Leftmost receiver location
# X_RR=6000 #Rightmost receiver location


distC=X_C-X_RL

def TimeShiftFre(xt,tt,t0):
    #positive t0 - advance
    #negative t0 - delay
    N=int(2**ceil(np.log2(len(xt)*2)))
    dt=tt[1]-tt[0]
    df=1/(N*dt)
    ff=np.arange(0,N*df,df)
    
    xtpad=np.zeros(N)
    xtpad[0:len(xt)]=xt
    Xfpad=fft(xtpad)
    Coe=np.exp(2j*pi*ff*t0)
    XfNewpad=Xfpad*Coe
    XfNewpad[int(N/2)+1:]=np.conjugate(XfNewpad[int(N/2-1):0:-1]) #conjugate symmetry
    XfNewpad[int(N/2)]=XfNewpad[int(N/2)].real #Midpoint, need to be real? 
    xtNewpad=np.real(ifft(XfNewpad))
    xtNew=xtNewpad[0:len(xt)]
    
    return xtNew # ,ff,Xfpad,XfNewpad,N

# def TimeShiftFreDispersion(xt,tt,t0Arr,rcoefArr):
def TaperCos(xt,tt,tb,te,Perc=0.01):
    # Taper using cosine function, only works for single trace,multiple traces - matrix multiplication TBD
    # xt - waveform/signal
    # tt - time/sampling domain
    # tb - start of tapering
    # te - end of tapering
    # Perc - Percentage of taper on one side, <0.5
    xt_taper=np.zeros(xt.shape)
    Ttaper=(te-tb)*Perc #time length of taper on one side
    for itmp, tmp in enumerate(tt):
        if tmp>=tb and tmp<tb+Ttaper:
            xt_taper[itmp]=xt[itmp]*(1-np.cos(np.pi*(tmp-tb)/Ttaper))/2
        elif tmp>=tb+Ttaper and tmp<=te-Ttaper:
            xt_taper[itmp]=xt[itmp]
        elif tmp>te-Ttaper and tmp<=te:
            xt_taper[itmp]=xt[itmp]*(1-np.cos(np.pi*(te-tmp)/Ttaper))/2
    return xt_taper

def TaperCos1(xt,tt,tb,te,Perc=0.1,Pts=60):
    # Taper using cosine function, only works for single trace,multiple traces - matrix multiplication TBD
    # xt - waveform/signal
    # tt - time/sampling domain
    # tb - start of tapering
    # te - end of tapering
    # Perc - Percentage of taper on one side
    # Pts - Points of taper on one side
    xt_taper=np.zeros(xt.shape)
    Ttaper=(te-tb)*Perc #time length of taper on one side
    for itmp, tmp in enumerate(tt):
        # if tmp>=tb-Ttaper and tmp<tb:
        #     xt_taper[itmp]=xt[np.argmin(np.abs(tt-tb))]*(1-np.cos(np.pi*(Ttaper+tmp-tb)/Ttaper))/2
        # if tmp>0 and tmp<tb: #Taper All the way to 0 Hz
        #     xt_taper[itmp]=xt[np.argmin(np.abs(tt-tb))]*(1-np.cos(np.pi*tmp/tb))/2
        
        if tmp>0 and tmp<tb: #Taper All the way to 0 Hz
            xt_taper[itmp]=xt[np.argmin(np.abs(tt-tb))]*(1-np.cos(np.pi*tmp/tb))/2
        elif tmp>=tb and tmp<=te:
            xt_taper[itmp]=xt[itmp]
        elif tmp>te and tmp<=te+Ttaper:
            xt_taper[itmp]=xt[np.argmin(np.abs(tt-te))]*(1-np.cos(np.pi*(Ttaper+te-tmp)/Ttaper))/2
    return xt_taper

def TaperCosM(xt,tt,xt0,tt0,Pts=60,Perc=0.1):
    xt_taper=np.zeros(xt.shape)
    Ttaper=(tt[1]-tt[0])*Pts
    tb=tt0[0]
    te=tt0[-1]
    for itmp, tmp in enumerate(tt):
        if tmp>=tb-Ttaper and tmp<tb:
            xt_taper[itmp]=xt0[0]*(1-np.cos(np.pi*(Ttaper+tmp-tb)/Ttaper))/2
        elif tmp>=tb and tmp<=te:
            xt_taper[itmp]=xt[itmp]
        elif tmp>te and tmp<=te+Ttaper:
            xt_taper[itmp]=xt0[-1]*(1-np.cos(np.pi*(Ttaper+te-tmp)/Ttaper))/2
    return xt_taper
    
def ToTrace(data,delta):
    tr=Trace()
    tr.data=data
    tr.stats.delta=delta
    return tr

def TimeShiftDisper_old(xt,tt,Freqs,t0Arr,r0Arr):
    # positive t0 - advance
    # negative t0 - delay
    # Freqs, must be consistent with t0Arr and r0Arr
    # t0Arr, time shift for each frequency (Freqs)
    # r0Arr, amplitude change for each frequency (Freqs)
    
    N=int(2**ceil(np.log2(len(xt)*2)))
    dt=tt[1]-tt[0]
    df=1/(N*dt)
    ff=np.arange(0,N*df,df)
    
    tck_t0=interpolate.splrep(Freqs,t0Arr)
    tck_r0=interpolate.splrep(Freqs,r0Arr)
    
    t0Arr_interp=interpolate.splev(ff,tck_t0); # works well inside Freqs, need filtering afterwards
    r0Arr_interp=interpolate.splev(ff,tck_r0); 
    # Filtering may not work as expected, 
    t0Arr_interp=TaperCos(t0Arr_interp,ff,np.min(Freqs)/1.2,np.max(Freqs)/0.8) #Taper in spectrum, not necessary, 
    r0Arr_interp=TaperCos(r0Arr_interp,ff,np.min(Freqs)/1.2,np.max(Freqs)/0.8) #filtering afterwards does the same work

    xtpad=np.zeros(N)
    xtpad[0:len(xt)]=xt
    Xfpad=fft(xtpad)
    Coe=r0Arr_interp*np.exp(2j*pi*ff*t0Arr_interp)
    XfNewpad=Xfpad*Coe
    XfNewpad[int(N/2)+1:]=np.conjugate(XfNewpad[int(N/2-1):0:-1]) #conjugate symmetry
    XfNewpad[int(N/2)]=XfNewpad[int(N/2)].real #Midpoint, need to be real?                                           
    xtNewpad=np.real(ifft(XfNewpad))
    xtNew=xtNewpad[0:len(xt)]
    
    return xtNew

def TimeShiftDisper(xt,tt,Freqs,t0Arr,r0Arr,tmpMode):
    # positive t0 - advance
    # negative t0 - delay
    # Freqs, must be consistent with t0Arr and r0Arr
    # t0Arr, time shift for each frequency (Freqs)
    # r0Arr, amplitude change for each frequency (Freqs)
    
    N=int(2**ceil(np.log2(len(xt)*2)))
    dt=tt[1]-tt[0]
    df=1/(N*dt)
    ff=np.arange(0,N*df,df)
    
    # xb=np.max([Freqs[~np.isnan(t0Arr)][0],Freqs[~np.isnan(r0Arr)][0]]) #should be the same, just in case
    # xe=np.min([Freqs[~np.isnan(t0Arr)][-1],Freqs[~np.isnan(r0Arr)][-1]])
    
    t0Arr1=t0Arr[~np.isnan(t0Arr)]
    Freqs1t=Freqs[~np.isnan(t0Arr)]
    r0Arr1=r0Arr[~np.isnan(r0Arr)]
    Freqs1r=Freqs[~np.isnan(r0Arr)]
    
    if len(Freqs1t)<=5: #TBD, not enough values, due to cut-off freq for higher modes
        print('Mode %d Not enough points for interpolation'%tmpMode)
        return np.zeros(xt.shape) 
    
    tck_t0=interpolate.splrep(Freqs1t,t0Arr1)
    tck_r0=interpolate.splrep(Freqs1r,r0Arr1)
        
    t0Arr_interp=interpolate.splev(ff,tck_t0); # works well inside Freqs, need filtering afterwards
    r0Arr_interp=interpolate.splev(ff,tck_r0); 
    # Filtering may not work as expected, tapering in frequency domain first would help. 
    t0Arr_interp=TaperCos(t0Arr_interp,ff,np.min(Freqs1t)/1.2,np.max(Freqs1t)/0.8) #Taper in spectrum, not necessary, 
    r0Arr_interp=TaperCos(r0Arr_interp,ff,np.min(Freqs1r)/1.2,np.max(Freqs1r)/0.8) #filtering afterwards does the same work
    
    xtpad=np.zeros(N)
    xtpad[0:len(xt)]=xt
    Xfpad=fft(xtpad)
    Coe=r0Arr_interp*np.exp(2j*pi*ff*t0Arr_interp)
    XfNewpad=Xfpad*Coe
    XfNewpad[int(N/2)+1:]=np.conjugate(XfNewpad[int(N/2-1):0:-1]) #conjugate symmetry
    XfNewpad[int(N/2)]=XfNewpad[int(N/2)].real #Midpoint, need to be real?                                           
    xtNewpad=np.real(ifft(XfNewpad))
    xtNew=xtNewpad[0:len(xt)]
    
    return xtNew


#%%
effileL=dirM+'/SRDER_L.TXT'
effileR=dirM+'/SRDER_R.TXT'
# effileL=dirM+'/SRDER_R.TXT'
# effileR=dirM+'/SRDER_L.TXT'
per0=150
per1=10
dp=-1
Pers=np.arange(per0,per1+dp,dp)
Freqs=1/Pers
NModes=10
# # Freqs=np.arange(freq0,freq1+df,df)
tcm=np.empty((NModes,len(Freqs)));tcm[:]=np.nan
rcm=np.empty((NModes,len(Freqs)));rcm[:]=np.nan

Tcm=np.empty((NModes,NModes,len(Freqs)));Tcm[:]=np.nan
Rcm=np.empty((NModes,NModes,len(Freqs)));Rcm[:]=np.nan

cArrL=np.empty((NModes,len(Freqs)));cArrL[:]=np.nan
cArrR=np.empty((NModes,len(Freqs)));cArrR[:]=np.nan

UzL=np.empty((NModes,200,len(Freqs)));UzL[:]=np.nan #200 -> Number of layers in Eigenfunction
UrL=UzL.copy();UrL[:]=np.nan
UzR=UzL.copy();UzR[:]=np.nan
UrR=UzL.copy();UrR[:]=np.nan

for idf, freq in tqdm(enumerate(Freqs)):
    # print(idf,freq)
    erL=EigenRead(freq,effileL)
    erR=EigenRead(freq,effileR)
    cArrL[:len(erL.Vph),idf]=erL.Vph#float(erL.Vph) #[0] for fundamental mode
    cArrR[:len(erR.Vph),idf]=erR.Vph#float(erR.Vph) #[0]
    # rcoef,tcoef=single_freq(freq,effileL,effileR)
    Csf=Csingle_freqA(freq,effileL,effileR)
    rcoef=Csf.rcoef
    tcoef=Csf.tcoef
    
    tcm[0:len(tcoef),idf]=tcoef.reshape(-1)
    rcm[0:len(rcoef),idf]=rcoef.reshape(-1)
    
    Tcm[0:Csf.Tcoef.shape[1],0:Csf.Tcoef.shape[0],idf]=Csf.Tcoef.T
    Rcm[0:Csf.Rcoef.shape[1],0:Csf.Rcoef.shape[0],idf]=Csf.Rcoef.T
    
    UzL[:erL.Nmed,:,idf]=erL.r2
    UrL[:erL.Nmed,:,idf]=erL.r1
    UzR[:erR.Nmed,:,idf]=erR.r2
    UrR[:erR.Nmed,:,idf]=erR.r1

tcm=Tcm[0]
rcm=Rcm[0]
#%%
# tmpn=1
# # plt.figure()
# # plt.imshow(UzL[tmpn])
# # plt.colorbar()

# plt.figure()
# plt.contourf(Freqs,zz,UzL[tmpn])
# plt.ylim([zz[-1],zz[0]])
# plt.colorbar()
#%%

# erL=EigenRead(1/52,effileL)
# erR=EigenRead(1/52,effileR)
# # CsfL=Csingle_freq(freq,effileL,effileR)
# CsfL=Csingle_freq(freq,effileR,effileL)


# u,d,vh=np.linalg.svd(CsfL.LHS,full_matrices=True)


#%%
# tmp_m=1
# plt.figure()
# plt.plot(erL.Deps,erL.Psi_xx[tmp_m,:],'r--^')
# plt.plot(erR.Deps,erR.Psi_xx[tmp_m,:],'b--^')

# plt.plot(erL.Deps,erL.r1[tmp_m,:],'r--^')
# plt.plot(erR.Deps,erR.r1[tmp_m,:],'b--^')

# plt.plot(erL.Deps,erL.Lambdas,'r--^')
# plt.plot(erR.Deps,erR.Lambdas,'b--^')

# plt.plot(erL.Deps,erL.Mus,'r--^')
# plt.plot(erR.Deps,erR.Mus,'b--^')

#%%
plt.close('all')

plt.figure()
for i in range(NModes):
   	mnum='Mode %d' %(i)
   	plt.plot(Freqs,cArrL[i],'o-',label=mnum)
plt.legend()
plt.ylabel('Vph (km/s)')
plt.xlabel('Freq (Hz)')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/10282022_LockedSummary/ModesL_Vph.png')

plt.figure()
for i in range(NModes):
   	mnum='Mode %d' %(i)
   	plt.plot(Freqs,rcm[i],'o-',label=mnum)
plt.legend()
# plt.ylabel('Vph (km/s)')
plt.xlabel('Freq (Hz)')
plt.title('Ref Coef')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/10282022_LockedSummary/ModesL_RefCoef.png')

plt.figure()
for i in range(NModes):
   	mnum='Mode %d' %(i)
   	plt.plot(Freqs,tcm[i],'o-',label=mnum)
plt.legend()
# plt.ylabel('Vph (km/s)')
plt.xlabel('Freq (Hz)')
plt.title('Trans Coef')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/10282022_LockedSummary/ModesL_TransCoef.png')

#%%
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
tmpSBD=0
SOURCE='L'
# SOURCE='R'
COMP='EHZ'
# Mode='0'
Mode='All'
dir_SK='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_LocWid/'
# dir_rec=dir_SK+'OUTPUT_SRP_SK_SOURCE'+SOURCE+'_SBD'+str(tmpSBD)+'km_WithRef/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Mode'+Mode+'_WithRef/'
dir_SPECFEM='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_LocWid/OUTPUT_SRP_SOURCEL_X6000kmZ600km_v2_44CPU/SAC/'

# stT=read(dir_SK+'OUTPUT_SRP_SOURCEL_Mw_WithRef/SAC/S1901.EHZ.sac')
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_SPCETmpMode'+Mode+'_WithRef/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_LockUB0_Mode'+Mode+'_WithRef/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_LockUB0TmpSPEC_Mode'+Mode+'_WithRef/'
suffix=dirM.split('/')[-1]
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Lock_'+suffix+'_Mode'+Mode+'_WithRef/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Lock_'+suffix+'_Mode'+Mode+'_WithRef_test0_10/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Lock_'+suffix+'_Mode'+Mode+'_WithRef_test0_10_RTcm/'
dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Lock_'+suffix+'_Mode'+Mode+'_WithRef0_10/'
# dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_Lock_'+suffix+'_Mode'+Mode+'_WithRef10_10_test4_nART_Taper/'



if not os.path.isdir(dir_rec):
    os.makedirs(dir_rec)
if not os.path.isdir(dir_rec+'/SAC'):
    os.makedirs(dir_rec+'/SAC')
    


Stas,Network = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=str,usecols=[0,1],unpack=True)
Dists=np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=float,usecols=[2],unpack=True)
Dists=Dists/1000-X_SL #m to km, substract source location


# dir_NormalModes='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/test_CPS330/ex11_ModeSumUB0/NormalModes/'
dir_NormalModes='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/test_CPS330/ex14A_Vs6p9_Dep900_M10/NormalModes/'

# dir_NormalModesL='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_LocWid/OUTPUT_SRP_SK_SOURCEL_SBD0km/'
Stas,Network = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=str,usecols=[0,1],unpack=True)
def ToSACModes(dir_rec,sta,dist1,SOURCE,Mode):
    st=Stream()
    for m in range(NModes):
        # if m==0:
        #     if SOURCE=='L':
        #         st+=read(dir_NormalModes+'S1501.EHZ.sac_1km')
        #     else:
        #         st+=read(dir_NormalModes+'S2501.EHZ.sac_1km')
        #     continue
        st+=read(dir_NormalModes+'150000100.ZEX_M'+str(m)+'_'+SOURCE)
        
        # st+=read(dir_NormalModes+'150000100.ZEX_M'+str(m)+'_'+SOURCE)
        # if SOURCE=='L':
        #     st+=read(dir_NormalModes+'S1501.EHZ.sac_1km')
        # else:
        #     st+=read(dir_NormalModes+'S2501.EHZ.sac_1km')
    
    uz=st[0].copy()
    Uz=st[0].copy()
    Uz.stats.sac.dist=dist1
    tb=uz.stats.sac.b
    te=uz.stats.sac.e
    dt=uz.stats.sac.delta
    tt=np.arange(tb,te+dt,dt)
    dist=uz.stats.sac.dist
    vmin=2.5
    vmax=5.
    tb_inc=dist/vmax
    te_inc=dist/vmin
    
    st1=st.copy()
    st1_inc=st.copy()
    st1_ref=st.copy()
    Uz.data=np.zeros(Uz.data.shape)
    for n in range(NModes):
        for m in range(NModes):                          
            st1[m].stats.sac.dist=dist1
            st1_inc[m].stats.sac.dist=dist1
            st1_ref[m].stats.sac.dist=dist1
            
            # data_inc=TaperCos(st[m].data,tt,tb_inc,te_inc)
            # if dist1<=distC: #behaviour exactly at boundary? 
            #     st1_inc[m].data=TimeShiftDisper(data_inc, tt, 
            #                                   Freqs,-(dist1-dist)/cArrL[m],np.ones(Freqs.shape),m)
            #     st1_ref[m].data=TimeShiftDisper(st1_inc[m].data,tt,Freqs,-2*(distC-dist1)/cArrL[m],-rcm[m],m)
            # else:
            #     st1_inc[m].data=TimeShiftDisper(data_inc, tt, 
            #                                   Freqs,-(distC-dist)/cArrL[m]-(dist1-distC)/cArrR[m],tcm[m],m)
            #     st1_ref[m].data=np.zeros(data_inc.shape)
            
            # data_inc=TaperCos(st[n].data,tt,tb_inc,te_inc)
            data_inc=st[n].data.copy()
            # data_inc=TaperCos(st[n].data,tt,tb_inc,te_inc)
            if dist1<=distC: #behaviour exactly at boundary? 
                st1_inc[m].data=TimeShiftDisper(data_inc, tt, 
                                              Freqs,-(dist1-dist)/cArrL[n],np.ones(Freqs.shape),m)
                st1_ref[m].data=TimeShiftDisper(st1_inc[m].data,tt,Freqs,-(distC-dist1)/cArrL[n]-(distC-dist1)/cArrL[m],-Rcm[n][m],m)
                if m: #high modes, before contrast, no incident wave
                    st1_inc[m].data=np.zeros(data_inc.shape)
            else: #transmitted wave actually
                st1_inc[m].data=TimeShiftDisper(data_inc, tt, 
                                              Freqs,-(distC-dist)/cArrL[n]-(dist1-distC)/cArrR[m],Tcm[n][m],m)
                st1_ref[m].data=np.zeros(data_inc.shape)
                
            st1[m].data=st1_inc[m].data+st1_ref[m].data
            
        # break
        
        if Mode=='All':
            nSum=NModes
        else:
            nSum=int(Mode)+1

        uz.stats.sac.dist=dist1
        uz.data=np.zeros(uz.data.shape)
        for m in range(nSum):
            uz.data+=st1[m].data
                
        Uz.data+=uz.data
        if n==0:
            uz0=uz.copy()
            break #fundamental incident wave!!!!

        
    Uz.write(dir_rec+'/SAC/'+sta+'.'+COMP+'.sac',format='SAC') #!!! Writing to disk
    # return tt, uz,data_inc,dist1,dist,st
    return tt, uz, st1_inc, st1_ref, st1, data_inc,st

#%%
Para_Pool=[]
for idx in range(len(Stas)):
    # dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_SPCETmpMode'+Mode+'_WithRef/'
    sta=Stas[idx]
    dist1=Dists[idx]
    if SOURCE =='R':
        dist1=4e3-dist1
    Para_Pool.append((dir_rec,sta,dist1,SOURCE,Mode))

#%%
pool=Pool(24)
results=pool.starmap(ToSACModes,Para_Pool)
pool.close()
pool.join()
print('process finished')


#%%
SOURCE='L'
m=0
# st=read(dir_NormalModes+'150000100.ZEX_M'+str(m)+'_'+SOURCE)
# uz=st[0]
# xt=uz.data

# tb=uz.stats.sac.b
# te=uz.stats.sac.e
# dt=uz.stats.sac.delta
# tt=np.arange(tb,te+dt,dt)

st=read(dir_NormalModes+'S1501.EHZ.sac')
uz=st[0]
dist=uz.stats.sac.dist
Nt=uz.stats.npts
vmin=2.5
vmax=5.
tb_inc=dist/vmax
te_inc=dist/vmin   
tb=uz.stats.sac.b
te=uz.stats.sac.e
dt=0.1
tt=np.arange(tb,tb+Nt*dt,dt)

xt=TaperCos1(uz.data,tt,tb_inc,te_inc) #Perc=0.05

dcZ=np.sum(xt)
dcR=np.sum(TaperCos1(read(dir_NormalModes+'S1501.EHR.sac')[0].data,tt,tb_inc,te_inc))

plt.figure()
plt.plot(tt,uz.data,'k')
plt.plot(tt,xt,'r--')
#%%
N=int(2**ceil(np.log2(len(xt)*2)))
dt=0.1 # dt=tt[1]-tt[0]
tt1=np.arange(0,N*dt,dt)
df=1/(N*dt)
ff=np.arange(0,N*df,df)

xx=np.arange(0,6e6+1,2.5e3)
zz=erL.Deps*1e3
x=3000 #km



xtpad=np.zeros(N)
xtpad[0:len(xt)]=xt
Xfpad=fft(xtpad)

#%%
def Interp_Taper(Ff,Freqs,ff,Pts=60,Perc=0.1):
    # interpolate from tt to tt1, with Taper at edge, exclude nan values
    # Ff - 1D array, could contain nan values
    # Freqs - original sampling
    # ff- new sampling
    
    Ffn=Ff[~np.isnan(Ff)]
    Freqsn=Freqs[~np.isnan(Ff)]
    if len(Freqsn)<=5:
        print('Warning: Not eough points for interpolation!')
        return np.zeros(ff.shape)
    
    tck=interpolate.splrep(Freqsn,Ffn)
    Ff_interp=interpolate.splev(ff,tck)
    # Ff1=TaperCos1(Ff_interp,ff,np.min(Freqsn),np.max(Freqsn))
    Ff1=TaperCosM(Ff_interp,ff,Ffn,Freqsn,Pts=60,Perc=0.1)
    
    return Ff1
#%%
UzL1=np.zeros((NModes,200,len(ff)))#200 -> Number of layers in Eigenfunction
UrL1=UzL1.copy()
UzR1=UzL1.copy()
UrR1=UzL1.copy()
for m in tqdm(range(NModes)):
    # if m !=1:
    #     continue
    for iz in tqdm(range(len(zz))):
        
        UzL1[m,iz,:]=Interp_Taper(UzL[m,iz,:],Freqs,ff)
        UrL1[m,iz,:]=Interp_Taper(UrL[m,iz,:],Freqs,ff)
        UzR1[m,iz,:]=Interp_Taper(UzR[m,iz,:],Freqs,ff)
        UrR1[m,iz,:]=Interp_Taper(UrR[m,iz,:],Freqs,ff)

    # break
#%%
tmpn=7
plt.figure()
plt.imshow(UzL[tmpn],aspect='auto')
# plt.xlim([0,])
plt.colorbar()

plt.figure()
plt.imshow(UzL1[tmpn],aspect='auto')
plt.xlim([0,700])
plt.colorbar()
#%%
iz=70
plt.figure()
plt.plot(Freqs,UzL[tmpn,iz],'k')
plt.plot(ff,UzL1[tmpn,iz],'r--')
#%%
# plt.figure()
# plt.plot(Freqs,t_delay,'k')
# plt.plot(ff,t_delay1,'r--')
#%%
dir_fp=dirM
t0=740 #s
x_meters=2800000
# WF=np.zeros((len(zz),len(xx)))
m=0 #Mode 0,1,...,9
TypeW='Inc' #'Inc' for incident (transmitted) or 'Ref' for reflected




# for jx,x in tqdm(enumerate(xx)):
#######################################################################################
name_fp='WF_'+str(t0)+'s'+'_Mode'+str(m)
x=x_meters/1000 # m to km
x0=2500 #km

############
# x0=2500 #km
dist=1500 #km x0-X_SL
dist1=x-X_SL # input x in meters
if dist1<=distC: #behaviour exactly at boundary? 
    td_inc=-(dist1-dist)/cArrL[0]
    td_ref=-(distC-dist)/cArrL[0]-(distC-dist1)/cArrL[m]
    
    r0_inc=np.ones(Freqs.shape)
    r0_ref=-rcm[m]
    
    if m: #higher modes, before contrast, no incident wave
        r0_inc=np.zeros(Freqs.shape)
else:
    td_inc=-(distC-dist)/cArrL[0]-(dist1-distC)/cArrR[m]
    td_ref=np.zeros(Freqs.shape)
    
    r0_inc=tcm[m]
    r0_ref=np.zeros(Freqs.shape) # Transmitted side, no reflected wave

r0_inc1=Interp_Taper(r0_inc,Freqs,ff)
td_inc1=Interp_Taper(td_inc,Freqs,ff)
r0_ref1=Interp_Taper(r0_ref,Freqs,ff)
td_ref1=Interp_Taper(td_ref,Freqs,ff)
Xfpad1_inc=r0_inc1*Xfpad*np.exp(2j*pi*ff*td_inc1) 
Xfpad1_ref=r0_ref1*Xfpad*np.exp(2j*pi*ff*td_ref1) 

if dist1<=distC:
    Uz1=UzL1[m]
    Ur1=UrL1[m]
else:
    Uz1=UzR1[m]
    Ur1=UrR1[m]
    
Xfpad1Z_inc=Uz1*Xfpad1_inc.reshape((1,-1));
Xfpad1R_inc=-1j*Ur1*Xfpad1_inc.reshape((1,-1));

Xfpad1Z_ref=Uz1*Xfpad1_ref.reshape((1,-1));
Xfpad1R_ref=1j*Ur1*Xfpad1_ref.reshape((1,-1)); # reflection, propagation direction reversed

Xfpad1Z_inc[:,0]=0;#dcZ
Xfpad1R_inc[:,0]=0;#dcR
Xfpad1Z_ref[:,0]=0;#dcZ
Xfpad1R_ref[:,0]=0;#dcR

### Snapshot
def Xfpad2wf(Xfpad1,ff,t0,dt,zz):
    N=len(ff)
    Xfpad1_t0=Xfpad1*exp(1j*2*pi*ff.reshape((1,-1))*t0)
    Xfpad1_t0[:,int(N/2)+1:]=np.conjugate(Xfpad1_t0[:,int(N/2-1):0:-1]) #conjugate symmetry
    Xfpad1_t0[:,int(N/2)]=Xfpad1_t0[:,int(N/2)].real #Midpoint, need to be real?
    wf=simps(Xfpad1_t0,dt*ff,axis=-1)
    wf=wf.real
    return wf

wfZ_inc=Xfpad2wf(Xfpad1Z_inc,ff,t0,dt,zz)
wfR_inc=Xfpad2wf(Xfpad1R_inc,ff,t0,dt,zz)
wfZ_ref=Xfpad2wf(Xfpad1Z_ref,ff,t0,dt,zz)
wfR_ref=Xfpad2wf(Xfpad1R_ref,ff,t0,dt,zz)

# with open (dir_fp+name_fp,'a+') as fp:
#     for iz,z in enumerate(zz):
#         fp.write('%f %f %.8E %.8E %.8E %.8E\n'%(z,x_meters,wfZ_inc[iz],wfR_inc[iz],wfZ_ref[iz],wfR_ref[iz]))
#%%
def func_wfSRP(Xfpad,x_meters,m,t0):
    name_fp='WF_'+str(t0)+'s'+'_Mode'+str(m)
    x=x_meters/1000 # m to km
    x0=2500 #km
    
    ############
    # x0=2500 #km
    dist=1500 #km x0-X_SL
    dist1=x-X_SL # input x in meters
    if dist1<=distC: #behaviour exactly at boundary? 
        td_inc=-(dist1-dist)/cArrL[0]
        td_ref=-(distC-dist)/cArrL[0]-(distC-dist1)/cArrL[m]
        
        r0_inc=np.ones(Freqs.shape)
        r0_ref=-rcm[m]
        
        if m: #higher modes, before contrast, no incident wave
            r0_inc=np.zeros(Freqs.shape)
    else:
        td_inc=-(distC-dist)/cArrL[0]-(dist1-distC)/cArrR[m]
        td_ref=np.zeros(Freqs.shape)
        
        r0_inc=tcm[m]
        r0_ref=np.zeros(Freqs.shape) # Transmitted side, no reflected wave
    
    r0_inc1=Interp_Taper(r0_inc,Freqs,ff)
    td_inc1=Interp_Taper(td_inc,Freqs,ff)
    r0_ref1=Interp_Taper(r0_ref,Freqs,ff)
    td_ref1=Interp_Taper(td_ref,Freqs,ff)
    Xfpad1_inc=r0_inc1*Xfpad*np.exp(2j*pi*ff*td_inc1) 
    Xfpad1_ref=r0_ref1*Xfpad*np.exp(2j*pi*ff*td_ref1) 
    
    if dist1<=distC:
        Uz1=UzL1[m]
        Ur1=UrL1[m]
    else:
        Uz1=UzR1[m]
        Ur1=UrR1[m]
        
    # Xfpad1Z_inc=Uz1*Xfpad1_inc.reshape((1,-1));Xfpad1Z_inc[0]=dcZ
    # Xfpad1R_inc=-1j*Ur1*Xfpad1_inc.reshape((1,-1));Xfpad1R_inc[0]=dcR
    
    # Xfpad1Z_ref=Uz1*Xfpad1_ref.reshape((1,-1));Xfpad1Z_ref[0]=dcZ
    # Xfpad1R_ref=1j*Ur1*Xfpad1_ref.reshape((1,-1));Xfpad1R_ref[0]=dcR # reflection, propagation direction reversed
    
    Xfpad1Z_inc=Uz1*Xfpad1_inc.reshape((1,-1));
    Xfpad1R_inc=-1j*Ur1*Xfpad1_inc.reshape((1,-1));
    
    Xfpad1Z_ref=Uz1*Xfpad1_ref.reshape((1,-1));
    Xfpad1R_ref=1j*Ur1*Xfpad1_ref.reshape((1,-1)); # reflection, propagation direction reversed
    
    # DC value at 0 Hz??? Small enough to be set to 0
    Xfpad1Z_inc[:,0]=0;#dcZ
    Xfpad1R_inc[:,0]=0;#dcR
    Xfpad1Z_ref[:,0]=0;#dcZ
    Xfpad1R_ref[:,0]=0;#dcR

    ### Snapshot
    def Xfpad2wf(Xfpad1,ff,t0,dt,zz):
        N=len(ff)
        Xfpad1_t0=Xfpad1*exp(1j*2*pi*ff.reshape((1,-1))*t0)
        Xfpad1_t0[:,int(N/2)+1:]=np.conjugate(Xfpad1_t0[:,int(N/2-1):0:-1]) #conjugate symmetry
        Xfpad1_t0[:,int(N/2)]=Xfpad1_t0[:,int(N/2)].real #Midpoint, need to be real?
        wf=simps(Xfpad1_t0,dt*ff,axis=-1)
        wf=wf.real
        return wf

    wfZ_inc=Xfpad2wf(Xfpad1Z_inc,ff,t0,dt,zz)
    wfR_inc=Xfpad2wf(Xfpad1R_inc,ff,t0,dt,zz)
    wfZ_ref=Xfpad2wf(Xfpad1Z_ref,ff,t0,dt,zz)
    wfR_ref=Xfpad2wf(Xfpad1R_ref,ff,t0,dt,zz)
    
    with open (dir_fp+name_fp,'a+') as fp:
        for iz,z in enumerate(zz):
            fp.write('%f %f %.8E %.8E %.8E %.8E\n'%(z,x_meters,wfZ_inc[iz],wfR_inc[iz],wfZ_ref[iz],wfR_ref[iz]))
#%%
Para_Pool=[]
for m in range(10):
    for x_meters in tqdm(xx):
        Para_Pool.append((Xfpad,x_meters,m,t0))
#%%

pool=Pool(120)
pool.starmap(func_wfSRP,Para_Pool)
pool.close()
pool.join()


#%%
def func_wf(x,tmp):
    name_fp='WF_'+str(t0)+'s'
    x=x/1000 # m to km
    x0=2500 #km
    t_delay=(x-x0)/cArrL[0]

    tck=interpolate.splrep(Freqs,t_delay)  
    t_delay1=interpolate.splev(ff,tck);
    t_delay1=TaperCos1(t_delay1,ff,np.min(Freqs),np.max(Freqs))
    
    Xfpad1=Xfpad*np.exp(2j*pi*ff*(-t_delay1))
    Xfpad1Z=UzL1[0]*Xfpad1.reshape((1,-1))
    Xfpad1Z_t0=Xfpad1Z*exp(1j*2*pi*ff.reshape((1,-1))*t0)
    
    Xfpad1R=-1j*UrL1[0]*Xfpad1.reshape((1,-1)) #1j to account for phase shift
    Xfpad1R_t0=Xfpad1R*exp(1j*2*pi*ff.reshape((1,-1))*t0)

    
    # for iz in tqdm(range(len(zz))):
    for iz in range(len(zz)):
        Xfpad1Z_t0[iz,int(N/2)+1:]=np.conjugate(Xfpad1Z_t0[iz,int(N/2-1):0:-1]) #conjugate symmetry
        Xfpad1Z_t0[iz,int(N/2)]=Xfpad1Z_t0[iz,int(N/2)].real #Midpoint, need to be real?  
        Xfpad1R_t0[iz,int(N/2)+1:]=np.conjugate(Xfpad1R_t0[iz,int(N/2-1):0:-1]) #conjugate symmetry
        Xfpad1R_t0[iz,int(N/2)]=Xfpad1R_t0[iz,int(N/2)].real #Midpoint, need to be real? 
    wfZ=simps(Xfpad1Z_t0,dt*ff,axis=-1)
    wfR=simps(Xfpad1R_t0,dt*ff,axis=-1)      
    
    with open (dir_fp+name_fp,'a+') as fp:
        for iz,z in enumerate(zz):
            fp.write('%f %f %.8E %.8E\n'%(z,x,wfZ[iz],wfR[iz]))
    return
#%%
Para_Pool=[]
for x in tqdm(xx):
    Para_Pool.append((x-X_SL,0))

#%% !!!!!!!!!!!! Rename existing files !!!!!!!!!!!!! 
pool=Pool(100)
pool.starmap(func_wf,Para_Pool)
pool.close()
pool.join()
#%%
# results=np.genfromtxt(dir_fp+'WF_500s',names=['a','b','c','d'])
# results.sort(order=['a','b'])
# WFZ=results['c'].reshape((len(zz),len(xx)))
# WFR=results['d'].reshape((len(zz),len(xx)))

# results=np.genfromtxt(dir_fp+'WF_740s_Mode9',names=['a','b','c','d','e','f'])
# results.sort(order=['a','b'])
# WFZ_inc=results['c'].reshape((len(zz),len(xx)))
# WFR_inc=results['d'].reshape((len(zz),len(xx)))
# WFZ_ref=results['e'].reshape((len(zz),len(xx)))
# WFR_ref=results['f'].reshape((len(zz),len(xx)))

WF=np.zeros((len(zz),len(xx)))
for m in range(1,10):
    results=np.genfromtxt(dir_fp+'WF_740s_Mode'+str(m),names=['a','b','c','d','e','f'])
    results.sort(order=['a','b'])
    WFZ_inc=results['c'].reshape((len(zz),len(xx)))
    WFZ_ref=results['e'].reshape((len(zz),len(xx)))
    WF=WF+WFZ_inc+WFZ_ref
    
#%%


# WF=np.sqrt(WFZ**2+WFR**2)
# WF=WFR

# WF=np.sqrt(WFZ_ref**2+WFR_ref**2)
# WF=np.sqrt(WFZ_inc**2+WFR_inc**2)

# WF=WFR_inc+WFR_ref
# WF=WFZ_inc+WFZ_ref
# WF=np.sqrt((WFZ_inc+WFZ_ref)**2+(WFR_inc+WFR_ref)**2)
# WF=WFR_ref

LEVELS=np.linspace(np.min(WF),np.max(WF),100)
# LEVELS=np.linspace(-3e-5,3e-5,100)
# LEVELS=np.linspace(-2e-5,2e-5,100)
# LEVELS=np.linspace(0,3e-5,100)
# LEVELS=np.linspace(-6e-7,5e-7,100)
# LEVELS=np.linspace(-5e-7,4e-7,100)
# LEVELS=np.linspace(-3e-5,3e-5) #x
# LEVELS=np.linspace(-3e-5,3.5e-5) #z
LEVELS=np.linspace(-6e-7,7e-7)

plt.figure(figsize=(11,4))
plt.contourf(xx,zz,WF,levels=LEVELS,cmap='jet',extend='both')
# plt.gca().invert_yaxis()
# plt.ylim([6e5,0])
plt.ylim([94e4,0])
plt.colorbar()
#%%
# Xfpad1[int(N/2)+1:]=np.conjugate(Xfpad1[int(N/2-1):0:-1]) #conjugate symmetry
# Xfpad1[:int(N/2)]=Xfpad1[int(N/2)].real
# plt.figure()
# # plt.plot(tt1,xtpad,'k')
# plt.plot(tt1,ifft(Xfpad1),'r--')
#%%
# plt.figure()
# plt.plot(ff,np.abs(Xfpad1Z_t0[0,:]))
# #%%
# plt.figure()
# plt.imshow(WF,aspect='auto')
# # plt.colorbar()
# #%%
# plt.figure()
# plt.imshow(UZ1,aspect='auto')
# plt.colorbar()
# #%%
# plt.figure()
# plt.plot(ff,t_delay1)
#%%
# plt.figure()
# plt.plot(Freqs,Vph,'k')
# plt.plot(ff,Vph1,'r--')


# #%%
# iz=0
# plt.figure()
# plt.plot(Freqs,UZ[iz,:],'k')
# plt.plot(ff,U_interp,'r--')

# #%%
# iz=0
# plt.figure()
# plt.plot(Freqs,UZ[iz,:],'k.-')
# plt.plot(ff,UZ1[iz,:],'r.--')
# plt.xlim()



#%%
# plt.figure()
# plt.plot(tt1,xtpad,'k')
# plt.plot(tt1-20,xt1,'r--')

# #%%
# plt.figure()
# plt.plot(ff,np.abs(Xfpad)/np.max(np.abs(Xfpad)),'k')
# plt.plot(ff,np.abs(Xfpad1)/np.max(np.abs(Xfpad1)),'r--')
# #%%
# plt.figure()
# plt.plot(tt,xt,'k')
# # plt.plot(tt,xt_1/np.max(xt_1),'r--')
# plt.plot(tt,xt_3.real,'b:')
#%% Demonstration of waveform tapering
# freqf0=1/100.
# freqf1=1/40.
# tt, uz, st1_inc, st1_ref, st1, data_inc,st=ToSACModes(dir_rec,'S1951',2050,SOURCE,Mode)
# st=read(dir_SPECFEM+'S2051.EHZ.sac')

# st.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_inc.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_ref.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# uz.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
#%%
'''
vmax=5.
vmin=2.5
# k=np.max(stT[0].data)/np.max(st1[0].data)



# plt.figure()
# plt.plot(tt,st[0].data,'k',label='Original')
# plt.plot(tt,data_inc,'r--',label='Incident')
# plt.grid()
# plt.xlabel('time (s)')
# plt.xlim([-60,1000])
# plt.ylim([-6e-5,4e-5])
# plt.axvline(x=uz.stats.sac.dist/vmax,color='gray',linestyle='--')
# plt.axvline(x=uz.stats.sac.dist/vmin,color='gray',linestyle='--')
# plt.legend()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/5-1-1_Inc.pdf',transparent=False)


# st.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_inc.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_ref.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# uz.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)

# plt.figure()
# plt.plot(tt,st[0].data,'k',linewidth=2,label='Original')
# plt.plot(tt,st1_inc[0].data,'r--',label='Incident')
# plt.plot(tt,st1_ref[0].data,'b--',label='reflected')
# # plt.plot(tt,uz.data,'yellow',linestyle=(0,(5,10)),label='Reconstructed')
# plt.plot(tt,uz.data,'g:',linewidth=4,label='Reconstructed')
# plt.axvline(x=uz.stats.sac.dist/vmax,color='gray',linestyle='--')
# plt.axvline(x=uz.stats.sac.dist/vmin,color='gray',linestyle='--')

# plt.grid()
# plt.xlabel('time (s)')
# # plt.xlim([-60,1000])
# # plt.ylim([-6e-5,4e-5])
# plt.legend()
# plt.xlim([600,900])
# plt.ylim([-4e-6,4e-6])

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/5-1-2_Ref.pdf',transparent=False)


# st.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_inc.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# st1_ref.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)
# uz.filter('bandpass',freqmin=freqf0,freqmax=freqf1,corners=4,zerophase=True)

plt.figure()
plt.plot(tt,st[0].data,'k',linewidth=2,label='Original')
# plt.plot(tt,st1_inc[0].data,'r--',label='Incident')
# plt.plot(tt,st1_ref[0].data,'b--',label='reflected')
# plt.plot(tt,uz.data,'yellow',linestyle=(0,(5,10)),label='Reconstructed')
plt.plot(tt,uz.data,'g:',linewidth=4,label='Transmitted')
plt.axvline(x=uz.stats.sac.dist/vmax,color='gray',linestyle='--')
plt.axvline(x=uz.stats.sac.dist/vmin,color='gray',linestyle='--')

plt.grid()
plt.xlabel('time (s)')
plt.xlim([-60,1000])
plt.ylim([-6e-5,4e-5])
plt.legend()
# plt.xlim([600,900])
# plt.ylim([-4e-6,4e-6])

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/5-1-3_Trans.pdf',transparent=False)

#%%


# plt.figure()
# plt.plot(tt,st[0].data,'k',linewidth=2,label='Original')
# # plt.plot(tt,st1_inc[0].data,'r--',label='Incident')
# # plt.plot(tt,st1_ref[0].data,'b--',label='reflected')
# # plt.plot(tt,uz.data,'yellow',linestyle=(0,(5,10)),label='Reconstructed')
# plt.plot(tt,uz.data,'g:',linewidth=4,label='Fund Mode')
# plt.plot(tt,uzR.data,'y--',linewidth=2,label='Locked Mode')
# plt.axvline(x=uz.stats.sac.dist/vmax,color='gray',linestyle='--')
# plt.axvline(x=uz.stats.sac.dist/vmin,color='gray',linestyle='--')

# plt.grid()
# plt.xlabel('time (s)')
# # plt.xlim([-60,1000])
# # plt.ylim([-6e-5,4e-5])
# plt.legend()
# plt.xlim([600,900])
# plt.ylim([-2e-6,2e-6])
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/5-2-1_Ref.pdf',transparent=False)


plt.figure()
plt.plot(tt,st[0].data,'k',linewidth=2,label='Original')
# plt.plot(tt,st1_inc[0].data,'r--',label='Incident')
# plt.plot(tt,st1_ref[0].data,'b--',label='reflected')
# plt.plot(tt,uz.data,'yellow',linestyle=(0,(5,10)),label='Reconstructed')
plt.plot(tt,uz.data,'g:',linewidth=4,label='Fund Mode')
plt.plot(tt,uzT.data,'y--',linewidth=2,label='Locked Mode')
plt.axvline(x=uz.stats.sac.dist/vmax,color='gray',linestyle='--')
plt.axvline(x=uz.stats.sac.dist/vmin,color='gray',linestyle='--')

plt.grid()
plt.xlabel('time (s)')
plt.xlim([-60,1000])
plt.ylim([-4e-5,4e-5])
plt.legend()
# plt.xlim([600,900])
# plt.ylim([-2e-6,2e-6])
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/5-2-2_Trans.pdf',transparent=False)
'''
#%%

# Para_Pool=[]
# for idx in range(len(Stas)):
#     # dir_rec=dir_SK+'OUTPUT_SRP_SOURCE'+SOURCE+'_SPCETmpMode'+Mode+'_WithRef/'
#     sta=Stas[idx]
#     dist1=Dists[idx]
#     if SOURCE =='R':
#         dist1=4e3-dist1
#     Para_Pool.append((dir_rec,sta,dist1,SOURCE,Mode))

# #%%
# pool=Pool(100)
# results=pool.starmap(ToSACModes,Para_Pool)
# pool.close()
# pool.join()
# print('process finished')






