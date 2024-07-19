#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:21:47 2023

@author: u1318104
"""

import matplotlib.pyplot as plt
from obspy import read
import numpy as np
from tqdm import tqdm
import matplotlib.gridspec as gridspec
#%%
dir_tmp='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_LocWid/'
dir_SPECFEM_L=dir_tmp+'OUTPUT_SOURCEL_txt/SAC/'
dir_NM0=dir_tmp+'OUTPUT_SRP_SOURCEL_Lock_wavefieldGJI_Mode0_WithRef0_10/SAC/'
dir_NM10=dir_tmp+'OUTPUT_SRP_SOURCEL_Lock_wavefieldGJI_Mode10_WithRef0_10/SAC/'
dir_NM9=dir_tmp+'OUTPUT_SRP_SOURCEL_Lock_wavefieldGJI_Mode9_WithRef0_10/SAC/'


def read_wfs(dir_SPECFEM_L):
    # stas=np.arange(1501,2501+1,50,dtype=int)
    stas=np.arange(1001,3001+1,100,dtype=int)
    
    for i,sta in enumerate(stas):
        tmp=read(dir_SPECFEM_L+'S%d'%sta+'.EHZ.sac')
        # tmp[0].data=tmp[0].data*3.7e-5/np.max(tmp[0].data)
        if i==0:
            st=tmp
        else:
            st+=tmp
    
    st.filter('bandpass',freqmin=1/100,freqmax=1/40,corners=4,zerophase=True)
    t0=st[0].stats.sac.b
    dt=0.1#st[0].stats.sac.delta
    nt=st[0].stats.sac.npts
    tt=np.arange(t0,t0+nt*dt,dt)
    
    ccut=1.5e-6
    # ccut=1e-6
    # ccut=7e-7
    # ccut=4e-6
    thresh=1e-7
    st0=st.copy()
    stp=st.copy()
    stn=st.copy()
    for i in tqdm(range(len(stp))):
        tmp=st0[i].data
        tmp[tmp>ccut]=ccut
        tmp[tmp<-ccut]=-ccut
        st0[i].data=tmp
        
        tmp=st0[i].data.copy()
        tmp[tmp<thresh]=np.nan
        stp[i].data=tmp
        
        tmp=st0[i].data.copy()
        tmp[tmp>-thresh]=np.nan
        stn[i].data=tmp
    
    st1=st0.copy()
    stp1=stp.copy()
    stn1=stn.copy()
    
    # cnorm=np.max(np.abs(st[0].data))*3
    cnorm=ccut*2.5
    
    for i in range(len(st1)):
        st1[i].data=st1[i].data/cnorm
        stp1[i].data=stp1[i].data/cnorm
        stn1[i].data=stn1[i].data/cnorm
    
    return stas, thresh, tt, st1, stp1, stn1

stas, thresh, tt, st1, stp1, stn1 = read_wfs(dir_SPECFEM_L)
_, _, _, st1NM0,stp1NM0, stn1NM0 = read_wfs(dir_NM0)
_, _, _, st1NM10,stp1NM10, stn1NM10 = read_wfs(dir_NM10)
_, _, _, st1NM9,stp1NM9, stn1NM9 = read_wfs(dir_NM9)
#%%
plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'Times New Roman',
    "font.size": 11,
    "figure.autolayout": True}) #'axes.linewidth': 0.8

sep=5
XTICKS=np.arange(0,len(stas),sep)
XTICKSL=(stas[0:len(stas):sep]-1+1000)

fig=plt.figure(figsize=(8,4))
gs=gridspec.GridSpec(1,1)
ax1=fig.add_subplot(gs[0,0])
for i in range(len(st1)):
    ax1.fill_betweenx(tt,i+stp1[i].data,i+thresh,color='red')
    ax1.fill_betweenx(tt,i+stn1[i].data,i-thresh,color='blue')
    ax1.plot(i+st1[i].data,tt,color='k',linewidth=1)
# plt.ylim([tt[-1],tt[0]])
ax1.set_ylim([1200,0])
ax1.set_xticks(XTICKS,XTICKSL)
ax1.set_ylabel('Time (s)')
ax1.set_xlabel('Distance (km)')

TRANSPARENT=False
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig3_wf_SPECFEM.png',transparent=TRANSPARENT,dpi=300)
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig3_wf_SPECFEM.eps',transparent=TRANSPARENT)
#%%
sep=5
XTICKS=np.arange(0,len(stas),sep)
XTICKSL=(stas[0:len(stas):sep]-1+1000)

fig=plt.figure(figsize=(6,9))
gs=gridspec.GridSpec(3,1)

ax1=fig.add_subplot(gs[0,0])
for i in range(len(st1)):
    ax1.fill_betweenx(tt,i+stp1NM10[i].data,i+thresh,color='red')
    ax1.fill_betweenx(tt,i+stn1NM10[i].data,i-thresh,color='blue')
    ax1.plot(i+st1NM10[i].data,tt,color='k',linewidth=1)
# plt.ylim([tt[-1],tt[0]])
ax1.set_ylim([1200,0])
ax1.set_xticks(XTICKS,XTICKSL)
ax1.set_ylabel('Time (s)')

ax2=fig.add_subplot(gs[1,0])
for i in range(len(st1)):
    ax2.fill_betweenx(tt,i+stp1NM0[i].data,i+thresh,color='red')
    ax2.fill_betweenx(tt,i+stn1NM0[i].data,i-thresh,color='blue')
    ax2.plot(i+st1NM0[i].data,tt,color='k',linewidth=1)
# plt.ylim([tt[-1],tt[0]])
ax2.set_ylim([1200,0])
ax2.set_xticks(XTICKS,XTICKSL)
ax2.set_ylabel('Time (s)')

ax3=fig.add_subplot(gs[2,0])
for i in range(len(st1)):
    ax3.fill_betweenx(tt,i+stp1NM9[i].data,i+thresh,color='red')
    ax3.fill_betweenx(tt,i+stn1NM9[i].data,i-thresh,color='blue')
    ax3.plot(i+st1NM9[i].data,tt,color='k',linewidth=1)
# plt.ylim([tt[-1],tt[0]])
ax3.set_ylim([1200,0])
ax3.set_xticks(XTICKS,XTICKSL)
ax3.set_ylabel('Time (s)')

ax3.set_xlabel('Distance (km)')

TRANSPARENT=False
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig6_wf_NM.png',transparent=TRANSPARENT,dpi=300)
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig6_wf_NM.eps',transparent=TRANSPARENT)
#%%
plt.figure()
# plt.plot(tt,st[0].data,'--')
# plt.fill_between(tt,st[0].data,0,color='red')
for i in range(len(st1)):
    plt.fill_betweenx(tt,i+stp1[i].data,i+thresh,color='red')
    plt.fill_betweenx(tt,i+stn1[i].data,i-thresh,color='blue')
    plt.plot(i+st1[i].data,tt,color='k',linewidth=1)
    # break
plt.ylim([tt[-1],tt[0]])
# plt.xlim([-1,22])
# plt.swap_axes()
plt.ylim([1000,0])
#%%


plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'Times New Roman',
    "font.size": 11,
    "figure.autolayout": True}) #'axes.linewidth': 0.8

# fig,axs=plt.subplots(2,2,figsize=(8.5,6),gridspec_kw={'height_ratios':[1,2]})#gridspec_kw=
fig=plt.figure(figsize=(8.5,6))
gs=gridspec.GridSpec(2,10)

# gs=
# plt.figure(figsize=(11,4))
# plt.contourf(X1,Z1,M,levels=LEVELS,cmap='jet',extend='both')

ax1=fig.add_subplot(gs[0,:5])