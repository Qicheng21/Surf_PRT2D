#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 06:07:27 2022

@author: u1318104
"""
'''
1PSI calculation after synthesizing waveforms
'''
import numpy as np
from tqdm import tqdm
from numpy.linalg import inv
import matplotlib.pyplot as plt
from obspy import read
from multiprocessing import Pool
from scipy import interpolate
import sys
import matplotlib.gridspec as gridspec

sys.path.insert(1,'/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi')

dir_1PSI='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/RayleighWaveGJI/'
# dir_1PSIR='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP_SmoothBD/OUTPUT_R_0kkmBD_dXsta10km/SAC/'

Stas,Network = np.loadtxt(dir_1PSI+'/DATA/STATIONS',dtype=str,usecols=[0,1],unpack=True)
XXs,ZZs = np.loadtxt(dir_1PSI+'/DATA/STATIONS',dtype=float,usecols=[2,3],unpack=True) 
#%%
Comps=['Z'] #!!! Only vertical component
# Comps=['Z','R'] # Vertical and Radial Components

def get_valT(per,perA,valA):
    per1=perA[perA<per][-1]
    val1=valA[perA<per][-1]
    per2=perA[perA>per][0]
    val2=valA[perA>per][0]
    val=val1+(val2-val1)/(per2-per1)*(per-per1)
    return val

def judge_per(per,tt1,tt2,err):
    while tt1-tt2>per/2:
        tt1=tt1-per
    while tt2-tt1>per/2:
        tt1=tt1+per
    if abs(tt1-tt2)<err:
        return True
    else:
        return False

def func_data0_0(per,dir_SAC):
        
    dis_cri=3*per #3 km/s ref vel, half wavelength
    # dis_cri=1e-5
    snr_cri=5
    phase_cri=per/8
    
    
    # dir_per=dir_Coda+'HV/'+str(per)+'sec_RV_snr'+str(snr_cri)+'_dis'+str(dis_cri)+'_phase'+str(phase_cri)
    # os.system('mkdir '+dir_per)
    fALL=[]
    # sta1='033'
    
    # per=11
    fPhV_RR=[]
    fPhV_ZZ=[]
    
    k=0
    for ista in tqdm(range(len(Stas))): #['033']: #HV for sta1
        sta1=Stas[ista]
    
        tmp_st=read(dir_SAC+sta1+'.EHZ.sac')
        dist=tmp_st[0].stats.sac.dist
        if dist < dis_cri: #distance criteria, 3 wavelength?
            # print(str(per)+'s\t\t'+str(dis_cri)+'km\n'+sta1+'-'+sta2+' '+str(dist)+'km')
            continue
            
        snr2Comps=[];snr4CompsR=[]
        phv2Comps=[];phv4CompsR=[]
        amp2Comps=[];amp4CompsR=[]
        # try:
        for comp in Comps:
            compN=comp
          
            fp=dir_SAC+sta1+'.EH'+compN+'.sac'
            fpR=dir_SAC+sta1+'.EHR.sac'
 
            kA,cpA,apA,gvA,phvA,ampA=np.loadtxt(fp+'_1_DISP.0',usecols=[0,1,2,3,4,5],dtype=float,unpack=True)
            phv2Comps.append(get_valT(per,apA,phvA))
            amp2Comps.append(get_valT(per,apA,ampA))
        
        tt2Comps=dist/np.array(phv2Comps)
        
        # if sta1<sta2: #%% SNR for both components and Phase Shift Criteria
        # if snr2Comps[0]>snr_cri and snr2Comps[1]>snr_cri and judge_per(per, tt2Comps[0]+per/4, tt2Comps[1], phase_cri):
        # if judge_per(per, tt2Comps[0]+per/4, tt2Comps[1], phase_cri):
        if True:
            # tmpHV=amp4Comps[0]/amp4Comps[2]
            # HV.append(tmpHV)
            # fRV.append('%.6f %s %.6f %.6f'%(tmpHV,sta2,lon2,lat2))
            # pass
            fPhV_ZZ.append('%d %.6f %s %s %.6f %.8f'%(k,phv2Comps[0],sta1,'EHZ',dist,amp2Comps[0]))
            # fPhV_RR.append('%d %.6f %s %s %.6f'%(k,phv4CompsR[0],sta1,'EHR',dist))
    
            k+=1 
    
    np.savetxt(dir_SAC+'/data'+str(per)+'.0s_0.0.txt',fPhV_ZZ,fmt='%s')
    # np.savetxt(dir_SAC+'/data'+str(per)+'.0s_0.0.R.txt',fPhV_RR,fmt='%s')
    
    return tt2Comps


# suffix='Lock'
# func_data0_0(50,dir_1PSI+'/OUTPUT_SRP_SOURCEL_'+suffix+'/SAC/')
#%%
Para_Pool=[]

for per in tqdm([20,30,40,50,60,70,80,90]):
    Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL1000km/SAC/'))
    Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5000km/SAC/'))
    # for mode in [0,10]:
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL1000km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5000km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))
    # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL800km/SAC/'))
    # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5200km/SAC/'))
    # for mode in [0,10]:
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL800km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5200km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))


#%%
# pool=Pool(28)
# results=pool.starmap(func_data0_0,Para_Pool)
# pool.close()
# pool.join()
# print('Process Finished!')

#%%
eps=1e-5

# else:
X_SL=1000 #Left Souce Location
X_SR=5000 #Right Source Location

X_C=3000 #Center or structural location

X_RL=1000 #Leftmost receiver location
X_RR=5000 #Rightmost receiver location



def correct_2pi_progressive(XX,TT,per):
    # need first two points have correct vel, difference can be explained by 2pi ambiguity
    # mono increase XX
    tmp_TT=TT.copy()
    tmp_c=(XX[1]-XX[0])/(TT[1]-TT[0])
    for ii in range(len(XX)-2):
        tmp_c=(XX[ii+1]-XX[ii])/(TT[ii+1]-TT[ii])
        t0=TT[ii+1]+(XX[ii+2]-XX[ii+1])/tmp_c # predicted arrival time for next point
        n1=(t0-TT[ii+2])/per-1/2.0 #lower bound for nC
        n2=(t0-TT[ii+2])/per+1/2.0 #upper bound for nC
        if n1>0:
            nC=int(n2)
        elif n2<0:
            nC=int(n1)
        else:
            nC=0
        TT[ii+2:]=TT[ii+2:]+nC*per
    
    return XX,TT


class C_Psi1:
    def __init__(self,dir_SAC, dir_SACR,per,dXSta,dXShift,Comp='Z'):
        
        if Comp=='Z':
            PhV,Dist=np.loadtxt(dir_SAC+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4],unpack=True)
            PhVR,DistR=np.loadtxt(dir_SACR+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4],unpack=True)
        elif Comp=='R':
            PhV,Dist=np.loadtxt(dir_SAC+'/data'+str(per)+'.0s_0.0.R.txt',dtype=float,usecols=[1,4],unpack=True)
            PhVR,DistR=np.loadtxt(dir_SACR+'/data'+str(per)+'.0s_0.0.R.txt',dtype=float,usecols=[1,4],unpack=True)
        else:
            print('Wrong Component!')
            return
        TT_SL=Dist/PhV #Source from left
        
        TT_SR=DistR/PhVR #Source from right
        
        #%    
        
        XXL=Dist+X_SL
        XXR=X_SR-DistR
        tmp1=np.arange(dXSta/2+dXShift,X_RR-X_C+eps,dXSta)
        tmp2=np.arange(-dXSta/2+dXShift,X_RL-X_C-eps,-dXSta)
        XXSta=np.concatenate(((X_C+tmp2)[::-1],X_C+tmp1))
        
        TTStaL=np.zeros(len(XXSta))
        TTStaR=np.zeros(len(XXSta))
        
        
        for i in range(len(XXSta)):
            tmp=np.argmin(np.abs(XXL-XXSta[i]))
            TTStaL[i]=TT_SL[tmp]
            tmp=np.argmin(np.abs(XXR-XXSta[i]))
            TTStaR[i]=TT_SR[tmp]
        
        tmp_TTStaL=TTStaL.copy()
        tmp_TTStaR=TTStaR.copy()
        
        tmp_XX,TTStaL=correct_2pi_progressive(XXSta,TTStaL,per)
        tmp_XX,TTStaR=correct_2pi_progressive(XXSta,TTStaR,per)
        
        VphStaL=(XXSta[1:]-XXSta[:-1])/(TTStaL[1:]-TTStaL[:-1]) #np.abs
        VphStaR=(XXSta[1:]-XXSta[:-1])/(-TTStaR[1:]+TTStaR[:-1]) #np.abs
        PSI1=100*(VphStaL-VphStaR)/((VphStaL+VphStaR)/2)
        
        XXStaVph=(XXSta[1:]+XXSta[:-1])/2 # mid of stations, represent location of the phase velocity measurements
        
        self.XXL=XXL
        self.XXR=XXR
        self.TT_SL=TT_SL
        self.TT_SR=TT_SR
        
        self.XXSta=XXSta
        self.TTStaL=TTStaL
        self.TTStaR=TTStaR
        self.XXStaVph=XXStaVph
        self.VphStaL=VphStaL
        self.VphStaR=VphStaR
        self.PSI1=PSI1
        
        self.VphStaL_interp, self.VphStaL_interp_smooth, self.VphStaL_smooth, self.XXSta_interp, self.XXStaVph_interp=self.SmoothVph(self.TTStaL)
        self.VphStaR_interp, self.VphStaR_interp_smooth, self.VphStaR_smooth, self.XXSta_interp, self.XXStaVph_interp=self.SmoothVph(self.TTStaR)
        
        self.PSI1_interp=100*(self.VphStaL_interp-self.VphStaR_interp)/((self.VphStaL_interp+self.VphStaR_interp)/2)
        self.PSI1_interp_smooth=100*(self.VphStaL_interp_smooth-self.VphStaR_interp_smooth)/((self.VphStaL_interp_smooth+self.VphStaR_interp_smooth)/2)
        self.PSI1_smooth=100*(self.VphStaL_smooth-self.VphStaR_smooth)/((self.VphStaL_smooth+self.VphStaR_smooth)/2)
    
    def SmoothVph(self,TTSta):
        XStaL=2500
        XStaR=3500
        # dxGrid=22
        dxGrid=70/3
        idx=np.where(np.logical_and(self.XXSta<XStaR,self.XXSta>XStaL))
        # interpolate using b-spines
        XXSta_interp=np.arange(XStaL,XStaR+dxGrid,dxGrid) #interpolated x axis
        tck=interpolate.splrep(self.XXSta[idx],TTSta[idx])
        TTSta_interp=interpolate.splev(XXSta_interp,tck); 
        
        # VphSta=(self.XXSta[1:]-self.XXSta[:-1])/np.abs(TTSta[1:]-TTSta[:-1])
        VphSta=(self.XXSta[1:]-self.XXSta[:-1])/np.abs(TTSta[1:]-TTSta[:-1])
        SlowSta=np.abs(TTSta[1:]-TTSta[:-1])/(self.XXSta[1:]-self.XXSta[:-1])
        # XXStaVph=(self.XXSta[1:]+self.XXSta[:-1])/2
        VphSta_interp=(XXSta_interp[1:]-XXSta_interp[:-1])/np.abs(TTSta_interp[1:]-TTSta_interp[:-1])
        XXStaVph_interp=(XXSta_interp[1:]+XXSta_interp[:-1])/2
        #3-point average
        kernel_size=3
        # kernel_size=5
        kernel=np.ones(kernel_size)/kernel_size
        # kernel=np.array([1,0,0,1,0,0,1])
        # kernel=np.array([1])
        
        VphSta_interp_smooth=np.convolve(VphSta_interp,kernel,mode='same')
        VphSta_interp_smooth[0]=VphSta_interp[0];VphSta_interp_smooth[-1]=VphSta_interp[-1]
        VphSta_smooth=1/np.convolve(SlowSta,kernel,mode='same')#np.convolve(VphSta,kernel,mode='same')
        VphSta_smooth[0]=VphSta[0];VphSta_smooth[-1]=VphSta[-1]
        return VphSta_interp, VphSta_interp_smooth, VphSta_smooth, XXSta_interp, XXStaVph_interp

class C_Psi1H:
    def __init__(self,dir_SAC, dir_SACR,per,dXSta,dXShift,Comp='Z'):
        
        if Comp=='Z':
            PhV,Dist,Amp=np.loadtxt(dir_SAC+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4,5],unpack=True)
            
            PhVR,DistR,AmpR=np.loadtxt(dir_SACR+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4,5],unpack=True)
        elif Comp=='R':
            PhV,Dist=np.loadtxt(dir_SAC+'/data'+str(per)+'.0s_0.0.R.txt',dtype=float,usecols=[1,4],unpack=True)
            PhVR,DistR=np.loadtxt(dir_SACR+'/data'+str(per)+'.0s_0.0.R.txt',dtype=float,usecols=[1,4],unpack=True)
        else:
            print('Wrong Component!')
            return
        TT_SL=Dist/PhV #Source from left
        
        TT_SR=DistR/PhVR #Source from right
        
        #%    
        
        XXL=Dist+X_SL
        XXR=X_SR-DistR
        tmp1=np.arange(dXSta/2+dXShift,X_RR-X_C+eps,dXSta)
        tmp2=np.arange(-dXSta/2+dXShift,X_RL-X_C-eps,-dXSta)
        XXSta=np.concatenate(((X_C+tmp2)[::-1],X_C+tmp1))
        
        TTStaL=np.zeros(len(XXSta))
        TTStaR=np.zeros(len(XXSta))
        AAStaL=np.zeros(len(XXSta))
        AAStaR=np.zeros(len(XXSta))
        
        for i in range(len(XXSta)):
            tmp=np.argmin(np.abs(XXL-XXSta[i]))
            TTStaL[i]=TT_SL[tmp]
            AAStaL[i]=Amp[tmp]
            tmp=np.argmin(np.abs(XXR-XXSta[i]))
            TTStaR[i]=TT_SR[tmp]
            AAStaR[i]=AmpR[tmp]
        
        tmp_TTStaL=TTStaL.copy()
        tmp_TTStaR=TTStaR.copy()
        
        tmp_XX,TTStaL=correct_2pi_progressive(XXSta,TTStaL,per)
        tmp_XX,TTStaR=correct_2pi_progressive(XXSta,TTStaR,per)
        
        VphStaL=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaL[1:]-TTStaL[:-1]) #
        VphStaR=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaR[1:]-TTStaR[:-1]) #
        # PSI1=100*(VphStaL-VphStaR)/((VphStaL+VphStaR)/2)
        
        tmpA2=np.zeros(len(XXSta)-1);
        tmpA2[1:-1]=(AAStaL[3:]+AAStaL[:-3]-AAStaL[2:-1]-AAStaL[1:-2])/(2*(XXSta[2:-1]-XXSta[1:-2])**2)
        tmpA0=np.zeros(len(XXSta)-1);
        tmpA0[1:-1]=(AAStaL[1:-2]+AAStaL[2:-1])/2
        VphStaL=1/np.sqrt((np.abs(TTStaL[1:]-TTStaL[:-1])/(XXSta[1:]-XXSta[:-1]))**2-tmpA2/tmpA0/(2*np.pi/per)**2)
        
        tmpA2=np.zeros(len(XXSta)-1);
        tmpA2[1:-1]=(AAStaR[3:]+AAStaR[:-3]-AAStaR[2:-1]-AAStaR[1:-2])/(2*(XXSta[2:-1]-XXSta[1:-2])**2)
        tmpA0=np.zeros(len(XXSta)-1);
        tmpA0[1:-1]=(AAStaR[1:-2]+AAStaR[2:-1])/2
        VphStaR=1/np.sqrt((np.abs(TTStaR[1:]-TTStaR[:-1])/(XXSta[1:]-XXSta[:-1]))**2-tmpA2/tmpA0/(2*np.pi/per)**2)
        
        PSI1=100*(VphStaL-VphStaR)/((VphStaL+VphStaR)/2)
        
        XXStaVph=(XXSta[1:]+XXSta[:-1])/2 # mid of stations, represent location of the phase velocity measurements
        
        self.XXL=XXL
        self.XXR=XXR
        self.TT_SL=TT_SL
        self.TT_SR=TT_SR
        
        self.XXSta=XXSta
        self.TTStaL=TTStaL
        self.TTStaR=TTStaR
        self.AAStaL=AAStaL
        self.AAStaR=AAStaR
        self.XXStaVph=XXStaVph
        self.VphStaL=VphStaL
        self.VphStaR=VphStaR
        self.PSI1=PSI1
        
        self.VphStaL_interp, self.VphStaL_interp_smooth, self.VphStaL_smooth, self.XXSta_interp, self.XXStaVph_interp=self.SmoothVph(self.TTStaL,self.AAStaL)
        self.VphStaR_interp, self.VphStaR_interp_smooth, self.VphStaR_smooth, self.XXSta_interp, self.XXStaVph_interp=self.SmoothVph(self.TTStaR,self.AAStaR)
        
        self.PSI1_interp=100*(self.VphStaL_interp-self.VphStaR_interp)/((self.VphStaL_interp+self.VphStaR_interp)/2)
        self.PSI1_interp_smooth=100*(self.VphStaL_interp_smooth-self.VphStaR_interp_smooth)/((self.VphStaL_interp_smooth+self.VphStaR_interp_smooth)/2)
        self.PSI1_smooth=100*(self.VphStaL_smooth-self.VphStaR_smooth)/((self.VphStaL_smooth+self.VphStaR_smooth)/2)
    
    def SmoothVph(self,TTSta,AASta):
        XStaL=2500
        XStaR=3500
        # dxGrid=22
        dxGrid=70/3
        idx=np.where(np.logical_and(self.XXSta<XStaR,self.XXSta>XStaL))
        # interpolate using b-spines
        XXSta_interp=np.arange(XStaL,XStaR+dxGrid,dxGrid) #interpolated x axis
        tck=interpolate.splrep(self.XXSta[idx],TTSta[idx])
        TTSta_interp=interpolate.splev(XXSta_interp,tck); 
        
        # VphSta=(self.XXSta[1:]-self.XXSta[:-1])/np.abs(TTSta[1:]-TTSta[:-1])
        tmpA2=np.zeros(len(self.XXSta)-1);
        tmpA2[1:-1]=(AASta[3:]+AASta[:-3]-AASta[2:-1]-AASta[1:-2])/(2*(self.XXSta[2:-1]-self.XXSta[1:-2])**2)
        tmpA0=np.zeros(len(self.XXSta)-1);
        tmpA0[1:-1]=(AASta[1:-2]+AASta[2:-1])/2
        SlowSta=np.sqrt((np.abs(TTSta[1:]-TTSta[:-1])/(self.XXSta[1:]-self.XXSta[:-1]))**2-tmpA2/tmpA0/(2*np.pi/per)**2)
        
        VphSta=(self.XXSta[1:]-self.XXSta[:-1])/np.abs(TTSta[1:]-TTSta[:-1])
        # SlowSta=np.abs(TTSta[1:]-TTSta[:-1])/(self.XXSta[1:]-self.XXSta[:-1])
        # XXStaVph=(self.XXSta[1:]+self.XXSta[:-1])/2
        VphSta_interp=(XXSta_interp[1:]-XXSta_interp[:-1])/np.abs(TTSta_interp[1:]-TTSta_interp[:-1])
        XXStaVph_interp=(XXSta_interp[1:]+XXSta_interp[:-1])/2
        #3-point average
        kernel_size=3
        # kernel_size=5
        kernel=np.ones(kernel_size)/kernel_size
        # kernel=np.array([1,0,0,1,0,0,1])
        # kernel=np.array([1])
        
        VphSta_interp_smooth=np.convolve(VphSta_interp,kernel,mode='same')
        VphSta_interp_smooth[0]=VphSta_interp[0];VphSta_interp_smooth[-1]=VphSta_interp[-1]
        VphSta_smooth=1/np.convolve(SlowSta,kernel,mode='same')#np.convolve(VphSta,kernel,mode='same')
        VphSta_smooth[0]=VphSta[0];VphSta_smooth[-1]=VphSta[-1]
        return VphSta_interp, VphSta_interp_smooth, VphSta_smooth, XXSta_interp, XXStaVph_interp


#%%
plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'Times New Roman',
    "font.size": 11,
    "figure.autolayout": True}) #'axes.linewidth': 0.8
#%%
#40s-70km; 60s-100km; 80s-140km
per=60
dXSta=70
dXShift=0
dXSmooth=50

dirSL=dir_1PSI+'OUTPUT_SRP_SOURCEL%dkm/SAC/'%X_SL
dirSR=dir_1PSI+'OUTPUT_SRP_SOURCER%dkm/SAC/'%X_SR

dirM0L=dir_1PSI+'OUTPUT_SRP_SOURCEL%dkm_Lock_wavefieldGJI_Mode0_WithRef0_10/SAC/'%X_SL
dirM0R=dir_1PSI+'OUTPUT_SRP_SOURCER%dkm_Lock_wavefieldGJI_Mode0_WithRef0_10/SAC/'%X_SR

dirM10L=dir_1PSI+'OUTPUT_SRP_SOURCEL%dkm_Lock_wavefieldGJI_Mode10_WithRef0_10/SAC/'%X_SL
dirM10R=dir_1PSI+'OUTPUT_SRP_SOURCER%dkm_Lock_wavefieldGJI_Mode10_WithRef0_10/SAC/'%X_SR

CPS_10km=C_Psi1H(dirSL,dirSR,per,10,dXShift,'Z')
CPS_70km=C_Psi1H(dirSL,dirSR,per,dXSta,dXShift,'Z')

CP0_10km=C_Psi1(dirM0L,dirM0R,per,10,dXShift,'Z')
CP0_70km=C_Psi1(dirM0L,dirM0R,per,dXSta,dXShift,'Z')

CP10_10km=C_Psi1(dirM10L,dirM10R,per,10,dXShift,'Z')
CP10_70km=C_Psi1(dirM10L,dirM10R,per,dXSta,dXShift,'Z')

# CPS_10km=C_Psi1_wl(dirSL,dirSR,per,10,dXShift,dXSmooth,'Z')
# CPS_70km=C_Psi1_wl(dirSL,dirSR,per,dXSta,dXShift,dXSmooth,'Z')

# CP0_10km=C_Psi1_wl(dirM0L,dirM0R,per,10,dXShift,dXSmooth,'Z')
# CP0_70km=C_Psi1_wl(dirM0L,dirM0R,per,dXSta,dXShift,dXSmooth,'Z')

# CP10_10km=C_Psi1_wl(dirM10L,dirM10R,per,10,dXShift,dXSmooth,'Z')
# CP10_70km=C_Psi1_wl(dirM10L,dirM10R,per,dXSta,dXShift,dXSmooth,'Z')

# plt.figure()
# plt.plot(CPS_10km.XXL,CPS_10km.TT_SL,'.-')


#%
fig=plt.figure(figsize=(8.5,9))
gs=gridspec.GridSpec(3,1)

ax1=fig.add_subplot(gs[0,0])
ax1.plot(CPS_10km.XXStaVph,CPS_10km.VphStaL,'k--^',markersize=4,linewidth=2,label='dX=10km')
# ax1.plot(CPS_70km.XXStaVph,CPS_70km.VphStaL,'r--^',markersize=4,linewidth=2,label='L, SPECFEM Real')
ax1.plot(CPS_70km.XXStaVph,CPS_70km.VphStaL_smooth,'r--^',markersize=4,linewidth=2,label='dX=70km')

ax1.set_xlim([2500,3500])
ax1.set_ylim([3.5,4.1])
ax1.set_ylabel('Phase Velocity (km/s)')
ax1.legend()
ax1.grid()

ax2=fig.add_subplot(gs[1,0])
ax2.plot(CPS_10km.XXStaVph,CPS_10km.VphStaR,'k--^',markersize=4,label='dX=10km')
# ax2.plot(CPS_70km.XXStaVph,CPS_70km.VphStaR,'r--^',markersize=4,label='L, SPECFEM Real')
ax2.plot(CPS_70km.XXStaVph,CPS_70km.VphStaR_smooth,'r--^',markersize=4,linewidth=2,label='dX=70km')

ax2.set_xlim([2500,3500])
ax2.set_ylim([3.5,4.1])
ax2.set_ylabel('Phase Velocity (km/s)')
ax2.legend()
ax2.grid()

ax3=fig.add_subplot(gs[2,0])
ax3.plot(CPS_10km.XXStaVph,CPS_10km.PSI1,'kX--',label='dX=10km')
ax3.plot(CPS_70km.XXStaVph,CPS_70km.PSI1_smooth,'rX-',label='dX=70km')
# ax3.plot(CPS_70km.XXStaVph,CPS_70km.PSI1,'rX-',label='SPECFEM Smooth')

ax3.set_xlim([2500,3500])
ax3.set_ylim([-5,5])
ax3.set_ylabel('$A_1$ (\%)')
ax3.set_xlabel('Distance (km)')
ax3.legend()
ax3.grid()

TRANSPARENT=False
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/FigS1_1psi_SPECFEM2D.png',transparent=TRANSPARENT,dpi=300)
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/FigS1_1psi_SPECFEM2D.eps',transparent=TRANSPARENT)

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig4A_1psi_SPECFEM2D.png',transparent=TRANSPARENT,dpi=300)
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig4A_1psi_SPECFEM2D.eps',transparent=TRANSPARENT)

# plt.xlim([2500,3500])
# plt.ylim([-5,5])
#%%
fig=plt.figure(figsize=(8.5,9))
gs=gridspec.GridSpec(3,1)
    
ax1=fig.add_subplot(gs[0,0])
# ax1.plot(CP0_10km.XXStaVph,CP0_10km.VphStaL,'-',color=np.array([144,238,144])/255,markersize=4,linewidth=1,label='dx=10km, Mode 0')
# ax1.plot(CP10_10km.XXStaVph,CP10_10km.VphStaL,'-',color=np.array([144,144,238])/255,markersize=4,linewidth=1,label='dx=10km, Mode 0-9')

ax1.plot(CP0_70km.XXStaVph,CP0_70km.VphStaL_smooth,'g--^',markersize=4,linewidth=2,label='dx=70km, Mode 0')
ax1.plot(CP10_70km.XXStaVph,CP10_70km.VphStaL_smooth,'b--^',markersize=4,linewidth=2,label='dx=70km, Mode 0-9')

ax1.set_xlim([2500,3500])
# ax1.set_ylim([3.5,4.1])
# ax1.set_ylim([3.5,4.1]) #3.4,4.0

ax1.set_ylabel('Phase Velocity (km/s)')
ax1.legend()
ax1.grid()

ax2=fig.add_subplot(gs[1,0])
# ax2.plot(CPS_10km.XXStaVph,CP0_10km.VphStaR,'-',color=np.array([144,238,144])/255,markersize=4,linewidth=1,label='dx=10km, Mode 0')
# ax2.plot(CPS_10km.XXStaVph,CP10_10km.VphStaR,'-',color=np.array([144,144,238])/255,markersize=4,linewidth=1,label='dx=10km, Mode 0-9')
ax2.plot(CPS_70km.XXStaVph,CP0_70km.VphStaR_smooth,'g--^',markersize=4,linewidth=2,label='dx=70km, Mode 0')
ax2.plot(CPS_70km.XXStaVph,CP10_70km.VphStaR_smooth,'b--^',markersize=4,linewidth=2,label='dx=70km, Mode 0-9')
ax2.set_xlim([2500,3500])
# ax2.set_ylim([3.5,4.1])
# ax2.set_ylim([3.5,4.1])

ax2.set_ylabel('Phase Velocity (km/s)')
ax2.legend()
ax2.grid()

ax3=fig.add_subplot(gs[2,0])
ax3.plot(CPS_70km.XXStaVph,CPS_70km.PSI1_smooth,'rX-',label='SPECFEM2D')
ax3.plot(CP0_70km.XXStaVph,CP0_70km.PSI1_smooth,'gX-',label='Mode 0')
ax3.plot(CP10_70km.XXStaVph,CP10_70km.PSI1_smooth,'bX-',label='Mode 0-9')
# ax3.plot(CPS_70km.XXStaVph,CPS_70km.PSI1,'rX-',label='SPECFEM Smooth')
# ax3.plot(CP0_70km.XXStaVph,CP0_70km.PSI1,'gX-',label='SPECFEM Smooth')
# ax3.plot(CP10_70km.XXStaVph,CP10_70km.PSI1,'bX-',label='SPECFEM Smooth')
ax3.set_xlim([2500,3500])
# ax3.set_ylim([-3,3])
ax3.set_ylabel('$A_1$ (\%)')
ax3.set_xlabel('Distance (km)')

ax3.legend()
ax3.grid()

TRANSPARENT=False

