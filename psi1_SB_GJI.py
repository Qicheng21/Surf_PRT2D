#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 06:07:27 2022

@author: u1318104
"""



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
    # per=40
    # DSB=0
    # dir_SAC=dir_1PSI+'/OUTPUT_SRP_SK_SOURCEL_SBD'+str(DSB)+'km/SAC/'
        
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
            # if comp != 'ZZ':
            #     compN='ro_'+comp
            # fp=dir_Coda+tmpsta1+'/COR_'+tmpsta1+'_'+tmpsta2+'.SAC_'+compN
            # fpR=dir_Coda+tmpsta1+'/COR_'+tmpsta1+'_'+tmpsta2+'.SAC_'+compN+'_r'
            fp=dir_SAC+sta1+'.EH'+compN+'.sac'
            fpR=dir_SAC+sta1+'.EHR.sac'
            # perA,snrPreA,snrFolA=np.loadtxt(fp+'_snr_rms.txt',usecols=[0,1,2],dtype=float,unpack=True) #QCZ
            # snr2Comps.append(get_valT(per,perA,snrPreA)) #QCZ
            
            # try:
            #     perA,snrPreA,snrFolA=np.loadtxt(fpR+'_snr_rms.txt',usecols=[0,1,2],dtype=float,unpack=True)
            #     snr4CompsR.append(get_valT(per,perA,snrFolA))
            # except:
            #     pass
            
            kA,cpA,apA,gvA,phvA,ampA=np.loadtxt(fp+'_1_DISP.0',usecols=[0,1,2,3,4,5],dtype=float,unpack=True)
            phv2Comps.append(get_valT(per,apA,phvA))
            amp2Comps.append(get_valT(per,apA,ampA))
            # try:
            #     kA,cpA,apA,gvA,phvA,ampA=np.loadtxt(fpR+'_1_DISP.0',usecols=[0,1,2,3,4,5],dtype=float,unpack=True)
            #     phv4CompsR.append(get_valT(per,apA,phvA))
            #     amp4CompsR.append(get_valT(per,apA,ampA))   
            # except:
            #     pass
        # except:
        #     continue
        
        # if sum(np.array(phv2Comps)==0):
        #     print(sta1)
        #     continue
        # else:
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
# per=60
    # for tmpLayer in tqdm(range(60,61)):
    # for suffix in ['WithRef','WithOutRef']:
    # for suffix in ['Mw','Mw0km']:
    # for suffix in ['Mw_WithRef','Mw_WithOutRef']:
    # for suffix in ['Mw','Mw0km','Mw_WithRef','Mw_WithOutRef']:
    # for suffix in ['Mw0km','Mw1km','Mw5km','Mw10km','Mw20km']:
    # for suffix in ['X6000kmZ600km','X8000kmZ600km','X8000kmZ1000km']:
    # for suffix in ['Mode0_WithRef','ModeAll_WithRef']:    
    # for suffix in ['Mode0_WithRef']:
    # for suffix in ['Mw_WithRef']:
    # for suffix in ['CPSTmpMode0_WithRef','SPCETmpMode0_WithRef']:
    # for suffix in ['LockUB0_ModeAll_WithRef','LockUB0_Mode0_WithRef']:
    # for suffix in ['LockUB0TmpSPEC_ModeAll_WithRef']:
    # for suffix in ['Lock_ex12B_Vs6p9_Dep900_ModeAll_WithRef','Lock_ex13A_Vs6p9_Dep900_ModeAll_WithRef','Lock_ex13E_Vs6p9_Dep1300_ModeAll_WithRef']:
    # for suffix in ['Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef0_10','Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef10_10']:
    # for suffix in ['Lock_ex12B_Vs6p9_Dep900_ModeAll_WithRef_test0_10']:
    # for suffix in ['Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef_test0_10_RTcm']:
    # for suffix in ['Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef0_10']:
    # for suffix in ['Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef10_10_test3_nART']:
    # for suffix in ['Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef10_10_test4_nART_Taper']:
    # for suffix in['Lock_wavefield_ModeAll_WithRef0_0']:

    #     # dir_SAC=dir_1PSI+'/SK_Rho10perc/OUTPUT_SRP_SK_SOURCEL_Layer'+str(tmpLayer)+'_dVs0.00_dVp0.00_dRho10.00_PerturbSideR'+'/SAC/'
    #     # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SK_SOURCEL_SBD'+str(DSB)+'km/SAC/'))
    #     # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SK_SOURCEL_SBD0km_'+suffix+'/SAC/'))
    #     # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SK_SOURCER_SBD0km_'+suffix+'/SAC/'))
    #     # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL_'+suffix+'/SAC/'))
    #     # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER_'+suffix+'/SAC/'))
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SOURCEL_txt/SAC/'))
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SOURCER_txt/SAC/'))
    
    # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL1000km/SAC/'))
    # Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5000km/SAC/'))
    # for mode in [0,10]:
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCEL1000km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))
    #     Para_Pool.append((per,dir_1PSI+'/OUTPUT_SRP_SOURCER5000km_Lock_wavefieldGJI_Mode%d_WithRef0_10/SAC/'%mode))
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
# X_SL=500 #Left Souce Location
# X_SR=5500 #Right Source Location

# X_SL=800 #Left Souce Location
# X_SR=5200 #Right Source Location

X_C=3000 #Center or structural location

X_RL=1000 #Leftmost receiver location
X_RR=5000 #Rightmost receiver location

# X_SL=2000 #Left Souce Location
# X_SR=6000 #Right Source Location
# X_C=4000 #Center or structural location

# X_RL=2000 #Leftmost receiver location
# X_RR=6000 #Rightmost receiver location


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


def cal_1psi(per,dXSta,dXShift,DSB=0):
    
    dir_SAC=dir_1PSI+'OUTPUT_SRP_SK_SOURCEL_SBD'+str(DSB)+'km/SAC/'
    dir_SACR=dir_1PSI+'OUTPUT_SRP_SK_SOURCER_SBD'+str(DSB)+'km/SAC/'
    
    # dir_SAC=dir_1PSI+'/SK/OUTPUT_SRP_SK_SOURCEL_Layer'+str(0)+'_dVs0.00_dVp0.00_dRho0.00_PerturbSideR'+'/SAC/'
    
    PhV,Dist=np.loadtxt(dir_SAC+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4],unpack=True)
    TT_SL=Dist/PhV #Source from left
    PhVR,DistR=np.loadtxt(dir_SACR+'/data'+str(per)+'.0s_0.0.txt',dtype=float,usecols=[1,4],unpack=True)
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
    
    # tmp_XX,TTStaL=correct_2pi_progressive(XXSta,TTStaL,per)
    # tmp_XX,TTStaR=correct_2pi_progressive(XXSta,TTStaR,per)
    
    VphStaL=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaL[1:]-TTStaL[:-1])
    VphStaR=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaR[1:]-TTStaR[:-1])
    PSI1=(VphStaL-VphStaR)/((VphStaL+VphStaR)/2)
    
    XXStaVph=(XXSta[1:]+XXSta[:-1])/2 # mid of stations, represent location of the phase velocity measurements
    
    return XXStaVph, VphStaL, VphStaR, PSI1

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

class C_Psi1_wl:
    def __init__(self,dir_SAC, dir_SACR,per,dXSta,dXShift,dXSmooth,Comp='Z'):
        
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
        XXStaVph=(XXSta[1:]+XXSta[:-1])/2
        
        TTStaL=np.zeros(len(XXSta))
        TTStaR=np.zeros(len(XXSta))
        
        VphStaL_smooth=np.zeros(len(XXStaVph))
        VphStaR_smooth=np.zeros(len(XXStaVph))
        for i in range(len(XXStaVph)):
            tmp1=np.argmin(np.abs(XXL-XXStaVph[i]+dXSmooth/2))
            tmp2=np.argmin(np.abs(XXR-XXStaVph[i]-dXSmooth/2))
            VphStaL_smooth[i]=(XXL[tmp2]-XXL[tmp1])/np.abs(TT_SL[tmp2]-TT_SL[tmp1])
            VphStaR_smooth[i]=(XXR[tmp2]-XXR[tmp1])/np.abs(TT_SR[tmp2]-TT_SR[tmp1])
        
            tmp=np.argmin(np.abs(XXL-XXSta[i]))
            TTStaL[i]=TT_SL[tmp]
            # tmp=np.argmin(np.abs(XXR-XXSta[i]))
            TTStaR[i]=TT_SR[tmp]
        PSI1_smooth=100*(VphStaL_smooth-VphStaR_smooth)/((VphStaL_smooth+VphStaR_smooth)/2)
        
        VphStaL=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaL[1:]-TTStaL[:-1])
        VphStaR=(XXSta[1:]-XXSta[:-1])/np.abs(TTStaR[1:]-TTStaR[:-1])
        PSI1=100*(VphStaL-VphStaR)/((VphStaL+VphStaR)/2)
        
        # VphStaL
        # VphStaL_smooth
        # VphStaR
        # VphStaR_smooth
        # PSI1_smooth
        
        self.XXStaVph=XXStaVph
        self.VphStaL=VphStaL
        self.VphStaL_smooth=VphStaL_smooth
        self.VphStaR=VphStaR
        self.VphStaR_smooth=VphStaR_smooth
        self.PSI1_smooth=PSI1_smooth
        
        self.XXL=XXL
        self.XXR=XXR
        self.TT_SL=TT_SL
        self.TT_SR=TT_SR
        
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
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig7_1psi_NM.png',transparent=TRANSPARENT,dpi=300)
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig7_1psi_NM.eps',transparent=TRANSPARENT)

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig7A_1psi_NM.png',transparent=TRANSPARENT,dpi=300)
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig7A_1psi_NM.eps',transparent=TRANSPARENT)
#%%
from SWRT_Zeng_ex13_LockedDeep import EigenRead
dir_CPS='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/test_CPS330/'
dirM=dir_CPS+'ex12B_Vs6p9_Dep900'
effileL=dirM+'/SRDER_L.TXT'
effileR=dirM+'/SRDER_R.TXT'
freq=1/per
erL=EigenRead(freq,effileL)
erR=EigenRead(freq,effileR)

#%%
# id1=np.where(CPS.XXSta<X_C)[0][-1]
# id2=np.where(CPS.XXSta>X_C)[0][0]
# T_C=(CPS.TTStaL[id1]+CPS.TTStaL[id2])/2
# TTStaL_DM=np.zeros(CPS.TTStaL.shape)
# for i,x in enumerate(CPS.XXSta):
#     if x<X_C:
#         TTStaL_DM[i]=CPS.TTStaL[i]-(x-X_C)/erL.Vph[0]-T_C
#     else:
#         TTStaL_DM[i]=CPS.TTStaL[i]-(x-X_C)/erR.Vph[0]-T_C

# plt.close('all')
plt.figure()
# plt.plot(CPS.XXSta,CPS.TTStaL,'r--^',markersize=4,label='Left Source')
# plt.plot(CPS.XXSta,CPS.TTStaR,'b--.',label='Right Source')
plt.plot(CP2.XXSta,CP2.TTStaL,'r--^',markersize=4,label='Left Source')
plt.plot(CP2.XXSta,CP2.TTStaR,'b--.',label='Right Source')
plt.legend()
plt.xlim([2500,3500])
plt.ylim([300,800])
# plt.axvline(x=X_C,color='gray')
plt.grid()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/SPECFEM_T_dX70km.pdf')

# plt.figure()
# plt.plot(CPS.XXSta, TTStaL_DM,'r-o')
# plt.xlim([2500,3500])
# plt.ylim([-1,1])
#%%
# dirRec=dir_1PSI+'/OUTPUT_SRP_SOURCEL_X6000kmZ600km_v2_44CPU/SAC/'
dirRec=dir_1PSI+'/OUTPUT_SRP_SOURCEL_Lock_ex14A_Vs6p9_Dep900_M10_ModeAll_WithRef0_10/SAC/'
st=read(dirRec+'S0001.EHZ.sac')
tb=st[0].stats.sac.b
dt=st[0].stats.sac.delta
nt=st[0].stats.sac.npts
tt=np.arange(tb,tb+nt*dt,dt)
#%%
st.clear()
plt.figure()
k=1e2*3
Dists=np.arange(1000,5001,100)
# for i, dis in tqdm(enumerate(CPS.XXSta)):
for i, dis in tqdm(enumerate(Dists)):

    staInt=str(int(dis+1-X_RL))
    staStr=(4-len(staInt))*'0'+str(staInt)
    nm=str(dis+1)
    print(staStr)
    st=read(dirRec+'S'+staStr+'.EHZ.sac')
    # print(st)
    plt.plot(tt,st[0].data/np.max(st[0].data)*k+dis,'k')

# plt.ylim([2500,3500])
# plt.xlim([300,900])'
plt.ylim([5500,500])
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/SPECFEM_Record.pdf')


#%%
plt.figure()

# plt.plot(np.linspace(1e3,3e3,10),erL.Vph[0]*np.ones(10),'k',label='1D Vph')
# plt.plot(np.linspace(3e3,5e3,10),erR.Vph[0]*np.ones(10),'k')

# plt.plot(CPS.XXStaVph,CPS.VphStaL,'r--^',markersize=4,label='L, SPECFEM Real')
# plt.plot(CPS.XXStaVph,CPS.VphStaR,'b--.',label='R, SPECFEM Real')


plt.plot(CPS.XXStaVph,CPS.VphStaL,'k--^',markersize=4,label='L, SPECFEM Real')
plt.plot(CPS.XXStaVph,CPS.VphStaR,'k--.',label='R, SPECFEM Real')

# plt.plot(CPR.XXStaVph,CPR.VphStaL,'k--^',label='L, SPECFEM Predict')
# plt.plot(CPR.XXStaVph,CPR.VphStaR,'k--.',label='R, SPECFEM Predict')
# plt.plot(CPR.XXStaVph,CPR.VphStaL,'r-^',label='L, SPECFEM Predict')
# plt.plot(CPR.XXStaVph,CPR.VphStaR,'b-.',label='R, SPECFEM Predict')

# plt.plot(CP0.XXStaVph,CP0.VphStaL,'r-^',label='L, Normal Modes')
# plt.plot(CP0.XXStaVph,CP0.VphStaR,'b-.',label='R, Normal Modes')

# plt.plot(CP0.XXStaVph,CP0.VphStaL,'r--^',markersize=5,label='L, Locked Mode 0')
# plt.plot(CP0.XXStaVph,CP0.VphStaR,'b-.',label='R, Locked Mode 0')

# plt.plot(CP0.XXStaVph,CP1.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP0.XXStaVph,CP1.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CP2.XXStaVph,CP2.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP2.XXStaVph,CP2.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CP4.XXStaVph,CP4.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP4.XXStaVph,CP4.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CP5.XXStaVph,CP5.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP5.XXStaVph,CP5.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CP2.XXStaVph,CP8.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP2.XXStaVph,CP8.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CP2.XXStaVph,CP10.VphStaL,'r-^',label='L, Locked Mode All')
# plt.plot(CP2.XXStaVph,CP10.VphStaR,'b-.',label='R, Locked Mode All')

plt.plot(CP2.XXStaVph,CP14.VphStaL,'r-^',label='L, Locked Mode All')
plt.plot(CP2.XXStaVph,CP14.VphStaR,'b-.',label='R, Locked Mode All')

# plt.plot(CPA.XXStaVph,CPA.VphStaL,'r-^')
# plt.plot(CPA.XXStaVph,CPA.VphStaR,'b-.')



# plt.ylim([])
plt.legend()
plt.grid()
plt.xlim([2500,3500])
plt.ylim([3.5,4.1])
# plt.ylim([3.2,4.3])

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/VphCP8_dxSta70km')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/2-2_Vph_dX70km.pdf')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/11112022_LockedUB0/Vph_'+str(per)+'s_12B1')


#%%
# plt.close('all')
# plt.figure(figsize=(10,8))
plt.figure()
# plt.plot(CPS.XXStaVph,CPS.PSI1,'ko-',label='SPECFEM Ori')
# plt.plot(CPS.XXStaVph_interp,CPS.PSI1_interp,'r^-',label='SPECFEM Interp')
# plt.plot(CP.XXStaVph_interp,CP.PSI1_interp_smooth,'bs-',label='SPECFEM Interp_Smooth')
plt.plot(CPS.XXStaVph,CPS.PSI1_smooth,'kX-',label='SPECFEM Smooth')
# plt.plot(CPS.XXStaVph,CPS.PSI1,'kX-',label='SPECFEM Smooth')

# plt.plot(CPRef1.XXStaVph,CPRef1.PSI1,'k--',label='ModeSum Ori')
# plt.plot(CPRef1.XXStaVph_interp,CPRef1.PSI1_interp,'r--',label='ModeSum Interp')
# plt.plot(CPRef1.XXStaVph_interp,CPRef1.PSI1_interp_smooth,'b--',label='ModeSum Interp_Smooth')
# plt.plot(CP1.XXStaVph,CP1.PSI1_smooth,'g--',label='ModeSum Smooth')
# plt.plot(CP0.XXStaVph,CP0.PSI1_smooth,'g--',label='ModeSum Smooth')
# plt.plot(CP1.XXStaVph,CP1.PSI1_smooth,'g--',label='ModeSum Smooth')
plt.plot(CP0.XXStaVph,CP0.PSI1_smooth,'bX:',label='ModeSum Fund Smooth')
# plt.plot(CP1.XXStaVph,CP1.PSI1_smooth,'rX--',label='ModeSum LockM10 Smooth')
# plt.plot(CP1.XXStaVph,CP1.PSI1,'g--',label='ModeSum LockM10 Smooth')

# plt.plot(CP2.XXStaVph,CP2.PSI1_smooth,'rX--',label='ModeSum LockM10 Smooth')
plt.plot(CP8.XXStaVph,CP14.PSI1_smooth,'rX--',label='ModeSum LockM10 Smooth')


plt.title('%d s'%per)
plt.legend()
plt.xlim([2500,3500])
plt.ylim([-5,5])
plt.grid()
plt.xlabel('Station Location (km)')
plt.ylabel('1-psi (%)')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/2-3_1psi_dX70km.pdf')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/Psi1CP9_dxSta70km')

# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/11112022_LockedUB0/PSI1_'+str(per)+'s_13A')
#%%
'''
# XStaL=2500
# XStaR=3500
# dxGrid=22
# idxL=np.where(np.logical_and(CP_Z.XXSta<XStaR,CP_Z.XXSta>XStaL))
# # interpolate using b-spines
# tck_TTStaL=interpolate.splrep(CP_Z.XXSta[idxL],CP_Z.TTStaL[idxL])
# XXSta_interp=np.arange(XStaL,XStaR+dxGrid,dxGrid)
# TTStaL_interp=interpolate.splev(XXSta_interp,tck_TTStaL); # works well inside Freqs, need filtering afterwards

# VphStaL=(CP_Z.XXSta[1:]-CP_Z.XXSta[:-1])/np.abs(CP_Z.TTStaL[1:]-CP_Z.TTStaL[:-1])
# VphStaL_interp=(XXSta_interp[1:]-XXSta_interp[:-1])/np.abs(TTStaL_interp[1:]-TTStaL_interp[:-1])
# XXStaVph=(CP_Z.XXSta[1:]+CP_Z.XXSta[:-1])/2
# XXStaVph_interp=(XXSta_interp[1:]+XXSta_interp[:-1])/2
# #3-point average
# kernel_size=3
# kernel=np.ones(kernel_size)/kernel_size
# VphStaL_interp_smooth=np.convolve(VphStaL_interp,kernel,mode='same')

# plt.close('all')
# plt.figure()
# plt.plot(CP.XXSta,CP.TTStaL,'ko-')
# plt.plot(CP.XXSta_interp,CP.TTStaL_interp,'r^--')

# plt.figure(figsize=(8,6))
plt.figure()
plt.plot(CP0.XXStaVph,CP0.VphStaL,'ko-',label='SPECFEM Ori')

# plt.plot(CPAll.XXStaVph,CPAll.VphStaL,'ko-',label='SPECFEM Ori')
# plt.plot(CP.XXStaVph_interp,CP.VphStaL_interp,'r^-',label='SPECFEM Interp')
# plt.plot(CP.XXStaVph_interp,CP.VphStaL_interp_smooth,'bs-',label='SPECFEM Interp_Smooth')
# plt.plot(CP.XXStaVph,CP.VphStaL_smooth,'gX-',label='SPECFEM Smooth')

# plt.plot(CPRef1.XXStaVph,CPRef1.VphStaL,'k--',label='ModeSum Ori')
# plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaL_interp,'r--',label='ModeSum Interp')
# plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaL_interp_smooth,'b--',label='ModeSum Interp_Smooth')
# plt.plot(CPRef1.XXStaVph,CPRef1.VphStaL_smooth,'g--',label='ModeSum Smooth')

plt.xlim([2500,3500])
plt.ylim([3.5,4.3])
# plt.ylim([3.7,4.1])
plt.grid()
plt.legend()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/08262022_VphSmooth/VphL_'+str(per)+'s',transparent=True)

# # plt.figure(figsize=(8,6))
# plt.figure()
# plt.plot(CP0.XXStaVph,CP0.VphStaR,'ko-',label='SPECFEM Ori')
# # plt.plot(CP.XXStaVph_interp,CP.VphStaR_interp,'r^-',label='SPECFEM Interp')
# # plt.plot(CP.XXStaVph_interp,CP.VphStaR_interp_smooth,'bs-',label='SPECFEM Interp_Smooth')
# # plt.plot(CP.XXStaVph,CP.VphStaR_smooth,'gX-',label='SPECFEM Smooth')

# # plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR,'k--',label='ModeSum Ori')
# # plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaR_interp,'r--',label='ModeSum Interp')
# # plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaR_interp_smooth,'b--',label='ModeSum Interp_Smooth')
# # plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR_smooth,'g--',label='ModeSum Smooth')

# plt.xlim([2500,3500])
# plt.ylim([3.5,4.3])
# # plt.ylim([3.7,4.1])
# plt.grid()
# plt.legend()
# # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/08262022_VphSmooth/VphR_'+str(per)+'s',transparent=True)

# plt.figure()
# plt.plot(CP.XXStaVph,CP.VphStaR,'ko-')
# # plt.plot(CP.XXStaVph_interp,CP.VphStaR_interp,'r^-')
# # plt.plot(CP.XXStaVph_interp,CP.VphStaR_interp_smooth,'bs-')
# # plt.plot(CP.XXStaVph,CP.VphStaR_smooth,'gX-')

# plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR,'k--')
# # plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaR_interp,'r--')
# # plt.plot(CPRef1.XXStaVph_interp,CPRef1.VphStaR_interp_smooth,'b--')
# # plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR_smooth,'g--')

# plt.xlim([2500,3500])
# # plt.ylim([3.5,4.3])
# plt.ylim([3.7,4.1])
# plt.grid()
#%%
plt.figure()
# plt.plot(CP.XXStaVph,CP.VphStaL,'r--^')
# plt.plot(CP.XXStaVph,CP.VphStaR,'b--.')
plt.plot(CP_Z.XXStaVph,CP_Z.VphStaL,'k-^',label='SL SPECFEM')
plt.plot(CP_Z.XXStaVph,CP_Z.VphStaR,'ko-',label='SR SPECFEM')

# plt.plot(CPRef1.XXStaVph,CPRef1.VphStaL,'r--^',label='SL Synthetic')
# plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR,'bo--',label='SR Synthetic')

# plt.plot(CP_R.XXStaVph,CP_R.VphStaL,'r--^',label='SL Synthetic')
# plt.plot(CP_R.XXStaVph,CP_R.VphStaR,'bo--',label='SR Synthetic')

# plt.plot(CPRef0.XXStaVph,CPRef0.VphStaL,'r:^')
# plt.plot(CPRef0.XXStaVph,CPRef0.VphStaR,'b:.')

plt.xlim([2500,3500])
# plt.xlim([3500,4500])

plt.ylim([3.7,4.1])
plt.ylim([3.5,4.3])

# plt.ylim([3.6,4.0])
# plt.ylim([3.5,3.9])
plt.grid()
plt.legend()

#%%
per=80
dXSta=30
dXShift=0
# dirL=dir_1PSI+'/OUTPUT_SRP_SK_SOURCEL_SBD0km'
# dirR=dir_1PSI+'/OUTPUT_SRP_SK_SOURCER_SBD0km'
dirL=dir_1PSI+'/OUTPUT_SRP_SOURCEL_Mw'
dirR=dir_1PSI+'/OUTPUT_SRP_SOURCER_Mw'
CP=C_Psi1(dirL+'/SAC/',dirR+'/SAC/',per,dXSta,dXShift)
CPRef1=C_Psi1(dirL+'_WithRef/SAC/',dirR+'_WithRef/SAC/',per,dXSta,dXShift)
CPRef0=C_Psi1(dirL+'_WithOutRef/SAC',dirR+'_WithOutRef/SAC/',per,dXSta,dXShift)
#%%
# plt.close('all')
plt.figure()
# plt.plot(CP.XXStaVph,CP.VphStaL,'r--^')
# plt.plot(CP.XXStaVph,CP.VphStaR,'b--.')
plt.plot(CP.XXStaVph,CP.VphStaL,'k-^',label='SL SPECFEM')
plt.plot(CP.XXStaVph,CP.VphStaR,'ko-',label='SR SPECFEM')

plt.plot(CPRef1.XXStaVph,CPRef1.VphStaL,'r--^',label='SL Synthetic')
plt.plot(CPRef1.XXStaVph,CPRef1.VphStaR,'bo--',label='SR Synthetic')

# plt.plot(CPRef0.XXStaVph,CPRef0.VphStaL,'r:^')
# plt.plot(CPRef0.XXStaVph,CPRef0.VphStaR,'b:.')

plt.xlim([2500,3500])
plt.ylim([3.5,4.3])

# plt.ylim([3.7,4.1])
# plt.ylim([3.6,4.0])
# plt.ylim([3.5,3.9])
plt.grid()
plt.legend()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/08112022_SRP_SyntheticWaveform/Compare_'+str(per)+'s_'+str(dXSta)+'km')

#%%
ArrCP=[]
ArrCPRef1=[]
ArrCPRef0=[]
Periods=[40,60,80]
for per in tqdm(Periods):
    dXSta=70
    dXShift=0
    dirL=dir_1PSI+'/OUTPUT_SRP_SK_SOURCEL_SBD0km'
    dirR=dir_1PSI+'/OUTPUT_SRP_SK_SOURCER_SBD0km'
    CP=C_Psi1(dirL+'/SAC/',dirR+'/SAC/',per,dXSta,dXShift)
    CPRef1=C_Psi1(dirL+'_WithRef/SAC/',dirR+'_WithRef/SAC/',per,dXSta,dXShift)
    CPRef0=C_Psi1(dirL+'_WithOutRef/SAC',dirR+'_WithOutRef/SAC/',per,dXSta,dXShift)
    
    ArrCP.append(CP)
    ArrCPRef1.append(CPRef1)
    ArrCPRef0.append(CPRef0)
#%%
# plt.close('all')
plt.figure()
ColorsCP=['r','g','b']
for i,iColor in enumerate(ColorsCP):
    # plt.plot(ArrCP[i].XXL,ArrCPRef1[i].TT_SL-ArrCPRef0[i].TT_SL,color=iColor,linestyle='-.',label=str(Periods[i])+' s')
    plt.plot(ArrCP[i].XXR,ArrCPRef1[i].TT_SR-ArrCPRef0[i].TT_SR,color=iColor,linestyle='-.',label=str(Periods[i])+' s')
    # plt.plot(ArrCP[i].XXL,ArrCP[i].TT_SL-ArrCPRef1[i].TT_SL,color=iColor,linestyle='-.',label=str(Periods[i])+' s')
    # plt.plot(ArrCP[i].XXR,ArrCP[i].TT_SR-ArrCPRef1[i].TT_SR,color=iColor,linestyle='-.',label=str(Periods[i])+' s')

plt.legend()
plt.xlim([2500,3500])
plt.xlabel('Station Location (km)')
plt.ylabel('Time Difference (s)')
plt.grid()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/08112022_SRP_SyntheticWaveform/\
# TimeDifference_SpecfemRef1_SL')
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/08112022_SRP_SyntheticWaveform/\
# TimeDifference_Ref1Ref0_SR')

# plt.plot(CP.XXL,CP.TT_SL,'r^--')
# plt.plot(CP.XXR,CP.TT_SR,'bo--')
# plt.plot(CP.XXL,CP.TT_SL,'k')
# plt.plot(CP.XXR,CP.TT_SR,'k')

# plt.plot(CPRef1.XXL,CPRef1.TT_SL,'r--')
# plt.plot(CPRef1.XXR,CPRef1.TT_SR,'r--')

# plt.plot(CPRef0.XXL,CPRef0.TT_SL,'b:')
# plt.plot(CPRef0.XXR,CPRef0.TT_SR,'b:')

'''
'''
#%%
DSBs=np.arange(11)*20
dXShifts=np.arange(-34,34+eps,1)
Pers=np.arange(40,100,10)
dXSta=70 #fixed station spacing, km

CriDisPsi1=500 #distance range criteria for 1psi measurements, km, one side
nCDP=int(CriDisPsi1/dXSta) #number of 1psi measurements allowed, one side
NCDP=2*nCDP+1 # 2 sides including middle point
PSI1s=np.zeros([len(DSBs),len(dXShifts),len(Pers),NCDP]) #15 psi measurements, range 3000 +- dXSta*7
for i in range(len(DSBs)):
    for j in tqdm(range(len(dXShifts))):
        for k in range(len(Pers)):
            XXStaVph, VphStaL, VphStaR, PSI1=cal_1psi(Pers[k],dXSta,dXShifts[j],DSBs[i])
            idxM=np.argmin(np.abs(XXStaVph-X_C))
            PSI1s[i,j,k,:]=PSI1[idxM-nCDP:idxM+nCDP+1]
            
#%%
i=10
j=34
LStyles=['r','g','b','c','m','y']
plt.rcParams.update({'font.size': 12})
n=1

for i in [0]:#[0,5,10]:
    for j in [0]:#[0,34,-1]:
        # plt.subplot(3,3,n)
# n+=1
        plt.figure(figsize=(6.4,4.8))
        for k in range(len(Pers)):
            XXStaVph, VphStaL, VphStaR, PSI1=cal_1psi(Pers[k],dXSta,dXShifts[j],DSBs[i])
            idxM=np.argmin(np.abs(XXStaVph-X_C))
            plt.plot(XXStaVph,PSI1*100,LStyles[k]+'.-',label=str(Pers[k])+' s')
        plt.legend()
        plt.ylim([-4,4])
        plt.xlim([2500,3500])
        plt.grid()
        plt.xlabel('Station Location (km)')
        plt.ylabel('% 1psi')
        plt.title('BD %skm, dXShift %skm'%(DSBs[i],dXShifts[j]))
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/06232022_SBDInv_CPS/BD%skmdXShift%skm.eps'%(DSBs[i],dXShifts[j]),transparent=True)
#%%
from def_colormap import colormap
cmnew=colormap()

iDSBs=3
jdXShifts=30
PSI1=PSI1s[iDSBs,jdXShifts,:,:]
# UnPSI1=np.random.normal(0,10,PSI1.shape)
# PSI1=PSI1*(1+UnPSI1/100)

Errs=np.zeros([len(DSBs),len(dXShifts)])
for i in range(len(DSBs)):
    for j in range(len(dXShifts)):
        Errs[i,j]=np.sqrt(np.sum((PSI1s[i,j,:,:]-PSI1)**2)/(len(Pers)*NCDP))

XX,YY=np.meshgrid(dXShifts,DSBs)
levels=np.linspace(0,np.max(Errs),100)
levelsC=np.linspace(0,np.max(Errs),10)

v=np.arange(0,1.1,0.1)
plt.figure()
plt.contourf(XX,YY,Errs*100,levels=levels*100,cmap=cmnew)
plt.colorbar(ticks=v,format='%.2f')

CS=plt.contour(XX,YY,Errs*100,v,colors='k',linestyles='--',linewidth=1,zorder=8)
# plt.clabel(CS,v, inline=True,fmt='%.1f')

plt.plot(dXShifts[jdXShifts],DSBs[iDSBs],'k*',markersize=10,label='True')
plt.plot(dXShifts[np.where(Errs==np.min(Errs))[1][0]],DSBs[np.where(Errs==np.min(Errs))[0][0]],'m^',markersize=4,label='MinMis')
plt.xlabel('dXShifts (km)')
plt.ylabel('DSBs (km)')
plt.legend()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/06232022_SBDInv_CPS/Err_BD%skmdXShift%skm.png'%(DSBs[i],dXShifts[j]),transparent=True)

# XXStaVph, VphStaL, VphStaR, PSI1=cal_1psi(60,70,0,0)
# idxM=np.argmin(np.abs(XXStaVph-X_C))
# Arr_PSI1_dXShift[idXShift]=PSI1[idxM]

#%%

iDSBs=3
jdXShifts=30

Nrand=1000
NdXShifts=np.zeros(Nrand)
NDSBs=np.zeros(Nrand)

for n in tqdm(np.arange(Nrand)):
    PSI1=PSI1s[iDSBs,jdXShifts,:,:]
    UnPSI1=np.random.normal(0,50,PSI1.shape)
    PSI1=PSI1*(1+UnPSI1/100)
    
    Errs=np.zeros([len(DSBs),len(dXShifts)])
    for i in range(len(DSBs)):
        for j in range(len(dXShifts)):
            Errs[i,j]=np.sqrt(np.sum((PSI1s[i,j,:,:]-PSI1)**2)/(len(Pers)*NCDP))
            
            
    NdXShifts[n]=dXShifts[np.where(Errs==np.min(Errs))[1][0]]
    NDSBs[n]=DSBs[np.where(Errs==np.min(Errs))[0][0]]
#%
print(np.mean(NdXShifts), np.sqrt(np.mean((NdXShifts-dXShifts[jdXShifts])**2)))
print(np.mean(NDSBs), np.sqrt(np.mean((NDSBs-DSBs[iDSBs])**2)))

'''


