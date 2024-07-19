#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 09:38:31 2022

@author: u1318104
"""

#%% Create input velocity model for the two sides
import os
import numpy as np
from scipy.integrate import simps
from tqdm import tqdm

from glob import glob
import sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi')
from ToParfile_10kmdz_SBD import ReadEMC, Brocher_rho_vp

# tmpDepSRP0, tmpVsSRP0, tmpVsSRP1=ReadEMC()
# vp,rho=Brocher_rho_vp(tmpVsSRP0)


dir_CPS='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/test_CPS330/'
dir_Datta='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/SWRT_Zeng/example_modL_rhslyr7/'

def CModelSRP(dir_save,Side='L'):
    # Create Model from Schmandt et al., 2015, two locs in and outside of Snake River Plain
    LayerN=np.arange(1,201)
    Depths=np.zeros(LayerN.shape)
    for i in range(len(LayerN)):
        if i==0:
            Depths[i]=0
        else:
            if Depths[i-1]<200:
                Depths[i]=Depths[i-1]+2
            # else:
            #     Depths[i]=Depths[i-1]+25
            elif Depths[i-1]<450:
                Depths[i]=Depths[i-1]+5
            else:
                Depths[i]=Depths[i-1]+10
           
                
    tmpDepSRP, tmpVsSRP0, tmpVsSRP1=ReadEMC() #Takes a while, need improvement TBD
    if Side=='L':
        tmpVsSRP=tmpVsSRP0
    elif Side=='R':
        tmpVsSRP=tmpVsSRP1
    tmpVpSRP,tmpRhoSRP=Brocher_rho_vp(tmpVsSRP)
    VsC0=5
    VsC1=12
    # VpC0,RhoC0=Brocher_rho_vp(VsC0)
    # VpC1,RhoC1=Brocher_rho_vp(VsC1)
    VpC0=VsC0*1.73
    RhoC0=3.5
    VpC1=VsC1*1.73
    RhoC1=5
    
    Rhos=np.zeros(Depths.shape)
    Vps=np.zeros(Depths.shape)
    Vss=np.zeros(Depths.shape)
    for n,dep in enumerate(Depths):
        idv=np.where(dep>=tmpDepSRP)[0][-1]
        Vss[n]=tmpVsSRP[idv]
        Vps[n]=tmpVpSRP[idv]
        Rhos[n]=tmpRhoSRP[idv]
        if dep<900:
        # if True:
            idv=np.where(dep>=tmpDepSRP)[0][-1]
            Vss[n]=tmpVsSRP[idv]
            Vps[n]=tmpVpSRP[idv]
            Rhos[n]=tmpRhoSRP[idv]
            
            tmpn=n
        # elif dep<900:
        else:
        
            # Vss[n]=6
            # Vps[n]=Vss[n]*1.73
            # Rhos[n]=3.5
            
            # Vss[n]=6+(n-tmpn)/(200-tmpn)*(7.3-6)
            # Vps[n]=10.8+(n-tmpn)/(200-tmpn)*(13.7-10.8)
            # Rhos[n]=4.39+(n-tmpn)/(200-tmpn)*(5.57-4.39)
            
            # Vss[n]=6+(n-tmpn)/(200-tmpn)*(7.3-6)
            # Vps[n]=10.8+(n-tmpn)/(200-tmpn)*(13.7-10.8)
            # Rhos[n]=4.39+(n-tmpn)/(200-tmpn)*(5.57-4.39)
            
            # Vs 5km/s to max realistic vel
            # Vss[n]=5+(n-tmpn)/(200-tmpn)*(7.3-5)
            # Vps[n]=9.2+(n-tmpn)/(200-tmpn)*(13.7-9.2)
            # Rhos[n]=3.75+(n-tmpn)/(200-tmpn)*(5.57-3.75)
            
            # Vss[n]=5+(n-tmpn-1)/(198-tmpn)*(7.3-5)
            # Vps[n]=9.2+(n-tmpn-1)/(198-tmpn)*(13.7-9.2)
            # Rhos[n]=3.75+(n-tmpn-1)/(198-tmpn)*(5.57-3.75)
            
            # Vss[n]=5+(n-tmpn-1)/(198-tmpn)*(7-5)
            # Vps[n]=9+(n-tmpn-1)/(198-tmpn)*(13-9)
            # Rhos[n]=3.5+(n-tmpn-1)/(198-tmpn)*(5.5-3.5)
            
            # Vss[n]=VsC0+(n-tmpn)/(200-tmpn)*(VsC1-VsC0)
            # Vps[n]=VpC0+(n-tmpn)/(200-tmpn)*(VpC1-VpC0)
            # Rhos[n]=RhoC0+(n-tmpn)/(200-tmpn)*(RhoC1-RhoC0)
            
            # Vss[n]=Vss[tmpn]+(n-tmpn)/(200-tmpn)*(7.3-Vss[tmpn])
            # Vps[n]=Vps[tmpn]+(n-tmpn)/(200-tmpn)*(13.7-Vps[tmpn])
            # Rhos[n]=Rhos[tmpn]+(n-tmpn)/(200-tmpn)*(5.57-Rhos[tmpn])
            
            # Vss[n]=12
            # Vps[n]=Vss[n]*1.73
            # Rhos[n]=5
            
            Vss[n]=6.9
            Vps[n]=Vss[n]*1.73
            Rhos[n]=5
            # Vss[n]=8
            # Vps[n]=Vss[n]*1.73
            # Rhos[n]=5
            

    Model=[]
    NLines=12
    # metadata for CPS model
    with open(dir_CPS+'model96','r') as fp:
        n=0
        for line in fp:
            if n==NLines:
                break
            Model.append(line.split('\n')[0])
            n+=1
            
    # Model information  
    for n in range(len(Depths)):
        if n==len(Depths)-1:
            H=Depths[n]-Depths[n-1]
        else:
            H=Depths[n+1]-Depths[n]
        # Model.append('%.1f %.2f %.2f %.2f 0.0 0.0 0.0 0.0 1.0 1.0'%(H,Vps[n],Vss[n],Rhos[n]))  
        Model.append('%.1f %.3f %.3f %.3f 0.0 0.0 0.0 0.0 1.0 1.0'%(H,Vps[n],Vss[n],Rhos[n]))  


    np.savetxt(dir_save+'/Model'+Side,Model,fmt='%s')
    
    return Depths,Vss,Vps,Rhos

def CModelH(dir_save,Side='L'):
    # Create Model for Hermmann 2013 as input
    # Current model from Datta 2018
    LayerN=np.arange(1,201)
    Depths=np.zeros(LayerN.shape)
    for i in range(len(LayerN)):
        if i==0:
            Depths[i]=1
        else:
            if Depths[i-1]<50:
                Depths[i]=Depths[i-1]+1
            elif Depths[i-1]<215:
                Depths[i]=Depths[i-1]+5
            else:
                Depths[i]=Depths[i-1]+10
        
    Rhos=np.zeros(Depths.shape)
    Vps=np.zeros(Depths.shape)
    Vss=np.zeros(Depths.shape)
    
    for n, dep in enumerate(Depths):
        # Two sides are different
        if Side=='L':
            if dep<=30:
                Rhos[n]=2.80
                Vps[n]=6.24
                Vss[n]=3.60
                continue
        if Side=='R':
            if dep<=5:
                Rhos[n]=3.00
                Vps[n]=6.58
                Vss[n]=3.80
                continue  
        # Two sides are the same
        if dep<=50:
            Rhos[n]=3.20
            Vps[n]=8.23
            Vss[n]=4.75
        elif dep<=150:
            Rhos[n]=3.05
            Vps[n]=7.53
            Vss[n]=4.35
        elif dep<=215:
            Rhos[n]=3.05
            Vps[n]=7.88
            Vss[n]=4.55
        else:
            Rhos[n]=3.00
            Vps[n]=8.31
            Vss[n]=4.80
        
        # elif dep<300:
        #     Rhos[n]=3.00
        #     Vps[n]=8.31
        #     Vss[n]=4.80
        # else:
        #     Rhos[n]=3.00
        #     Vps[n]=8.31
        #     Vss[n]=4.80
    
    Model=[]
    NLines=12
    # metadata for CPS model
    with open(dir_CPS+'model96','r') as fp:
        n=0
        for line in fp:
            if n==NLines:
                break
            Model.append(line.split('\n')[0])
            n+=1
            
    # Model information  
    for n in range(len(Depths)):
        if n==0:
            H=Depths[n]
        else:
            H=Depths[n]-Depths[n-1]
        Model.append('%.1f %.2f %.2f %.2f 0.0 0.0 0.0 0.0 1.0 1.0'%(H,Vps[n],Vss[n],Rhos[n]))  

    np.savetxt(dir_save+'/Model'+Side,Model,fmt='%s')
    
    return Depths
    
# dirM=dir_CPS+'ex8_ModeSumWave'
# dirM=dir_CPS+'ex9_ModeSumLock'
# dirM=dir_CPS+'ex11_ModeSumUB0'
# dirM=dir_CPS+'ex13D_Vs8_Dep2000'
# dirM=dir_CPS+'ex14B_Vs6p9_Dep900_M20'
dirM=dir_CPS+'ex14A_Vs6p9_Dep900_M10'

if not os.path.isdir(dirM):
    os.makedirs(dirM)

# if not os.path.isfile(dirM+'/ModelL'):
# # if True:
#     Depths=CModelSRP(dirM,'L')
#     CModelSRP(dirM,'R')    

#%% Get Rayleigh Wave Eigenfunctions from Herrmann 2013
# freq0=0.01
# freq1=0.15
# df=0.005
# FARR=[]
# PARR=[]
# for freq in np.arange(freq0,freq1+df,df):
#     FARR.append(freq)
#     PARR.append(1/freq)
per0=10
per1=150
dp=1
FARR=[]
PARR=[]
for per in np.arange(per0,per1+dp,dp):
    PARR.append(per)
    FARR.append(1/per)

NModes=10
if not (os.path.isfile(dirM+'/SRDER_L.TXT') and os.path.isfile(dirM+'/SRDER_R.TXT')):
# if True:
    np.savetxt(dirM+'/FARR.txt',FARR,fmt='%.6f')
    np.savetxt(dirM+'/PARR.txt',PARR,fmt='%.6f')
    
    runCPS=['#!/bin/csh']
    # runCPS.append('sprep96 -M ModelL -DT 5.0 -NPTS 2 -L -R -NMOD 7 -PER 50')
    # runCPS.append('sprep96 -M ModelL -R -NMOD 7 -FREQ 0.15')
    # runCPS.append('sprep96 -M ModelL -R -NMOD 7 -FARR FARR.txt ')
    runCPS.append('sprep96 -M $1 -R -NMOD %d -PARR PARR.txt '%NModes)
    runCPS.append('sdisp96 -v')
    runCPS.append('sregn96 -HS 10 -HR 0 -DE -NOQ') #-HS 0
    runCPS.append('sdpder96 -R -TXT -K 2 -XLEN 3.0 -X0 1.0 -YLEN 4.0')
    
    np.savetxt(dirM+'/runCPS.csh',runCPS,fmt='%s')
    
    os.chdir(dirM)
    os.system('csh runCPS.csh ModelL')
    os.system('mv SRDER.TXT SRDER_L.TXT')
    os.system('csh runCPS.csh ModelR')
    os.system('mv SRDER.TXT SRDER_R.TXT')
#%%
import gzip
import matplotlib.pyplot as plt
# plt.close('all')

# Freqs=np.arange(freq0,freq1+df,df)
Pers=np.arange(per0,per1+dp,dp)


freq=0.15 # Input Frequency
effileL=dirM+'/SRDER_L.TXT'
effileR=dirM+'/SRDER_R.TXT'

eps=1e-7
# NModes=10
Modes=np.arange(0,NModes)

class EigenRead:
    def __init__(self,freq,effile):
        Hs=[]
        Vps=[]
        Vss=[]
        Rhos=[]
        Vph=np.empty(NModes); Vph[:]=np.nan
        
        flagModel=False
        flagT=False
        flagL=False
        
        #Start Reading Eigenfunction file
        fpCPS = open(effile,'r') #dirM global variable
        for line in fpCPS:
            # Read Model Information
            if 'H(km)' in line:
                flagModel=True
                continue        
            if flagModel:
                tmp=line.split()
                if tmp:#Non-empty string 
                    Hs.append(float(tmp[1]))
                    Vps.append(float(tmp[2]))
                    Vss.append(float(tmp[3]))
                    Rhos.append(float(tmp[4]))
                else:
                    Hs=np.array(Hs)
                    Vps=np.array(Vps)
                    Vss=np.array(Vss)
                    Rhos=np.array(Rhos)
                    NLayers=len(Hs) 
                    flagModel=False
                    # Initialize Eigenfunction arrays 
                    Ur=np.empty((NModes,NLayers)); Ur[:]=np.nan
                    Uz=np.empty((NModes,NLayers)); Uz[:]=np.nan
                    Tr=np.empty((NModes,NLayers)); Tr[:]=np.nan
                    Tz=np.empty((NModes,NLayers)); Tz[:]=np.nan
                continue
            # Read Eigenfunctions at particular periods on all available modes
        
            if 'RAYLEIGH WAVE' in line:
                mode=int(line.split()[-1])
                continue
            if 'T =' in line:
                tmpT=float(line.split()[2])
                if round(1/tmpT,4)==round(freq,4): #! Be careful about this
                    flagT=True
                    Vph[mode]=float(line.split()[5])
                    continue
            if flagT:
                if line.split()[0]=='1':
                    flagL=True
                    n=0
            if flagT and flagL:
                vals_line=[float(i) for i in line.split()]
                Ur[mode,n]=vals_line[1]
                Tr[mode,n]=vals_line[2]
                Uz[mode,n]=vals_line[3]
                Tz[mode,n]=vals_line[4]
                n+=1
        
                if n==NLayers:
                    flagT=False
                    flagL=False
        #End Reading Eigenfunction file
        fpCPS.close() 
        
        Ks=2*np.pi*freq/Vph #modes
        Mus=Rhos*(Vss**2) # Be Careful about the Unit!!! Consistent with original unit, just use it.
        Lambdas=Rhos*(Vps**2)-2*Mus # Be Careful about the Unit!!!
        
        # Eq 7.26 and 7.27 from Aki & Richards 2009 (Z-down positive) for derivation of Surface Wave Reflection Transmission coefficients
        # Eigenfunctions from Herrmann 2013 (Z-up positive)
        r1=-Ur #Convention Difference of Vertical direction between Herrmann 2013 and Aki & Richards 2009.
        r2=Uz
        r3=-Tr 
        r4=Tz
        
        Psi_xx=np.empty((NModes,NLayers)); Psi_xx[:]=np.nan
        for mode in Modes:
            Psi_xx[mode]=Lambdas*(r4[mode]-Ks[mode]*Lambdas*r1[mode])/(Lambdas+2*Mus)+Ks[mode]*(Lambdas+2*Mus)*r1[mode]
        
        Deps=np.array([np.sum(Hs[:i]) for i in range(len(Hs))])
        # Exclude Higher modes that does not exist
        Nmed=NModes
        for mode in Modes:
            if sum(np.isnan(r1[mode])):
                Nmed=mode
                break
        
        self.Deps=Deps
        self.Vss=Vss
        self.Vps=Vps
        self.Mus=Mus
        self.Lambdas=Lambdas
        self.Rhos=Rhos
        self.totm=NModes
        self.Ks=Ks[:Nmed]
        self.r1=r1[:Nmed] #To be consistent with format of read_egnfile_per
        self.r2=r2[:Nmed]
        self.r3=r3[:Nmed]
        self.r4=r4[:Nmed]
        self.Psi_xx=Psi_xx[:Nmed]
        self.Nmed=Nmed
        self.Vph=Vph[:Nmed]
        # self.utmat=None
        # self.ttmat=None
        
        # self.final_result(sum(~np.isnan(Ks)),True)
        


    # return r1[:Nmed], r2[:Nmed], r3[:Nmed], r4[:Nmed], Psi_xx[:Nmed], Deps,Nmed

def single_freq(freq,effileL,effileR):
    erL=EigenRead(freq,effileL)
    erR=EigenRead(freq,effileR)
    r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
    r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
    Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
    
    P=np.zeros((Nmed1,Nmed2))
    S0=np.zeros((Nmed2,1))
    LHS=np.zeros((Nmed1+Nmed2,Nmed2+Nmed1))
    RHS=np.zeros((Nmed1+Nmed2,1))
    for m in np.arange(Nmed1):
        for n in np.arange(Nmed2):
            Jm=2*simps(-r5_med1[m]*r1_med1[m]+r3_med1[m]*r2_med1[m],Deps)
            I_r5r1_12=simps(r5_med1[m]*r1_med2[n],Deps)
            I_r5r1_21=simps(r5_med2[n]*r1_med1[m],Deps)
            I_r3r2_12=simps(r3_med1[m]*r2_med2[n],Deps)
            I_r2r3_12=simps(r2_med1[m]*r3_med2[n],Deps)
            if m==0:
                Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                S0[n]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
            P[m][n]=(I_r5r1_12-I_r5r1_21+I_r3r2_12-I_r2r3_12)/Jm
    LHS[:Nmed1,:Nmed2]=P
    LHS[Nmed1:,Nmed2:]=-P.T #Transpose
    LHS[Nmed1:,:Nmed2]=np.eye(Nmed2)
    LHS[:Nmed1,Nmed2:]=np.eye(Nmed1)
    
    RHS[Nmed1:]=S0
    soln=np.linalg.solve(LHS,RHS)
    tcoef=soln[:Nmed2]
    rcoef=soln[Nmed2:]  
    return rcoef, tcoef#,LHS

class Csingle_freq:
    def __init__(self,freq,effileL,effileR):
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
        r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
        Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
        
        P=np.zeros((Nmed1,Nmed2))
        S0=np.zeros((Nmed2,1))
        LHS=np.zeros((Nmed1+Nmed2,Nmed2+Nmed1))
        RHS=np.zeros((Nmed1+Nmed2,1))
        for m in np.arange(Nmed1):
            for n in np.arange(Nmed2):
                Jm=2*simps(-r5_med1[m]*r1_med1[m]+r3_med1[m]*r2_med1[m],Deps)
                I_r5r1_12=simps(r5_med1[m]*r1_med2[n],Deps)
                I_r5r1_21=simps(r5_med2[n]*r1_med1[m],Deps)
                I_r3r2_12=simps(r3_med1[m]*r2_med2[n],Deps)
                I_r2r3_12=simps(r2_med1[m]*r3_med2[n],Deps)
                if m==0:
                    Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                    S0[n]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
                P[m][n]=(I_r5r1_12-I_r5r1_21+I_r3r2_12-I_r2r3_12)/Jm
                              
        LHS[:Nmed1,:Nmed2]=P
        LHS[Nmed1:,Nmed2:]=-P.T #Transpose
        LHS[Nmed1:,:Nmed2]=np.eye(Nmed2)
        LHS[:Nmed1,Nmed2:]=np.eye(Nmed1)
        
        RHS[Nmed1:]=S0
        soln=np.linalg.solve(LHS,RHS)
        tcoef=soln[:Nmed2]
        rcoef=soln[Nmed2:]  
        self.rcoef=rcoef
        self.tcoef=tcoef
        self.LHS=LHS
        
        

class Csingle_freq1:
    def __init__(self,freq,effileL,effileR):
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
        r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
        Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
        
        P=np.zeros((Nmed1,Nmed2))
        S0=np.zeros((Nmed2,1))

        for m in np.arange(Nmed1):
            for n in np.arange(Nmed2):
                Jm=2*simps(-r5_med1[m]*r1_med1[m]+r3_med1[m]*r2_med1[m],Deps)
                I_r5r1_12=simps(r5_med1[m]*r1_med2[n],Deps)
                I_r5r1_21=simps(r5_med2[n]*r1_med1[m],Deps)
                I_r3r2_12=simps(r3_med1[m]*r2_med2[n],Deps)
                I_r2r3_12=simps(r2_med1[m]*r3_med2[n],Deps)
                if m==0:
                    Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                    S0[n]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
                P[m][n]=(I_r5r1_12-I_r5r1_21+I_r3r2_12-I_r2r3_12)/Jm
                
                # if m==0 and n==Nmed2-1:
                #     self.I_r5r1_12=I_r5r1_12
                #     self.I_r5r1_21=I_r5r1_21
                #     self.I_r3r2_12=I_r3r2_12
                #     self.I_r2r3_12=I_r2r3_12
                #     self.Jm=Jm
                
           
        LHS=np.zeros((Nmed1+Nmed2-2,Nmed2+Nmed1-2))
        RHS=np.zeros((Nmed1+Nmed2-2,1))
    
        LHS[:Nmed1-1,:Nmed2-1]=P[:-1,:-1]
        LHS[Nmed1-1:,Nmed2-1:]=-P[:-1,:-1].T #Transpose
        LHS[Nmed1-1:,:Nmed2-1]=np.eye(Nmed2-1)
        LHS[:Nmed1-1,Nmed2-1:]=np.eye(Nmed1-1)
        
        RHS[Nmed1-1:]=S0[:-1]
        soln=np.linalg.solve(LHS,RHS)
        tcoef=soln[:Nmed2-1]
        rcoef=soln[Nmed2-1:]  
        
        # self.P=P
        # self.LHS=LHS
        # self.RHS=RHS
        self.rcoef=rcoef
        self.tcoef=tcoef
        
        # self.r1_med1=r1_med1
        # self.r2_med1=r2_med1
        # self.r3_med1=r3_med1
        # self.r4_med1=r4_med1
        # self.r5_med1=r5_med1
        # self.Deps_med1=Deps_med1
        # self.Nmed1=Nmed1
        
        # self.r1_med2=r1_med2
        # self.r2_med2=r2_med2
        # self.r3_med2=r3_med2
        # self.r4_med2=r4_med2
        # self.r5_med2=r5_med2
        # self.Deps_med2=Deps_med2
        # self.Nmed2=Nmed2
        
class Csingle_freq2:
    def __init__(self,freq,effileL,effileR):
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
        r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
        Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
        
        VphF=99 #phase velocity filter
        Nmed1=np.sum(erL.Vph<=VphF)
        Nmed2=np.sum(erR.Vph<=VphF)
        
        
        P=np.zeros((Nmed1,Nmed2))
        S0=np.zeros((Nmed2,1))

        for m in np.arange(Nmed1):
            for n in np.arange(Nmed2):
                Jm=2*simps(-r5_med1[m]*r1_med1[m]+r3_med1[m]*r2_med1[m],Deps)
                I_r5r1_12=simps(r5_med1[m]*r1_med2[n],Deps)
                I_r5r1_21=simps(r5_med2[n]*r1_med1[m],Deps)
                I_r3r2_12=simps(r3_med1[m]*r2_med2[n],Deps)
                I_r2r3_12=simps(r2_med1[m]*r3_med2[n],Deps)
                if m==0:
                    Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                    S0[n]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
                P[m][n]=(I_r5r1_12-I_r5r1_21+I_r3r2_12-I_r2r3_12)/Jm
                
                # if m==0 and n==Nmed2-1:
                #     self.I_r5r1_12=I_r5r1_12
                #     self.I_r5r1_21=I_r5r1_21
                #     self.I_r3r2_12=I_r3r2_12
                #     self.I_r2r3_12=I_r2r3_12
                #     self.Jm=Jm
                
           
        LHS=np.zeros((Nmed1+Nmed2,Nmed2+Nmed1))
        RHS=np.zeros((Nmed1+Nmed2,1))
    
    
        LHS[:Nmed1,:Nmed2]=P
        LHS[Nmed1:,Nmed2:]=-P.T #Transpose
        LHS[Nmed1:,:Nmed2]=np.eye(Nmed2)
        LHS[:Nmed1,Nmed2:]=np.eye(Nmed1)
        
        RHS[Nmed1:]=S0
        soln=np.linalg.solve(LHS,RHS)
       
        tcoef=soln[:Nmed2]
        rcoef=soln[Nmed2:]  
        
        self.rcoef=rcoef
        self.tcoef=tcoef

        self.r1_med1=r1_med1
        self.r2_med1=r2_med1
        self.r3_med1=r3_med1
        self.r4_med1=r4_med1
        self.r5_med1=r5_med1
        self.Deps_med1=Deps_med1
        self.Nmed1=Nmed1
        
        self.r1_med2=r1_med2
        self.r2_med2=r2_med2
        self.r3_med2=r3_med2
        self.r4_med2=r4_med2
        self.r5_med2=r5_med2
        self.Deps_med2=Deps_med2
        self.Nmed2=Nmed2

class Csingle_freqA:
    def __init__(self,freq,effileL,effileR):
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
        r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
        Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
        
        VphF=99 #phase velocity filter
        Nmed1=np.sum(erL.Vph<=VphF)
        Nmed2=np.sum(erR.Vph<=VphF)
        
        
        P=np.zeros((Nmed1,Nmed2))
        S0=np.zeros((Nmed2,1))
        S=np.zeros((Nmed2,Nmed1))
        
        for m in np.arange(Nmed1):
            for n in np.arange(Nmed2):
                Jm=2*simps(-r5_med1[m]*r1_med1[m]+r3_med1[m]*r2_med1[m],Deps)
                I_r5r1_12=simps(r5_med1[m]*r1_med2[n],Deps)
                I_r5r1_21=simps(r5_med2[n]*r1_med1[m],Deps)
                I_r3r2_12=simps(r3_med1[m]*r2_med2[n],Deps)
                I_r2r3_12=simps(r2_med1[m]*r3_med2[n],Deps)
                # if m==0:
                #     Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                #     S0[n]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
                Jn=2*simps(-r5_med2[n]*r1_med2[n]+r3_med2[n]*r2_med2[n],Deps)
                S[n][m]=(-I_r5r1_12-I_r5r1_21+I_r3r2_12+I_r2r3_12)/Jn
                
                P[m][n]=(I_r5r1_12-I_r5r1_21+I_r3r2_12-I_r2r3_12)/Jm
                
                # if m==0 and n==Nmed2-1:
                #     self.I_r5r1_12=I_r5r1_12
                #     self.I_r5r1_21=I_r5r1_21
                #     self.I_r3r2_12=I_r3r2_12
                #     self.I_r2r3_12=I_r2r3_12
                #     self.Jm=Jm
                
           
        LHS=np.zeros((Nmed1+Nmed2,Nmed2+Nmed1))
        RHS=np.zeros((Nmed1+Nmed2,Nmed1))
    
    
        LHS[:Nmed1,:Nmed2]=P
        LHS[Nmed1:,Nmed2:]=-P.T #Transpose
        LHS[Nmed1:,:Nmed2]=np.eye(Nmed2)
        LHS[:Nmed1,Nmed2:]=np.eye(Nmed1)
        
        RHS[Nmed1:,:]=S
        soln=np.linalg.solve(LHS,RHS)
       
        Tcoef=soln[:Nmed2,:]
        Rcoef=soln[Nmed2:,:]  
        
        tcoef=Tcoef[:,0]
        rcoef=Rcoef[:,0]
        
        self.Rcoef=Rcoef
        self.Tcoef=Tcoef
        
        self.rcoef=rcoef
        self.tcoef=tcoef

        self.r1_med1=r1_med1
        self.r2_med1=r2_med1
        self.r3_med1=r3_med1
        self.r4_med1=r4_med1
        self.r5_med1=r5_med1
        self.Deps_med1=Deps_med1
        self.Nmed1=Nmed1
        
        self.r1_med2=r1_med2
        self.r2_med2=r2_med2
        self.r3_med2=r3_med2
        self.r4_med2=r4_med2
        self.r5_med2=r5_med2
        self.Deps_med2=Deps_med2
        self.Nmed2=Nmed2

class Csingle_freqPQ:
    def __init__(self,freq,effileL,effileR):
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= erL.r1, erL.r2, erL.r3, erL.r4, erL.Psi_xx, erL.Deps, erL.Nmed
        r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= erR.r1, erR.r2, erR.r3, erR.r4, erR.Psi_xx, erR.Deps, erR.Nmed
        # Cri_dep=200
        # r1_med1=r1_med1[:,Deps_med1<Cri_dep]
        # r2_med1=r2_med1[:,Deps_med1<Cri_dep]
        # r3_med1=r3_med1[:,Deps_med1<Cri_dep]
        # r5_med1=r5_med1[:,Deps_med1<Cri_dep]
        # Deps_med1=Deps_med1[Deps_med1<Cri_dep]
        # r1_med2=r1_med2[:,Deps_med2<Cri_dep]
        # r2_med2=r2_med2[:,Deps_med2<Cri_dep]
        # r3_med2=r3_med2[:,Deps_med2<Cri_dep]
        # r5_med2=r5_med2[:,Deps_med2<Cri_dep]
        # Deps_med2=Deps_med2[Deps_med2<Cri_dep]
        
        Deps=Deps_med1 #Requires same depth sampling, OR Interpolate TBD
        
        VphF=99 #phase velocity filter
        Nmed1=np.sum(erL.Vph<=VphF)
        Nmed2=np.sum(erR.Vph<=VphF)
        
        
        LHS=np.diag(np.ones(Nmed1+Nmed2))
        RHS=np.zeros((Nmed1+Nmed2,Nmed1))
        
        Sii_A=np.empty(Nmed1)
        Sjj_A=np.empty(Nmed2)
        for i in range(Nmed1):
            Sii_A[i]=2*simps(-r5_med1[i]*r1_med1[i]+r3_med1[i]*r2_med1[i],Deps)
        for j in range(Nmed2):
            Sjj_A[j]=2*simps(-r5_med2[j]*r1_med2[j]+r3_med2[j]*r2_med2[j],Deps)
            
        for i in np.arange(Nmed1):
            for j in np.arange(Nmed2):
                
                
                I_r5r1_12=simps(r5_med1[i]*r1_med2[j],Deps)
                I_r5r1_21=simps(r5_med2[j]*r1_med1[i],Deps)
                I_r3r2_12=simps(r3_med1[i]*r2_med2[j],Deps)
                I_r3r2_21=simps(r3_med2[j]*r2_med1[i],Deps)
                
                Pij=I_r5r1_12+I_r3r2_12-I_r5r1_21-I_r3r2_21
                Sij=-I_r5r1_12+I_r3r2_12-I_r5r1_21+I_r3r2_21
                
                Sii=Sii_A[i]
                Sjj=Sjj_A[j]
                
                LHS[i,Nmed1+j]=Pij/Sii
                LHS[Nmed1+j,i]=-Pij/Sjj
                
                RHS[Nmed1+j,i]=Sij/Sjj

           
        # LHS=np.zeros((Nmed1+Nmed2,Nmed2+Nmed1))
        # RHS=np.zeros((Nmed1+Nmed2,Nmed1))
    
    
        # LHS[:Nmed1,:Nmed2]=P
        # LHS[Nmed1:,Nmed2:]=-P.T #Transpose
        # LHS[Nmed1:,:Nmed2]=np.eye(Nmed2)
        # LHS[:Nmed1,Nmed2:]=np.eye(Nmed1)
        
        # RHS[Nmed1:,:]=S
        soln=np.linalg.solve(LHS,RHS)
       
        # Tcoef=soln[:Nmed2,:]
        # Rcoef=soln[Nmed2:,:]  
        Rcoef=soln[:Nmed1,:]
        Tcoef=soln[Nmed1:,:]
        
        tcoef=Tcoef[:,0]
        rcoef=Rcoef[:,0]
        
        self.Rcoef=Rcoef
        self.Tcoef=Tcoef
        
        self.rcoef=rcoef
        self.tcoef=tcoef

        self.r1_med1=r1_med1
        self.r2_med1=r2_med1
        self.r3_med1=r3_med1
        self.r4_med1=r4_med1
        self.r5_med1=r5_med1
        self.Deps_med1=Deps_med1
        self.Nmed1=Nmed1
        
        self.r1_med2=r1_med2
        self.r2_med2=r2_med2
        self.r3_med2=r3_med2
        self.r4_med2=r4_med2
        self.r5_med2=r5_med2
        self.Deps_med2=Deps_med2
        self.Nmed2=Nmed2

#%%

dirM1=dir_CPS+'ex12B_Vs6p9_Dep900/'
dirM2=dir_CPS+'ex13A_Vs6p9_Dep900/'

dirM1=dir_CPS+'ex14A_Vs6p9_Dep900_M10/'
# dirM1=dir_CPS+'ex14B_Vs6p9_Dep900_M20/'

effileL1=dirM1+'/SRDER_L.TXT'
effileL2=dirM1+'/SRDER_R.TXT'
# effileL1=dirM1+'/SRDER_R.TXT'
# effileL2=dirM2+'/SRDER_R.TXT'

erL1=EigenRead(1/80,effileL1)
erL2=EigenRead(1/80,effileL2)
# erL1.Deps
# erL1.Vss

def DepVsPlot(erL1,dz):
    
    DepsL1=np.zeros(len(erL1.Deps)*2)
    VssL1=np.zeros(len(erL1.Deps)*2)
    for ii in range(len(erL1.Deps)):
        if ii==0:
            DepsL1[ii]=erL1.Deps[ii]
        else:
            DepsL1[2*ii-1]=erL1.Deps[ii]
            DepsL1[2*ii]=erL1.Deps[ii]
        VssL1[2*ii]=erL1.Vss[ii]
        VssL1[2*ii+1]=erL1.Vss[ii]
    DepsL1[-1]=DepsL1[-2]+dz
    return DepsL1, VssL1

DepsL1,VssL1=DepVsPlot(erL1,10)
DepsL2,VssL2=DepVsPlot(erL2,10)
if __name__=='__main__':
    plt.figure()
    plt.plot(VssL1,DepsL1,'-',color='red',markersize=2,label='Model L')
    plt.plot(VssL2,DepsL2,'-',color='blue',markersize=1,label='Model R')
    plt.legend()
    # plt.ylim([950,0])
    # plt.ylim([DepsL1[-1],0])
    # plt.ylim([200,0])
    # plt.xlim([1.3,4.7])
    plt.ylim([200,0])
    plt.xlim([2.5,5])
    plt.grid()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/1-2_ModelLR.pdf',transparent=False)
#%%
nx=200
depM=200
M_LR=np.zeros((np.sum(erL1.Deps<=depM),nx))
for i,dep in enumerate(erL1.Deps):
    if dep<=depM:
        M_LR[i][:int(nx/2)]=erL1.Vss[i]
        M_LR[i][int(nx/2):]=erL2.Vss[i]

#%%
if __name__=='__main__':
    plt.figure()
    plt.imshow(M_LR,cmap='Greys')
    plt.xticks([])
    plt.yticks([])
    # # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/VelBlock.pdf',transparent=False)

#%%

Fps=glob(dir_CPS+'ex1??_Vs*')

# for fp in tqdm(Fps):
fp=Fps[0]
suffix=fp.split('/')[-1]
suffix='ex12B_Vs6p9_Dep900'
suffix='ex8_ModeSumWave'
dirM=dir_CPS+suffix

dirM=dirM1

#%% Fig 5 Reflection and Transmission Coefficients
dirM=dir_CPS+'ex14A_Vs6p9_Dep900_M10/'
# dirM=dir_CPS+'ex12D_Vs8_Dep900'

if __name__=='__main__':
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
    tcmL=np.empty((NModes,len(Freqs)));tcmL[:]=np.nan
    rcmL=tcmL.copy()
    tcmR=tcmL.copy()
    rcmR=tcmL.copy()
    
    cArrL=tcmL.copy()
    cArrR=tcmL.copy()
    
    for idf, freq in tqdm(enumerate(Freqs)):
        # print(idf,freq)
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        cArrL[:len(erL.Vph),idf]=erL.Vph#float(erL.Vph) #[0] for fundamental mode
        cArrR[:len(erR.Vph),idf]=erR.Vph#float(erR.Vph) #[0]
        # rcoef,tcoef=single_freq(freq,effileL,effileR)
        # Csf=Csingle_freqA(freq,effileL,effileR)
        CsfL=Csingle_freqPQ(freq,effileL,effileR)
        CsfR=Csingle_freqPQ(freq,effileR,effileL)

        rcoefL=CsfL.rcoef
        tcoefL=CsfL.tcoef
        rcoefR=CsfR.rcoef
        tcoefR=CsfR.tcoef
        
        tcmL[0:len(tcoefL),idf]=tcoefL.reshape(-1)
        rcmL[0:len(rcoefL),idf]=rcoefL.reshape(-1)
        tcmR[0:len(tcoefR),idf]=tcoefR.reshape(-1)
        rcmR[0:len(rcoefR),idf]=rcoefR.reshape(-1)
        
        # Tcm[0:Csf.Tcoef.shape[1],0:Csf.Tcoef.shape[0],idf]=Csf.Tcoef.T
        # Rcm[0:Csf.Rcoef.shape[1],0:Csf.Rcoef.shape[0],idf]=Csf.Rcoef.T
#%% Fig5 Plot
if __name__=='__main__':
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": 'Times New Roman',
        "font.size": 11,
        "figure.autolayout": True}) #'axes.linewidth': 0.8
    
    XTICKS=[0,0.02,0.04,0.06,0.08,0.10]
    XLIM=[0,0.105]
    # YLIMR=[-0.008,0.022]
    # YTICKR=[-0.005,0.000,0.005,0.010,0.015,0.020]

    # YLIMT=[-0.06,1.20]
    # YTICKT=[0.0,0.2,0.4,0.6,0.8,1.0,1.2]
    
    YLIMT=[-0.8,1.20]
    # YTICKT=[-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2]
    YTICKT0=[0.8,1.2]
    YTICKT1=[-0.02,-0.01,0.0,0.01,0.02]
    
    YLIMR=[-0.03,0.03]
    YTICKR=[-0.02,-0.01,0.0,0.01,0.02]
    # YLIMR=[-0.025,0.025]
    # YTICKR=[-0.02,-0.01,0.0,0.01,0.02]

    
    kcorm=0.04 #y location coef for corner mark
    # plt.close('all')
    # fig=plt.figure(figsize=(8.5,6))
    fig=plt.figure(figsize=(8.5,7))

    # plt.tight_layout()
    # fig.tight_layout()
    ax1=plt.subplot(2,2,1)
    for i in range(NModes):
       	mnum='Mode %d' %(i)
       	plt.plot(Freqs,rcmL[i],'o-',label=mnum)
        # if i==2:
        #     break
    # plt.legend(ncol=2)
    # plt.ylabel('Vph (km/s)')
    # plt.xlabel('Freq (Hz)')
    # ax1.get_xaxis().set_visible(False)
    plt.xlim([0,0.105])
    plt.xticks(XTICKS,"")
    plt.ylim(YLIMR)
    plt.yticks(YTICKR)
    plt.ylabel('Coefficient')
    # plt.grid()
    plt.title('Reflection Coefficients (L to R)')
    plt.text(0,ax1.get_ylim()[1]*(1+kcorm)-ax1.get_ylim()[0]*kcorm,'(a)')
    
    # ax2=plt.subplot(2,2,2)
    # for i in range(NModes):
    #    	mnum='Mode %d' %(i)
    #    	plt.plot(Freqs,tcmL[i],'o-',label=mnum)
    #     # if i==2:
    #     #     break
    # # plt.legend(ncol=2)
    # # plt.ylabel('Vph (km/s)')
    # # plt.xlabel('Freq (Hz)')
    # plt.xlim([0,0.105])
    # plt.xticks(XTICKS,"")
    # plt.ylim(YLIMT)
    # plt.yticks(YTICKT)
    # plt.title('Transmission Coefficients (1 to 2)')
    # # plt.text(0,1.18,'(b)')
    # plt.text(0,ax2.get_ylim()[1]*(1+kcorm)-ax2.get_ylim()[0]*kcorm,'(b)')
    
    ax2=plt.subplot(2,2,2)
    for i in range(1):
       	mnum='Mode %d' %(i)
       	ax2.plot(Freqs,tcmL[i],'o-',label=mnum)
        # if i==2:
        #     break
    ax2.plot(XLIM,0.8*np.ones(len(XLIM)),'r--')
    ax2.set_xlim([0,0.105])
    ax2.set_xticks(XTICKS,"")
    ax2.set_ylim(YLIMT)
    ax2.set_yticks(YTICKT0)
    ax2.tick_params(axis='y',labelcolor='r')
    # plt.grid()

    ax2.set_title('Transmission Coefficients (L to R)')
    
    ax2t=ax2.twinx()
    for i in range(NModes):
        mnum='Mode %d' %(i)
        ax2t.plot(Freqs,tcmL[i],'o-',label=mnum)
    ax2t.set_ylim(YLIMR)
    ax2t.set_yticks(YTICKT1)
    # plt.legend(ncol=2)
    # plt.ylabel('Vph (km/s)')
    # plt.xlabel('Freq (Hz)')
    
    # plt.text(0,1.18,'(b)')
    ax2.text(0,ax2.get_ylim()[1]*(1+kcorm)-ax2.get_ylim()[0]*kcorm,'(b)')
    
    
    
    ax3=plt.subplot(2,2,3)
    for i in range(NModes):
       	mnum='Mode %d' %(i)
       	plt.plot(Freqs,rcmR[i],'o-',label=mnum)
        # if i==2:
        #     break
    # plt.ylabel('Vph (km/s)')
    # plt.legend(ncol=2,loc='lower right',borderpad=0.5)

    plt.xlim([0,0.105])
    plt.xticks(XTICKS)
    plt.ylim(-np.array(YLIMR)[::-1])
    plt.yticks(-np.array(YTICKR))
    plt.ylabel('Coefficient')
    plt.xlabel('Freq (Hz)')
    plt.title('Reflection Coefficients (R to L)')
    # plt.grid()

    # plt.text(0,0.009,'(c)')
    plt.text(0,ax3.get_ylim()[1]*(1+kcorm)-ax3.get_ylim()[0]*kcorm,'(c)')
    
    
    ax4=plt.subplot(2,2,4)
    # for i in range(NModes):
    #    	mnum='Mode %d' %(i)
    #    	plt.plot(Freqs,tcmR[i],'o-',label=mnum)
    #     # if i==2:
    #     #     break
    # # plt.legend(ncol=2)
    # # plt.ylabel('Vph (km/s)')
    # plt.xlim([0,0.105])
    # plt.xticks(XTICKS)
    # plt.xlabel('Freq (Hz)')
    # plt.ylim(YLIMT)
    # plt.yticks(YTICKT)
    # plt.title('Transmission Coefficients (2 to 1)')
    # # plt.text(0,1.18,'(d)')
    # plt.text(0,ax4.get_ylim()[1]*(1+kcorm)-ax4.get_ylim()[0]*kcorm,'(d)')
    
    for i in range(1):
       	mnum='Mode %d' %(i)
       	ax4.plot(Freqs,tcmR[i],'o-',label=mnum)
        # if i==2:
        #     break
    ax4.plot(XLIM,0.8*np.ones(len(XLIM)),'r--')

    ax4.set_xlim([0,0.105])
    ax4.set_xticks(XTICKS)
    ax4.set_xlabel('Freq (Hz)')
    ax4.set_ylim(YLIMT)
    ax4.set_yticks(YTICKT0)
    ax4.tick_params(axis='y',labelcolor='r')

    ax4.set_title('Transmission Coefficients (R to L)')
    
    ax4t=ax4.twinx()
    for i in range(NModes):
        mnum='Mode %d' %(i)
        ax4t.plot(Freqs,tcmR[i],'o-',label=mnum)
    ax4t.set_ylim(-np.array(YLIMR)[::-1])
    ax4t.set_yticks(YTICKT1)
    # ax4t.grid()

    plt.legend(ncol=2,loc='lower right',borderpad=0.2)

    # plt.ylabel('Vph (km/s)')
    # plt.xlabel('Freq (Hz)')
    
    # plt.text(0,1.18,'(b)')
    ax4.text(0,ax2.get_ylim()[1]*(1+kcorm)-ax2.get_ylim()[0]*kcorm,'(d)')
    
    
    TRANSPARENT=False
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig5_RT_Coefficients_R1.png',transparent=TRANSPARENT,dpi=300)
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/Figures_1Psi/Fig5_RT_Coefficinets_R1.eps',transparent=TRANSPARENT)

#%% Single Side
tmpside='L'
# tmpside='R'
# dirM=dir_CPS+'ex12A_Vs6p9_Dep600'
# dirM=dir_CPS+'ex12B_Vs6p9_Dep900'
# dirM=dir_CPS+'ex12C_Vs8_Dep600'
# dirM=dir_CPS+'ex12D_Vs8_Dep900'

# dirM=dir_CPS+'ex13A_Vs6p9_Dep900'

dirM=dir_CPS+'ex14A_Vs6p9_Dep900_M10/'
# dirM=dir_CPS+'ex14B_Vs6p9_Dep900_M20/'

if __name__=='__main__':
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
    
    for idf, freq in tqdm(enumerate(Freqs)):
        # print(idf,freq)
        erL=EigenRead(freq,effileL)
        erR=EigenRead(freq,effileR)
        cArrL[:len(erL.Vph),idf]=erL.Vph#float(erL.Vph) #[0] for fundamental mode
        cArrR[:len(erR.Vph),idf]=erR.Vph#float(erR.Vph) #[0]
        # rcoef,tcoef=single_freq(freq,effileL,effileR)
        Csf=Csingle_freqA(freq,effileL,effileR)
        # Csf=Csingle_freqPQ(freq,effileL,effileR)

        rcoef=Csf.rcoef
        tcoef=Csf.tcoef
        
        tcm[0:len(tcoef),idf]=tcoef.reshape(-1)
        rcm[0:len(rcoef),idf]=rcoef.reshape(-1)
        
        Tcm[0:Csf.Tcoef.shape[1],0:Csf.Tcoef.shape[0],idf]=Csf.Tcoef.T
        Rcm[0:Csf.Rcoef.shape[1],0:Csf.Rcoef.shape[0],idf]=Csf.Rcoef.T
        
        
#%%
if __name__=='__main__':
    # plt.close('all')
    
    ModeP=0
    for ModeP in range(10):
        if ModeP==0:
            pass
        else:
            continue
        # plt.close('all')
        plt.figure()
        for i in range(NModes):
           	mnum='Mode %d' %(i)
           	plt.plot(Freqs,Rcm[ModeP][i],'o-',label=mnum)
            # if i==2:
            #     break
        plt.legend()
        # plt.ylabel('Vph (km/s)')
        plt.xlabel('Freq (Hz)')
        plt.title('Ref Coef Mode%d'%ModeP)
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/'+Prefix+'FundModes'+tmpside+'_Ref_'+suffix+'.pdf')  
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/Ref_Mode'+str(ModeP))
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/Ref_Mode20_'+str(ModeP)+'_R')

        plt.figure()
        for i in range(NModes):
           	mnum='Mode %d' %(i)
           	plt.plot(Freqs,Tcm[ModeP][i],'o-',label=mnum)
            # if i==2:
            #     break
        plt.legend()
        # plt.ylabel('Vph (km/s)')
        plt.xlabel('Freq (Hz)')
        plt.title('Trans Coef Mode%d'%ModeP)
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/'+Prefix+'FundModes'+tmpside+'_Ref_'+suffix+'.pdf')  
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/Trans_Mode'+str(ModeP))
        # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/01062023_AllModes10/Trans_Mode20_'+str(ModeP)+'_R')

    

#%%
Prefix='Fund_'
if __name__=='__main__':
    plt.figure()
    for i in range(NModes):
       	mnum='Mode %d' %(i)
       	plt.plot(Freqs,cArrL[i],'o-',label=mnum)
        # if i==2:
        #     break
    plt.legend()
    plt.ylabel('Vph (km/s)')
    plt.xlabel('Freq (Hz)')
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/'+Prefix+'FundModes'+tmpside+'_Vph_'+suffix+'.pdf')
    
    plt.figure()
    for i in range(NModes):
       	mnum='Mode %d' %(i)
       	plt.plot(Freqs,rcm[i],'o-',label=mnum)
        # if i==2:
        #     break
    plt.legend()
    # plt.ylabel('Vph (km/s)')
    plt.xlabel('Freq (Hz)')
    plt.title('Ref Coef')
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/'+Prefix+'FundModes'+tmpside+'_Ref_'+suffix+'.pdf')  
    
    plt.figure()
    for i in range(NModes):
       	mnum='Mode %d' %(i)
       	plt.plot(Freqs,tcm[i],'o-',label=mnum)
        # if i==2:
        #     break
        
    plt.legend()
    # plt.ylabel('Vph (km/s)')
    plt.xlabel('Freq (Hz)')
    plt.title('Trans Coef')
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/12022022_AGU2022/'+Prefix+'FundModes'+tmpside+'_Trans_'+suffix+'.pdf')
#%%
'''
if __name__=='__main__':
    per=52
    # dirM=dir_CPS+'ex12A_Vs6p9_Dep600'
    dirM=dir_CPS+'ex12B_Vs6p9_Dep900'
    # dirM=dir_CPS+'ex12C_Vs8_Dep600'
    # dirM=dir_CPS+'ex12D_Vs8_Dep900'
    # dirM=dir_CPS+'ex13A_Vs6p9_Dep900'
    # dirM=dir_CPS+'ex13E_Vs6p9_Dep1300'
    # dirM=dir_CPS+'ex13B_Vs6p9_Dep2000'
    
    
    effileL=dirM+'/SRDER_L.TXT'
    effileR=dirM+'/SRDER_R.TXT'
    
    
    # Csf=Csingle_freq2(1/per,effileL,effileR)
    Csf=Csingle_freq2(1/per,effileR,effileL)
    
    n=Csf.Nmed2-1
    # n=4
    plt.figure()
    plt.plot(Csf.Deps_med2,Csf.r1_med2[n],label='r1')
    plt.plot(Csf.Deps_med2,Csf.r2_med2[n],label='r2')
    plt.plot(Csf.Deps_med2,Csf.r3_med2[n],label='r3')
    plt.plot(Csf.Deps_med2,Csf.r4_med2[n],'--',label='r4')
    plt.plot(Csf.Deps_med2,Csf.r5_med2[n],label='r5')
    plt.axvline(x=600,color='k',linestyle=':')
    plt.legend()
    plt.title('Period=%d, Mode=%d'%(per,n))
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/11112022_LockedUB0/Eig_'+dirM.split('/')[-1]+'M%d'%n)
'''
#%% Continuity of Displacement and traction at structural contrast, ignoring body wave contribution
if False:
    freq=0.01
    r1_med1,r2_med1,r3_med1,r4_med1,r5_med1, Deps_med1, Nmed1= EigenRead(freq,effileL)
    r1_med2,r2_med2,r3_med2,r4_med2,r5_med2, Deps_med2, Nmed2= EigenRead(freq,effileR)
    
    rcoef, tcoef=single_freq(freq)
    print(rcoef,tcoef)
    # plt.figure()
    # plt.imshow(LHS)
    # plt.figure()
    plt.close('all')
    plt.figure(figsize=(10,8))
    plt.subplot(2,2,1)
    plt.plot(r2_med1[0],Deps_med1,'b-',label='inc')
    plt.plot(np.sum(rcoef*r2_med1,0),Deps_med1,'r-',label='refl')
    plt.plot(np.sum(tcoef*r2_med2,0),Deps_med1,'g-',label='trans')
    plt.plot(r2_med1[0]-np.sum(rcoef*r2_med1,0)-np.sum(r2_med2*tcoef,0),Deps_med1,'k-',label='inc+refl-trans')
    plt.grid()
    plt.ylim([1375,0])
    plt.legend()
    plt.ylabel('Depth (km)')
    # plt.xlabel('Relative Amplitude')
    plt.title('Uz')
    
    # plt.figure(figsize=(10,8))
    plt.subplot(2,2,2)
    plt.plot(r3_med1[0],Deps_med1,'b-',label='inc')
    plt.plot(np.sum(rcoef*r3_med1,0),Deps_med1,'r-',label='refl')
    plt.plot(np.sum(tcoef*r3_med2,0),Deps_med1,'g-',label='trans')
    plt.plot(r3_med1[0]+np.sum(rcoef*r3_med1,0)-np.sum(r3_med2*tcoef,0),Deps_med1,'k-',label='inc+refl-trans')
    plt.grid()
    plt.ylim([1375,0])
    plt.legend()
    # plt.ylabel('Depth (km)')
    # plt.xlabel('Relative Amplitude')
    plt.title('Txz')
    
    
    plt.subplot(2,2,3)
    plt.plot(r1_med1[0],Deps_med1,'b-',label='inc')
    plt.plot(np.sum(rcoef*r1_med1,0),Deps_med1,'r-',label='refl')
    plt.plot(np.sum(tcoef*r1_med2,0),Deps_med1,'g-',label='trans')
    plt.plot(r1_med1[0]+np.sum(rcoef*r1_med1,0)-np.sum(r1_med2*tcoef,0),Deps_med1,'k-',label='inc+refl-trans')
    plt.grid()
    plt.ylim([1375,0])
    plt.legend()
    plt.ylabel('Depth (km)')
    plt.xlabel('Relative Amplitude')
    plt.title('Ux')
    
    plt.subplot(2,2,4)
    plt.plot(r5_med1[0],Deps_med1,'b-',label='inc')
    plt.plot(np.sum(rcoef*r5_med1,0),Deps_med1,'r-',label='refl')
    plt.plot(np.sum(tcoef*r5_med2,0),Deps_med1,'g-',label='trans')
    plt.plot(r5_med1[0]-np.sum(rcoef*r5_med1,0)-np.sum(r5_med2*tcoef,0),Deps_med1,'k-',label='inc+refl-trans')
    plt.grid()
    plt.ylim([1375,0])
    plt.legend()
    # plt.ylabel('Depth (km)')
    plt.xlabel('Relative Amplitude')
    plt.title('Txx')
