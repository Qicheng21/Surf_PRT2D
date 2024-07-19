#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:53:15 2022

@author: u1318104
"""
import os 
from glob import glob
import numpy as np

# dir_fp='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/Rayleigh_wave_SRP/SK_Vs0perc/'
dir_fp='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/RayleighWaveGJI/'
os.chdir(dir_fp)

# FPs=glob('O*SOURCE*SideR')
# FPs=glob('OUTPUT_SRP_SK_SOURCE*_SBD*km')
# FPs=glob('OUTPUT_SRP_SOURCE*_Mw*')
# FPs=glob('OUTPUT_SRP_SOURCE*_X*kmZ*km')
# FPs=glob('OUTPUT_SRP_SOURCE*_X6000kmZ600km_v3*')
# FPs=glob('OUTPUT_SRP_SOURCE*_Mw10km_NormalModes')
# FPs=glob('OUTPUT_SRP_SOURCE*_Mode*_WithRef')
# FPs=glob('OUTPUT_SRP_SOURCE*_LockMode*_WithRef')
# FPs=glob('OUTPUT_SRP_SOURCE*_LockModeAll_WithRef')
# FPs=glob('OUTPUT_SRP_SOURCE*_Lock_ex*_Mode*_WithRef10_10_test4*')
# FPs=glob('OUTPUT_SRP_SOURCE*_Lock_ex*_Mode*_WithRef_test*_*RTcm')
# FPs=glob('OUTPUT_SRP_SOURCE*_Mw10km_X12000kmZ1200km')
# FPs=glob('OUTPUT_SRP_SOURCE*sem1d12kkm_*')
FPs=glob('OUTPUT_SRP_*200km*')



dir_SK='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/specfem2d/EXAMPLES/RayleighWaveGJI/'
Stas,Network = np.loadtxt(dir_SK+'/DATA/STATIONS',dtype=str,usecols=[0,1],unpack=True)
#%%
# dir_FTAN='/uufs/chpc.utah.edu/common/home/flin-group5/qicheng/ANT_Adj/aftani_c_pgl_output_amp.Utah/'
runsac=['sac <<END']

# nR=0 #read data files
# MaxR=9 #maximum read files
# for sta in Stas:
#     if nR==0:
#         runsac.append('r ')
#     if nR< MaxR:
#         runsac[-1]+=sta+'.EHZ.sac '
#         nR+=1
#     else:
#         runsac[-1]+=sta+'.EHZ.sac '
#         runsac.append('ch b -60')
#         runsac.append('wh')
#         nR=0
for i in range(5):
    runsac.append('r S'+str(i)+'???.EHZ.sac')
    runsac.append('ch b -60')
    # runsac.append('ch b 0')
    runsac.append('wh')                          

    # runsac.append('r *EHR.sac')
    runsac.append('r S'+str(i)+'???.EHR.sac')
    runsac.append('ch b -60')
    # runsac.append('ch b 0')
    runsac.append('wh')
runsac.append('quit')
runsac.append('END')

for fp in FPs:
    os.chdir(dir_fp)

    np.savetxt(fp+'/SAC/'+'touchsac.csh',runsac,fmt='%s')
    
    os.system('cp ./tool_FTAN_SNR_EqFC.csh '+fp+'/SAC/')
    
    os.chdir(dir_fp+fp+'/SAC/')
    os.system('csh touchsac.csh')
    os.system('csh tool_FTAN_SNR_EqFC.csh &')
    