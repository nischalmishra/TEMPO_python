# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:47:56 2019

@author: nmishra
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import pandas as pd


# How does RN change with coadds

#data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC\saved_quads\saved_hdf_input_files\Read_Noise'
#
#data_path_all = sorted([each for each in os.listdir(data_dir)
#                            if each.endswith('.h5')])
#all_data = len(data_path_all)
#
#quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
#
#
#COLOR = ['Red', 'Blue', 'Red', 'Magenta']
#for k in range(0, 4):
#    all_tsoc = []
#    for i in range( 0, len(data_path_all)):
#       
#        hdf_name = os.path.join(data_dir, data_path_all[i])
#        file = h5py.File(hdf_name, 'r')
#        full_frame = file.get(quads[k])    
#        tsoc = full_frame[2:1030, 1034:1056]
#          
#        
#        even_detector_bias = tsoc[:, 1::2]
#        even_detector_bias = even_detector_bias[:, 4:-1]
#        #print(np.mean(even_detector_bias))
#        
#        all_tsoc.append(even_detector_bias)
#        
#    all_tsoc = np.array(all_tsoc)
#    
#    
#    num_coadds = all_tsoc.shape[0]
#    all_RN = []
#    for j in range(0, num_coadds):
#        tsoc_RN = sum(all_tsoc[0:j+1, :, :])/j+1  
#        print(tsoc_RN.shape)
#        
#        RN = np.var(tsoc_RN)
#        mean_RN = np.sqrt(np.mean(RN))
#        all_RN.append(mean_RN)
#        
#    #all_RN = np.array(all_RN)
#    #print(all_RN.shape)
#    
#    plt.plot(list(range(1, 101)), all_RN,'b.', color=COLOR[k], label= quads[k])
#    plt.title('Read Noise Estimates Vs. Coadds ', fontsize=12)
#    plt.ylabel('Readn Noise Estimate (RN)',  fontsize=12)
#    plt.xlabel('Number of Coadds')
#    plt.legend(loc='best')
#    plt.grid(linestyle=':')
   
# How does read noise change with time
    
    
data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\saved_hdf_input_files\Dark'


# read the fpe temp
dframe = pd.read_excel(r'C:\Users\nmishra\Desktop\FPE_Temp_TVAC_VIS_CAL.xlsx', sheetname='Sheet2',
                        index_col=None, header=None)

fpe_temp = dframe.values
data_path_all = ([each for each in os.listdir(data_dir)
                            if each.endswith('.h5')])
all_data = len(data_path_all)
quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']

COLOR = ['Red', 'Blue', 'Green', 'Magenta']
for k in range(0, 4):
    all_RN = []
    for i in range( 0, len(data_path_all)):    
        hdf_name = os.path.join(data_dir, data_path_all[i])
        file = h5py.File(hdf_name, 'r')
        full_frame = file.get(quads[k])    
        tsoc = full_frame[2:1030, 1034:1056]   
        even_detector_bias = tsoc[:, 1::2]
        even_detector_bias = even_detector_bias[:, 4:-1]/63    
        all_tsoc = np.array(even_detector_bias)
        #all_tsoc = all_tsoc[:, :, 2:]
        RN = np.var(all_tsoc)
        
        all_RN.append(np.sqrt(RN))
    
    #plt.figure(); 
    plt.plot(fpe_temp, all_RN,'b.', label= quads[k], color=COLOR[k])
    text1 = 'RN = ' + str(round(np.sqrt(np.mean(RN)),3))+' DN'
    #plt.hist((RN), 100, label= text1,  facecolor='m', alpha=0.75)
    plt.title('Read Noise Estimate Vs FPE Temperature\n (TEMPO@TVac, 93 ms int. time, 63 coadds)', fontsize=12)
    plt.ylabel('Read Noise (DN)',  fontsize=12)
    plt.xlabel('FPA Temperature (Deg C)')
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.show()


    
    
    

