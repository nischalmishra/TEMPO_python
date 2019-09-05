# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:47:56 2019

@author: nmishra
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt

data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC\saved_quads\saved_hdf_input_files\Read_Noise'

data_path_all = sorted([each for each in os.listdir(data_dir)
                            if each.endswith('.h5')])
all_data = len(data_path_all)

quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
all_tsoc = []
coadds=1
for i in range( 0, len(data_path_all), coadds):
    if coadds==1:
         hdf_name = os.path.join(data_dir, data_path_all[i])
         file = h5py.File(hdf_name, 'r')
         full_frame = file.get(quads[3])    
         tsoc = full_frame[2:1030, 1034:1056]
    elif coadds==2:
        hdf_name1  = os.path.join(data_dir, data_path_all[i])
        file1 = h5py.File(hdf_name1, 'r')
        full_frame1 = file1.get(quads[3])    
        tsoc1 = full_frame1[2:1030, 1034:1056]
        hdf_name2  = os.path.join(data_dir, data_path_all[i+1])
        file2 = h5py.File(hdf_name2, 'r')
        full_frame2 = file2.get(quads[3])    
        tsoc2 = full_frame1[2:1030, 1034:1056]
        tsoc = (tsoc1+tsoc2)/2   
    
    even_detector_bias = tsoc[:, 1::2]
    even_detector_bias = even_detector_bias[:, 4:-1]
    #print(np.mean(even_detector_bias))
    
    all_tsoc.append(even_detector_bias)
    
all_tsoc = np.array(all_tsoc)
all_tsoc = all_tsoc[:, :, 2:]
print(all_tsoc.shape)

RN = np.var(all_tsoc, axis=0)
#print(RN.shape)
#cc
#RN = np.reshape(RN,[RN.shape[0]*RN.shape[1]])

text1 = 'RN = ' + str(round(np.sqrt(np.mean(RN)),3))+' DN'
plt.hist((RN), 100, label= text1,  facecolor='m', alpha=0.75)
plt.title('Histogram of variance of Trailing Serial Overclocks (TVac), ' + quads[3], fontsize=12)
plt.ylabel('Frequncy',  fontsize=12)
plt.xlabel('Variance')
plt.legend(loc='best')
plt.grid(linestyle=':')
plt.show()

    
    
    

