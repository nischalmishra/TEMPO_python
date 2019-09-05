# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 08:20:10 2019

@author: nmishra
"""

import numpy as np
import h5py

hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
file = h5py.File(hdf_name, 'r')
prnu = file.get('prnu')
prnu = np.array(prnu).transpose()
quad_d = prnu[2:1030, 10:1034]
quad_c = prnu[2:1030, 1078:2102]
quad_a = prnu[1062:2090, 10:1034]
quad_b = prnu[1062:2090, 1078:2102]

prnu_map_UV = np.concatenate((quad_d, quad_c), axis=1)
prnu_map_VIS = np.concatenate((quad_a, quad_b), axis=1)
prnu_map = np.concatenate((prnu_map_VIS, prnu_map_UV), axis=0)

prnu_file_name = r'C:\Users\nmishra\Desktop\SAO_deliverable\TEMPO_PRNU.csv'
np.savetxt(prnu_file_name, prnu_map,fmt='%1.5f',  delimiter=",")

cc

file_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\Outlier_median\save_masks_logical_OR\to_SAO_init'
mask_a = np.genfromtxt(file_path + '/' + 'quad_A_outlier_mask.csv',
                           delimiter=',')
mask_b = np.genfromtxt(file_path + '/' + 'quad_B_outlier_mask.csv',
                           delimiter=',')
mask_c = np.genfromtxt(file_path + '/' + 'quad_C_outlier_mask.csv',
                           delimiter=',')
mask_d = np.genfromtxt(file_path + '/' + 'quad_D_outlier_mask.csv',
                           delimiter=',')

outliers_UV =  np.concatenate((mask_d, mask_c), axis=1)
outliers_VIS = np.concatenate((mask_a, mask_b), axis=1)
outlier_mask = np.concatenate((outliers_VIS, outliers_UV), axis=0)

outlier_mask_name = r'C:\Users\nmishra\Desktop\SAO_deliverable\Outlier_mask.csv'
np.savetxt(outlier_mask_name,outlier_mask,fmt='%d',  delimiter=",")
