# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:22:32 2017

@author: nmishra
"""
import numpy as np
import matplotlib.pyplot as plt
outlier_mask = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_outlier_mask.csv', delimiter=',')

rows, cols = outlier_mask.shape
mask_a = outlier_mask[0: int(rows/2), 0:int(cols/2)]
mask_b = outlier_mask[int(rows/2):int(rows), 0:int(cols/2)]
mask_d = outlier_mask[0:int(rows/2), int(cols/2):int(cols)]
mask_c = outlier_mask[int(rows/2):int(rows), int(cols/2):int(cols)]

# In order to get it to spectrometer config do the following steps
#final_mask_a = mask_a
#final_mask_b = np.fliplr(mask_b)
##final_mask_c = np.rot90(mask_c, 2)
#final_mask_d = np.flipud(mask_d)

lower_quad = np.concatenate((mask_d , np.fliplr(mask_c)), axis=1)
upper_quad = np.concatenate((np.flipud(mask_a), np.rot90(mask_b, 2)), axis=1)
outlier_mask_full_frame_spectrometer = np.concatenate((lower_quad, upper_quad), axis=0)
plt.imshow(outlier_mask_full_frame_spectrometer, cmap='bwr', origin='lower')
print(outlier_mask_full_frame_spectrometer.shape)
plt.show()
cc

final_mask_name = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_spectrometer_outlier_mask.csv'
np.savetxt(final_mask_name, np.array(outlier_mask_full_frame_spectrometer), delimiter=",")


            