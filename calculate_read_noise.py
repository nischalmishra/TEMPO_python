# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:11:23 2017

@author: nmishra
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

The purpose of this code is to do some analysis on how to calculate the variance
of co-added frames. Variance calculation is the integral part of many tempo
requirement verification. Note this will also help calculate read noise

"""

import os
import numpy as np
import scipy.io
#from scipy.io.idl import readsav
import matplotlib.pyplot as plt
import random
from random import randint
from scipy.io.idl import readsav                                  
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter                                  

                                 
                                 

def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def filter_outlier_median(quads, ndims, nx_quad, ny_quad):    
    sigma = 5.
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < sigma]
    return outlier_filtered_data

def perform_bias_subtraction (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    ndims, nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[[0]*ndims]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[:, :, ::2]
    avg_bias_even = np.mean(even_detector_bias, axis=2)  
    odd_detector_bias = trailing_overclocks[:, :, 1::2]
    avg_bias_odd = np.mean(odd_detector_bias, axis=2)
    even_detector_active_quad = active_quad[:, :, ::2]     
    odd_detector_active_quad = active_quad[:, :, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, :,None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, :, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ndims, ny_quad, nx_quad))    
    bias_subtracted_quad[:,:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:,:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad

def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPA_Bias_Vs_Temp\-25C_Darks_Data\Script_Data\saved_quads'
    data_path_all = [each for each in os.listdir(file_path) if each.endswith('_0.dat.sav')]   
    #file_to_check = random.choice(data_path_all)
    #file_to_check = r'SPI_20160916022643568_PT_VS_TEMPERATURE_-20C_LIGHT_210_FT6_NOM_INT_99999.dat.sav'
    file_to_check = data_path_all[0]
    print(file_to_check)
    #cc
    data_path_name_split = file_to_check.split('_')    
    int_time = round(int(data_path_name_split[-1].split('.')[0])) #for integ_sweep 
    #temp = data_path_name_split[5]
    #print(temp)
    data_file = os.path.join(file_path, file_to_check)
    IDL_variable = readsav(data_file)    
    all_quads_frame = IDL_variable.q  
    quads = ['QuadA', 'QuadB',' QuadC', 'QuadD' ]
    color = ['blue','red','green','magenta']
    for i in range(0, 4):
        quad_D_all_frame =  all_quads_frame[:,int(i),:,:]
        active_quad_D = quad_D_all_frame[:, 4:1028, 10:1034] 
        bias_val_D = quad_D_all_frame[:, 4:1028, 1034:1056]
        avg_bias_D = np.mean(bias_val_D, axis=2)
        bias_subtracted_quad_D = active_quad_D - avg_bias_D[:,:, None]    
        
             
        
        # let us plot the histogram of the variance of each of the active pixels
        label1 = str(quads[i])
        variance_all_active_pixels = np.var(bias_subtracted_quad_D, axis=0)
        label2 = 'mean variance = '+str(round(np.mean(variance_all_active_pixels),3))
        print(round(np.mean(variance_all_active_pixels),3))
        nx_quad, ny_quad = variance_all_active_pixels.shape
        plt.figure(figsize= (10,8))
        plt.hist(np.reshape(variance_all_active_pixels,(nx_quad*ny_quad,1)), 30, normed=0, facecolor=color[i],alpha=1, label=label2)
        title = 'Histogram of Variance of' +label1+' (' + 'int.time = ' +str(int_time) + ' microsecs)' 
        plt.title(title,fontsize=14, fontweight="bold" )
        plt.grid(True, linestyle=':')
        plt.xlabel('Temporal Noise (Variance)', fontsize=14, fontweight="bold")
        plt.ylabel('Frequency', fontsize=14, fontweight="bold")
        plt.legend(loc='best')
        ax = plt.gca()
        #plt.xlim(0, 300)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPA_Bias_Vs_Temp\-25C_Darks_Data\Script_Data\saved_quads\read_noise'
        if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
        plt.savefig(save_dir+'/'+quads[i]+'_all_variance.png', dpi=100)
        plt.close('all')
        
        
           
        
        
        
    
    
    
if __name__ == "__main__":
    main()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        