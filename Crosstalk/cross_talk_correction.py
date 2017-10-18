# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 08:43:30 2017

@author: nmishra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_outlier_mask():
    outlier_mask= np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_outlier_mask_2_sigma.csv', delimiter=',')
    quad_A = outlier_mask[0:1024, 0:1024]
    quad_B = outlier_mask[1024:, 0:1024]
    quad_C = outlier_mask[1024:, 1024:]
    quad_D = outlier_mask[0:1024:, 1024:]
    outlier_mask_final = [quad_A, quad_B, quad_C, quad_D]
    return outlier_mask_final
    
    """
    This function reads the outlier_mask
    """
    


def perform_bias_subtraction_ave (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, ::2]
    # remove outliers 
    # First 4 hot lines in even and odd
    # last odd lne in odd
    even_detector_bias = even_detector_bias[:, 4:]
    avg_bias_even = np.mean(even_detector_bias, axis=1)  
    odd_detector_bias = trailing_overclocks[:, 1::2]
    odd_samples = odd_detector_bias[:, 4:]
    rows, cols = odd_samples.shape   
    odd_detector_bias = odd_samples[:, 0:cols-1] 
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]     
    odd_detector_active_quad = active_quad[:, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (nx_quad, ny_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad
    
   
def perform_smear_subtraction(active_quad, int_time):
    # the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    tFT = 8*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)
    #print(smear_factor.shape)
    #cc
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad 
    

def create_image(image_data, title):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='bwr', origin='lower')    
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)    
    plt.show()    
    plt.close('all')   
    

  
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Bar_target\Integration_Sweep\op_int'
    save_file_path = r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test' 
    outlier_mask = read_outlier_mask()    
    all_int_files = [each for each in os.listdir(file_path) \
                     if each.endswith('dat.sav')]  
       
    for data_files in all_int_files:
            data_file = os.path.join(file_path, data_files)            
            IDL_variable = readsav(data_file)
            data_path_name_split = data_files.split('_')
            int_time = int(data_path_name_split[-1].split('.')[0]) 
            all_full_frame = IDL_variable.q            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
            string1 = 'Integ_time_'
            string2 = 'Int.time = '
            all_active_quad =  []
            for i in range(0, 4):
                quad_full_frame = all_full_frame[:, i, :, :]
                avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
                active_quad = avg_quad[4:1028, 10:1034]
                tsoc = avg_quad[4:1028, 1034:1056]  
                bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc)              
                mask_array = np.ma.masked_array(bias_subtracted_quad, mask = outlier_mask[i])
                smear_subtracted_quad = perform_smear_subtraction(mask_array, int_time)
                all_active_quad.append(smear_subtracted_quad)
            
            quad_A = all_active_quad[0]
            quad_B = all_active_quad[1]
            quad_C = all_active_quad[2]
            quad_D = all_active_quad[3]
            
            cross_talk_corrected_D = quad_D - 0.0014*quad_C
            title = 'Quad D Image (Before cross talk correction)'
            create_image(quad_D, title)
            title = 'Quad D Image (After cross talk correction)'
            create_image(cross_talk_corrected_D, title)
            
                
if __name__ == "__main__":
    main()