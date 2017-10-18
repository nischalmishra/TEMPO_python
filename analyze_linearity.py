# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""

import copy
import os
import numpy as np
import csv
import pandas as pd
#from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from make_quads_from_raw_fpe import make_quads_from_IDL_sav_file,\
                                    make_full_frame_from_raw_fpe,\
                                    full_frame_to_quads,\
                                    find_active_region_each_quads,\
                                    find_trailing_overclock_pixels,\
                                    find_leading_overclock_pixels,\
                                    find_smear_pixels,\
                                    bias_corrected_quads
from analytical_functions import plot_full_frame_image,\
                                 plot_each_quad,\
                                 plot_hist_image,\
                                 plot_hist_each_quad
                                 #calculate_std_dev,\
                                 #calculate_fft_full_frame,\
                                 #calculate_fft_each_quad,\

from outlier_detection import identify_saturation,\
                              reject_outlier_median,\
                              reject_outlier_mean,\
                              create_outlier_mask,\
                              create_final_mask
def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads'
    data_path_all = [each for each in os.listdir(file_path) if each.endswith('.sav')]
    int_time = [ ]
    all_med_quad_A = [ ]
    all_med_quad_B = [ ]
    all_med_quad_C = [ ]
    all_med_quad_D = [ ]
    unct_lin_quad_A = [ ]
    unct_lin_quad_B = [ ]
    unct_lin_quad_C = [ ]
    unct_lin_quad_D = [ ]
    dframe1 = pd.DataFrame()
    dframe2 = pd.DataFrame()
    #names_to_check = ['SHORT', 'NOM', 'LONG']
    names_to_check = ['SHORT', 'NOM','OP', 'LONG']
    for data_path in data_path_all:
        print(data_path)
        all_full_frame, full_frame, collection_type, frames, integration_time = \
        make_quads_from_IDL_sav_file(file_path, data_path, names_to_check)
        int_time.append(integration_time)
        
        # Let's start the analysis for quad A. This analysis is based on Ball
        # linearity SER
        #let's find the uncertainty
        #Steps include finding uncertainty of the mean of 100 frames
        quad_A = all_full_frame[:,0,:,:]
        active_quad_A = quad_A[:, 4:1028, 10:1034]
        quad_A_active = copy.deepcopy(active_quad_A)
        bias_A = quad_A[:, 4:1028, 1034:1056]
        avg_bias_A = np.mean(bias_A, axis=2)
        bias_subtracted_quad_A = quad_A_active - avg_bias_A[:,:,None]
        median_quad_A = np.median(bias_subtracted_quad_A)   
        stdev_quad_A = np.mean(np.apply_over_axes(np.std, bias_subtracted_quad_A ,(1,2)))/10
        unct_quad_A = (stdev_quad_A)      
          
        fractional_unct_quad_A = unct_quad_A/median_quad_A        
        all_med_quad_A.append(median_quad_A)
        #unct_lin_quad_A.append(fractional_unct_quad_A )
        #print(all_med_quad_A)
        #for dark current
        unct_lin_quad_A.append(unct_quad_A)
       
        
        
        # For quad B
        quad_B = all_full_frame[:,1,:,:]
        active_quad_B = quad_B[:, 4:1028, 10:1034]
        bias_B = quad_B[:, 4:1028, 1034:1056]
        avg_bias_B = np.mean(bias_B, axis=2)
        bias_subtracted_quad_B = active_quad_B - avg_bias_B[:,:,None]
        median_quad_B = np.median(bias_subtracted_quad_B)        
        all_med_quad_B.append(median_quad_B)
         
        unct_quad_B = (stdev_quad_B)
        fractional_unct_quad_B = unct_quad_B/median_quad_B        
        #all_med_quad_B.append(median_quad_B)
        #unct_lin_quad_B.append(fractional_unct_quad_B )
        #for dark current
        unct_lin_quad_B.append(unct_quad_B)
        #print(all_med_quad_B)
       
        
        
        # For quad C
        quad_C = all_full_frame[:,2,:,:]
        active_quad_C = quad_C[:, 4:1028, 10:1034]
        bias_C = quad_C[:, 4:1028, 1034:1056]
        avg_bias_C = np.mean(bias_C, axis=2)
        bias_subtracted_quad_C = active_quad_C - avg_bias_C[:,:,None]
        median_quad_C = np.median(bias_subtracted_quad_C)
        
        all_med_quad_C.append(median_quad_C)        
        stdev_quad_C = np.mean(np.apply_over_axes(np.std,bias_subtracted_quad_C ,(1,2)))/10
        unct_quad_C = (stdev_quad_C)       
        fractional_unct_quad_C = unct_quad_C/median_quad_C
        
        #all_med_quad_C.append(median_quad_C)
        #unct_lin_quad_C.append(fractional_unct_quad_C )
        #for dark current
        unct_lin_quad_C.append(unct_quad_C)
        #print(all_med_quad_C)
       
        
        
        quad_D = all_full_frame[:,3,:,:]
        active_quad_D = quad_D[:, 4:1028, 10:1034]
        bias_D = quad_D[:, 4:1028, 1034:1056]
        avg_bias_D = np.mean(bias_D, axis=2)
        bias_subtracted_quad_D = active_quad_D - avg_bias_D[:,:,None]
        median_quad_D = np.median(bias_subtracted_quad_D)   
        stdev_quad_D = np.mean(np.apply_over_axes(np.std,bias_subtracted_quad_D , (1, 2)))/10
        unct_quad_D = (stdev_quad_D)/10
        fractional_unct_quad_D = unct_quad_D/median_quad_D        
        all_med_quad_D.append(median_quad_D)
        #unct_lin_quad_D.append(fractional_unct_quad_D )
        #for dark current
        unct_lin_quad_D.append(unct_quad_D)
#    print('A:', all_med_quad_A)
#    print('B:', all_med_quad_B)
#    print('C:',  all_med_quad_C)
#    print('D:', all_med_quad_D)
      
    
            
    dframe1 = pd.DataFrame(
                {'Integration' : int_time,
                 'Avg. Quad A' : all_med_quad_A,
                 'Avg. Quad B' : all_med_quad_B,
                 'Avg. Quad C' : all_med_quad_C,
                 'Avg. Quad D' : all_med_quad_D,
                 })
    dframe2 = pd.DataFrame(
                {'Integration' : int_time,
                 'Unct. Quad A' :  unct_lin_quad_A,
                 'Unct. Quad B' :  unct_lin_quad_B,
                 'Unct. Quad C' :  unct_lin_quad_C,
                 'Unct. Quad D' :  unct_lin_quad_D,
                 })
    
    dframe1.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Dark_Current_Stability\Dark_Current_mean.csv')
    dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Dark_Current_Stability\Dark_Current_Unct.csv')
if __name__ == "__main__":
    main()