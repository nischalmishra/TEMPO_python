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
from scipy.io.idl import readsav
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

def filter_outlier_median(quads):
    ndims, nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 5]
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
    
     
    
    
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPA_Bias_Vs_Temp'
    temp_folder = os.listdir(file_path)
    data_folder = r'Script_Data\saved_quads'
    for files in range(0, len(temp_folder)):        
        saved_quads = os.path.join( file_path, temp_folder[files], data_folder)
        
        nominal_int_files = [each for each in os.listdir(saved_quads) \
                             if each.endswith('.sav')]
        save_dir = r'quad_images'
        plot_dir = os.path.join(saved_quads, save_dir)
        all_int_time = [ ]
        all_med_quad_A = [ ]
        all_med_quad_B = [ ]
        all_med_quad_C = [ ]
        all_med_quad_D = [ ]
        all_var_quad_A = [ ]
        all_var_quad_B = [ ]
        all_var_quad_C = [ ]
        all_var_quad_D = [ ]
        dframe1 = pd.DataFrame()
       
        
        if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        for data_files in nominal_int_files:           
            data_path_name_split = data_files.split('_')    
            #integration_type = [x for x in names_to_check if x in data_path_name_split]
            #collection_type = integration_type[0]
            #frames_int = data_path_name_split[-1]            
            int_time = round(int(data_path_name_split[-1].split('.')[0])) #for integ_sweep             
            data_file = os.path.join(file_path, saved_quads, data_files)  
            IDL_variable = readsav(data_file)   
           # print(data_file)
            #cc
            all_full_frame = IDL_variable.q  
            all_int_time.append(int_time)
            
            quad_A = all_full_frame[:, 0, :, :]
            active_quad_A = quad_A[:, 4:1028, 10:1034]
            quad_A_active = copy.deepcopy(active_quad_A)
            bias_A = quad_A[:, 4:1028, 1034:1056]
            bias_subtracted_quad_A = perform_bias_subtraction (active_quad_A, bias_A)      
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_A)
            var_quad_A = np.var(outlier_filtered_data)     
            median_quad_A = np.median(outlier_filtered_data)            
            all_med_quad_A.append(median_quad_A)
            all_var_quad_A.append(var_quad_A)
           
            
            
            # For quad B
            quad_B = all_full_frame[:, 1, :, :]
            active_quad_B = quad_B[:, 4:1028, 10:1034]
            bias_B = quad_B[:, 4:1028, 1034:1056]
            bias_subtracted_quad_B = perform_bias_subtraction (active_quad_B, bias_B)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_B) 
            var_quad_B = np.var(outlier_filtered_data)                 
            median_quad_B = np.median(bias_subtracted_quad_B)           
            all_med_quad_B.append(median_quad_B)
            all_var_quad_B.append(var_quad_B)
            
            
            # For quad C
            quad_C = all_full_frame[:,2,:,:]
            active_quad_C = quad_C[:, 4:1028, 10:1034]
            bias_C = quad_C[:, 4:1028, 1034:1056]
            bias_subtracted_quad_C = perform_bias_subtraction (active_quad_C, bias_C)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_C) 
            var_quad_C = np.var(outlier_filtered_data)                 
            median_quad_C = np.median(bias_subtracted_quad_C)           
            all_med_quad_C.append(median_quad_C)
            all_var_quad_C.append(var_quad_C)
            
            
            quad_D = all_full_frame[:,3,:,:]
            active_quad_D = quad_D[:, 4:1028, 10:1034]
            bias_D = quad_D[:, 4:1028, 1034:1056]
            bias_subtracted_quad_D = perform_bias_subtraction (active_quad_D, bias_D)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_D) 
            var_quad_D = np.var(outlier_filtered_data)                 
            median_quad_D = np.median(bias_subtracted_quad_D)           
            all_med_quad_D.append(median_quad_D)
            all_var_quad_D.append(var_quad_D)
            
             #let's save quad images too
             
#            quad_a = np.mean(quad_A, axis=0)            
#            quad_b = np.fliplr(np.mean(quad_B, axis=0))
#            quad_c = np.rot90(np.mean(quad_C, axis=0), 2)
#            #print(quad_c.shape)
#            quad_d = np.flipud(np.mean(quad_D, axis=0))            
#            lower_quads = np.concatenate((quad_a, quad_b), axis=1)
#            upper_quads = np.concatenate((quad_d, quad_c), axis=1)
#            quads = np.concatenate((lower_quads, upper_quads), axis=0)
#            plt.figure()
#            image = plt.imshow(quads, cmap='nipy_spectral', origin='lower', interpolation='none')
#            plt.grid(b=False)
#            cbar = plt.colorbar(image)
#            title = 'Full Frame Quad Image' + ' (Int.time = '+ str(int_time)+')' 
#            plt.title(title, fontsize=14)
#            plt.xlabel('# of spatial pixels', fontsize=12)
#            plt.ylabel('# of spectral pixels', fontsize=12)
#            plt.grid(False)   
#            plt.savefig(plot_dir+'/'+ data_files+'.png',dpi=100,bbox_inches="tight")
#            plt.close('all')
       
        print('A:', all_med_quad_A)
        print('B:', all_med_quad_B)
        print('C:',  all_med_quad_C)
        print('D:', all_med_quad_D)
         
        dframe1 = pd.DataFrame(
                    {'Int_time.' : all_int_time,
                     'Avg. Quad A' : all_med_quad_A,
                     'Avg. Quad B' : all_med_quad_B,
                     'Avg. Quad C' : all_med_quad_C,
                     'Avg. Quad D' : all_med_quad_D,
                     'Avg. Var_Quad_A': all_var_quad_A,
                     'Avg. Var_Quad_B': all_var_quad_B,
                     'Avg. Var_Quad_C': all_var_quad_C,
                     'Avg. Var_Quad_D': all_var_quad_D,
                     
                     })
                
        dframe1.to_csv(plot_dir +'/'+'Median_DN_vs_Var.csv')
        #dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test\Unct_Quad_A_flooded.csv')
        
if __name__ == "__main__":
    main()