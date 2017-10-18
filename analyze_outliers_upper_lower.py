# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

The purpose of this code is to identify outliers on TEMPO CCDs. To begin with,
it takes the raw focal plane plane data in binay format.
This bascailly takes Dave F. IDL variables and does the analysis in python
The hardcoded variables are standard to tempo instruments.
The tests are based on both the photon transfer data as well as well as dark
current data

"""

import os
import numpy as np
import scipy.io
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
#*****************************************************************************
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
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads'
    data_path_all = [each for each in os.listdir(file_path) if each.endswith('.sav')]
    #names_to_check = ['SHORT', 'NOM', 'LONG']
    names_to_check = ['SHORT', 'NOM','OP', 'LONG']

    mean_mask_all_lower = []
    median_mask_all_lower = []    
    mean_mask_all_upper = []
    median_mask_all_upper = []
    saturation_mask_upper = []
    saturation_mask_lower = []
    for data_path in data_path_all:

    
    # Now let's find the integtration type
        print(data_path)
        full_frame, collection_type, frames, int_time = \
        make_quads_from_IDL_sav_file(file_path, data_path, names_to_check)
        #Let's fetch the quads now
        #quad_A = full_frame[0:ny_quad, 0:nx_quad]
        active_quad_A = np.array(full_frame[4:1028, 10:1034])
        #quad_B = full_frame[ny_quad:, 0:nx_quad]
        active_quad_B = np.array(full_frame[4:1028, 1078:2102])
        #quad_C = full_frame[ny_quad:, nx_quad:]
        active_quad_C = np.array(full_frame[1064:2088, 1078:2102])
        #quad_D = full_frame[ny_quad:, 0:nx_quad]
        active_quad_D = np.array(full_frame[1064:2088, 10:1034])
        lower_quads = np.concatenate((active_quad_A, active_quad_B), axis=1)
        upper_quads = np.concatenate((active_quad_D, active_quad_C), axis=1)
        active_quads =[lower_quads, upper_quads]
       
        
        #quads = ['Lower_quads','Upper_quads']
        
        # Now let's save the full_frame images
        plot_dir = file_path
        image_dir = 'Quad_images'
        save_quad_images = os.path.join(plot_dir, image_dir)
        
        if not os.path.exists(save_quad_images):
            os.makedirs(save_quad_images)
        image_title = 'Full Frame Quad Image'
        plot_full_frame_image(full_frame, image_title, collection_type, frames,
                              int_time, save_quad_images)
        # Create list of directories To save histograms before and after
        #outlier rejections
        
        
        mean_plot = r'Outlier_mean'+'/'+ 'Lower_quads'
        median_plot = r'Outlier_median'+'/'+ 'Lower_quads'
        folder_name_hist = 'Hist_plot'
        folder_name_mask = 'Mask_plot'
        folder_name_sat = 'Saturation_plot'
        
        save_mean_hist = os.path.join(plot_dir, mean_plot, folder_name_hist)
        if not os.path.exists(save_mean_hist):
            os.makedirs(save_mean_hist)

        save_median_hist = os.path.join(plot_dir, median_plot, folder_name_hist)
        if not os.path.exists(save_median_hist):
            os.makedirs(save_median_hist)

        save_mean_mask = os.path.join(plot_dir, mean_plot, folder_name_mask)
        if not os.path.exists(save_mean_mask):
            os.makedirs(save_mean_mask)

        save_median_mask = os.path.join(plot_dir, median_plot, folder_name_mask)
        if not os.path.exists(save_median_mask):
            os.makedirs(save_median_mask)
         
        save_sat_mask = os.path.join(plot_dir, median_plot, folder_name_sat)
        if not os.path.exists(save_sat_mask):
            os.makedirs(save_sat_mask)

        # Now let's call the outlier functions. First identify saturated
        #bits and then identify the 'hot and cold' outlier pixels
        
       
        title = 'Histogram of Quads A & B'
     
            #title = 'Histogram of Quads C & D'
        
        sat_pixels = identify_saturation(lower_quads,
                                                     save_sat_mask,
                                                     collection_type,
                                                     frames, int_time)            
        
        #outlier_filt_mean, outlier_mean = reject_outlier_mean(active_quads[k])

        outlier_filt_med, outlier_med = reject_outlier_median(lower_quads)

           

        # Ok now lets' plot the histograms to understand what we see.
        # First the mean approach
        print('ok, let us plot the histograms now')
        
        #plot_hist_image(active_quads[k], title, collection_type, frames,
                        #outlier_filt_mean, outlier_mean, int_time,
                       # save_mean_hist)
        
        
        
        # Secondly the median absolute difference approach
        plot_hist_image(lower_quads, title, collection_type, frames,
                        outlier_filt_med, outlier_med, int_time,
                        save_median_hist)
        #outlier_med= None
        # Ok, lets now create a binary mask of outlier pixels
        # First with the mean based approach


#            mean_mask = create_outlier_mask(active_quads[k], outlier_mean[0],
#                                            collection_type,
#                                            frames, save_mean_mask)
#            if k == 0:
#                mean_mask_all_lower.append(mean_mask)
#                
#            else :
#                mean_mask_all_upper.append(mean_mask)
#    
#            # Secondly, the same for median ased approach
      
        #if len(outlier_med) == 1:
         #   outlier_med = 0
        #else:
            #outlier_med = outlier_med[1]
        median_mask = create_outlier_mask(lower_quads, outlier_med[1],
                                          collection_type, frames,
                                          save_median_mask)
        lower_quads = None
        outlier_med = None

        median_mask_all_lower.append(median_mask[:])
        median_mask=None
        saturation_mask_lower.append(sat_pixels[:])
        
        
    
    # Ok now, lets' or all the median masks
    
    print(np.array(median_mask_all_lower).shape)    
    final_mask_median_lower = np.bitwise_or.reduce(median_mask_all_lower[:])

    plt.figure()
    plt.imshow(np.invert(final_mask_median_lower), cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title('Binary outlier mask', fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(r'C:\Users\nmishra\Desktop'+'/1.png', dpi=100, bbox_inches="tight")       
    cc
    
    plt.figure()
    plt.imshow(np.invert(median_mask_all_lower[1]), cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title('Binary outlier mask', fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(r'C:\Users\nmishra\Desktop'+'/2.png', dpi=100, bbox_inches="tight")          
            
                
    
     

if __name__ == "__main__":
    main()
