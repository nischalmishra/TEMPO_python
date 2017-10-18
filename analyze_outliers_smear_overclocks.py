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
                              create_final_mask,\
                              create_ORed_mask
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
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads'
    data_path_all = [each for each in os.listdir(file_path) if each.endswith('.sav')]
    #names_to_check = ['SHORT', 'NOM', 'LONG']
    names_to_check = ['SHORT', 'NOM','OP', 'LONG']

    
    median_mask_all_A = []
    median_mask_all_B = []
    median_mask_all_C = []
    median_mask_all_D = []
    saturation_mask_A = []
    saturation_mask_B = []
    saturation_mask_C = []
    saturation_mask_D = []
    for data_path in data_path_all:
    
    
        print(data_path)
        full_frame, collection_type, frames,int_time = \
        make_quads_from_IDL_sav_file(file_path, data_path, names_to_check)
        # all_full_frame contains all 100 frames
        
        #Let's fetch the quads now
        #quad_A = full_frame[0:ny_quad, 0:nx_quad]
        #active_quad_A = np.array(full_frame[4:1028, 10:1034])
        smear_A = np.array(full_frame[ 1028:1046, 10:1034])
        #quad_B = full_frame[ny_quad:, 0:nx_quad]
        #active_quad_B = np.array(full_frame[4:1028, 1078:2102])
        smear_B = np.array(full_frame[ 1028:1046, 1078:2102])
        #quad_C = full_frame[ny_quad:, nx_quad:]
        #active_quad_C = np.array(full_frame[1064:2088, 1078:2102])
        smear_D = np.array(full_frame[1046:1064, 10:1034])
        #quad_D = full_frame[ny_quad:, 0:nx_quad]
        #active_quad_D = np.array(full_frame[1064:2088, 10:1034])
        smear_C = np.array(full_frame[ 1046:1064, 1078:2102])
        #lower_quads = np.concatenate((active_quad_A, active_quad_B), axis=1)
        #upper_quads = np.concatenate((active_quad_D, active_quad_C), axis=1)
        #active_quads =[lower_quads, upper_quads]   
        # for overclocks
        
        active_quads =[smear_A, smear_B, smear_C, smear_D]   
            
        quads = ['Quad A','Quad B', 'Quad C', 'Quad D']
        
        # Now let's save the full_frame images
        plot_dir = file_path+'/Smear_Overclocks'
        if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
  
        # Create list of directories To save histograms before and after
        #outlier rejections
        
        for k in np.arange(0, 4):
            #mean_plot = r'Outlier_mean_tsoc'+'/'+ quads[k]
            median_plot = r'Outlier_median_smear'+'/'+ quads[k]
            folder_name_hist = 'Hist_plot'
            folder_name_mask = 'Mask_plot'
            folder_name_sat = 'Saturation_plot'
            
            
    
            save_median_hist = os.path.join(plot_dir, median_plot, folder_name_hist)
            if not os.path.exists(save_median_hist):
                os.makedirs(save_median_hist)
    
            
    
            save_median_mask = os.path.join(plot_dir, median_plot, folder_name_mask)
            if not os.path.exists(save_median_mask):
                os.makedirs(save_median_mask)
             
            save_sat_mask = os.path.join(plot_dir, median_plot, folder_name_sat)
            if not os.path.exists(save_sat_mask):
                os.makedirs(save_sat_mask)
            
            # Now let's call the outlier functions. First identify saturated
            #bits and then identify the 'hot and cold' outlier pixels
            
            if k == 0:
                title = 'Histogram of Smear Overclocks of Quad A'
            elif k==1:
                title = 'Histogram of Smear Overclocks of Quad B'
            elif k==2:
                title = 'Histogram of Smear Overclocks of Quad C'
            else:
                title = 'Histogram of Smear Overclocks of Quad C'
            
            sat_pixels = identify_saturation(active_quads[k],
                                                         save_sat_mask,
                                                         collection_type,
                                                         frames, int_time)            
            
            #outlier_filt_mean, outlier_mean = reject_outlier_mean(active_quads[k])
    
            outlier_filt_med, outlier_med = reject_outlier_median(active_quads[k])
    
               
    
            # Ok now lets' plot the histograms to understand what we see.
            # First the mean approach
            print('ok, let us plot the histograms now')
            
            #plot_hist_image(active_quads[k], title, collection_type, frames,
                            #outlier_filt_mean, outlier_mean, int_time,
                           # save_mean_hist)
            
            
            
            # Secondly the median absolute difference approach
            plot_hist_image(active_quads[k], title, collection_type, frames,
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
          
            if len(outlier_med) == 1:
                outlier_med = 0
            else:
                outlier_med = outlier_med[1]
                
            if k == 0:
                title = 'Binary Outlier Mask of Quad A Smear overclocks'
            elif k==1:
                title = 'Binary Outlier Mask of Quad B Smear overclocks'
            elif k == 2:
                title = 'Binary Outlier Mask of Quad C Smear overclocks'
            else:
                title = 'Binary Outlier Mask of Quad D Smear overclocks'    
                
                
            median_mask = create_outlier_mask(active_quads[k], outlier_med,
                                              title, collection_type, frames,
                                              save_median_mask)
            
            
            active_quads[k] = None
            outlier_med = None
                        
                
            if k == 0:
                median_mask_all_A.append(median_mask[:])
                median_mask=None
                saturation_mask_A.append(sat_pixels[:])
                sat_pixels = None
            elif k==1 :
                median_mask_all_B.append(median_mask[:])
                median_mask = None
                saturation_mask_B.append(sat_pixels[:])
                sat_pixels = None
                
            elif k==2 :
                median_mask_all_C.append(median_mask[:])
                median_mask = None
                saturation_mask_C.append(sat_pixels[:])
                sat_pixels = None      
                   
            elif k==3 :
                median_mask_all_D.append(median_mask[:])
                median_mask = None
                saturation_mask_D.append(sat_pixels[:])
                sat_pixels = None      
       
   #*************************************************************************** 
    # The following function creates final mask based on 80% occurence criteria
    
    cc
    final_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads\Outlier_median\save_masks'
    if not os.path.exists(final_path):
                os.makedirs(final_path)
    quad_name = 'lower_quads'
    title = 'Quads A & B outlier mask (80% criteria)'
    final_mask_lower = create_final_mask(median_mask_all_lower, quad_name, title, final_path)
    quad_name = 'upper_quads'
    title = 'Quads C & D outlier mask (80% criteria)'
    final_mask_upper = create_final_mask(median_mask_all_upper, quad_name, title, final_path)
    final_outlier_mask = np.concatenate((final_mask_lower, final_mask_upper), 
                                         axis=0)
    plt.figure()
    plt.imshow(final_outlier_mask, cmap='bwr', 
               interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    nx_quad, ny_quad = np.array(final_outlier_mask).shape
    final_outliers = np.array(np.where(np.reshape(final_outlier_mask,
                                          (nx_quad*ny_quad, 1))==0)).shape[1]     
    plt.title('TEMPO outlier Mask (outliers = '+ str(final_outliers)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    
    final_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads\Outlier_median\save_masks'
    if not os.path.exists(final_path):
                os.makedirs(final_path)
    plt.savefig(final_path+'/'+'final_mask.png', dpi=100, bbox_inches="tight")
   #*************************************************************************** 

   #*************************************************************************** 
    # if we want to 'OR" all the mask, use the following function
    quad_name = 'lower_quads_or'  
    title = 'Quads A & B outlier mask (Logical OR)'
    final_mask_lower_ored= create_ORed_mask(median_mask_all_lower, quad_name, title, final_path)
    quad_name = 'upper_quads_ored'
    title = 'Quads C & D outlier mask (Logical OR)'
    final_mask_upper_ored = create_ORed_mask(median_mask_all_lower, quad_name, title, final_path)
    final_outlier_mask_ored = np.concatenate((final_mask_lower_ored,
                                              final_mask_upper_ored), axis=0)
    plt.figure()
    plt.imshow(np.invert(final_outlier_mask_ored), cmap='bwr', 
               interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    nx_quad, ny_quad = np.array(final_outlier_mask_ored).shape
    final_outliers_ored = np.array(np.where(np.reshape(final_outlier_mask_ored,
                                          (nx_quad*ny_quad, 1))==1)).shape[1]     
    plt.title('TEMPO Outlier Mask (outliers = '+ str(final_outliers_ored)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)    
    
    plt.savefig(final_path+'/'+'final_mask_ored.png', dpi=100, bbox_inches="tight")
    #**************************************************************************
    
    
    # TO DO : Add validation work 
if __name__ == "__main__":
    main()
