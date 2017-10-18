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
                                    make_quads_from_fits_file,\
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
    Tme main function    """
   
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark_Current_Radiation_testing\Pre_rad'
    data_path_all = [each for each in os.listdir(file_path) if each.endswith('.fits')]
    collection_type = 'Pred rad Test'

    median_mask_all_quad_A = []
    median_mask_all_quad_B = []     
    for data_path in data_path_all: 
        
        full_frame, int_time, temp, frames =  make_quads_from_fits_file(file_path, data_path)
        # Radiation testing was done only to 1 CCD ( quads A and B)        
        #There are 16 frames each for Quads A and B. CLoser inspection suggested
        #that you can average these 16 frames and work on average
        # Let's fetch the quads now
        ndimes, nx_quad, ny_quad = full_frame.shape
        quad_A = np.median(full_frame[0:ndimes:2, :, :], axis=0)
       
        
        active_quad_A = np.array(quad_A[4:1028, 10:1034])
        
        quad_B = np.median(full_frame[1:ndimes:2, :, :], axis=0)        
        active_quad_B = np.array(quad_B[4:1028, 10:1034])               
        quads_name = ['Quad A','Quad B'] 
        lower_quads = [active_quad_A, active_quad_B]
        lower_quads_all = [quad_A, quad_B]
        plot_dir = file_path
        
       
        # Create list of directories To save histograms before and after
        #outlier rejections
        
        for k in np.arange(0,2):
            image_dir = r'Quad_images' +'/'+ quads_name[k]            
            median_plot = r'Outlier_median'+'/'+ quads_name[k]
            folder_name_hist = 'Hist_plot'
            folder_name_mask = 'Mask_plot'
            folder_name_sat = 'Saturation_plot'
            
            save_quad_images = os.path.join(plot_dir, image_dir)        
            if not os.path.exists(save_quad_images):
                os.makedirs(save_quad_images)
            
            
    
            save_median_hist = os.path.join(plot_dir, median_plot, folder_name_hist)
            if not os.path.exists(save_median_hist):
                os.makedirs(save_median_hist)
    
            
            save_median_mask = os.path.join(plot_dir, median_plot, folder_name_mask)
            if not os.path.exists(save_median_mask):
                os.makedirs(save_median_mask)
             
            save_sat_mask = os.path.join(plot_dir, median_plot, folder_name_sat)
            if not os.path.exists(save_sat_mask):
                os.makedirs(save_sat_mask)
    
            
            # Lets plt the frame first
            if k == 0:
                image_title = 'Full Frame Quad A Image'
            else:
                image_title = 'Full Frame Quad B Image'
            
            int_time = str(int_time)
   
            title = image_title + ' ('+ collection_type+')' + '\n Int. time = ' + int_time + \
                                    ', Temp = ' +temp   
                        
            plt.figure()
            image = plt.imshow(lower_quads_all[k], cmap='nipy_spectral', interpolation='none', 
                       origin='lower')
            cbar = plt.colorbar(image)
            plt.title(title, fontsize=14)
            plt.xlabel('# of spatial pixels', fontsize=12)
            plt.ylabel('# of spectral pixels', fontsize=12)
            plt.grid(False)
    
            plt.savefig(save_quad_images+'/'+ quads_name[k]+'_'+\
                    frames+'.png',dpi=100,bbox_inches="tight")
            plt.close('all')
            
            
            
            # Now let's call the outlier functions. First identify saturated
            #bits and then identify the 'hot and cold' outlier pixels
            
            if k == 0:
                title = 'Histogram of Quad A'
            else:
                title = 'Histogram of Quad B'
            
#            sat_pixels = identify_saturation(lower_quads[k],
#                                                         save_sat_mask,
#                                                         collection_type,
#                                                         frames, int_time)            
            
            outlier_filt_med, outlier_med = reject_outlier_median(lower_quads[k])
          # Ok now lets' plot the histograms to understand what we see.
            # First the mean approach
            print('ok, let us plot the histograms now')           
            # Secondly the median absolute difference approach
            plot_hist_image(lower_quads[k], title, collection_type, frames,
                            outlier_filt_med, outlier_med, int_time,
                            save_median_hist)
                      
            if len(outlier_med) == 1:
                outlier_med = 0
            else:
                outlier_med = outlier_med[1]
             
                
            if k == 0:
                title = 'Binary Outlier Mask of Quad A'
            else:
                title = 'Binary Outlier Mask of Quad B'   
                
            median_mask = create_outlier_mask(lower_quads[k], outlier_med,
                                              title, collection_type, frames,
                                              save_median_mask)
            lower_quads[k] = None
            outlier_med = None
    
            if k == 0:
                median_mask_all_quad_A.append(median_mask[:])
                median_mask=None
                #saturation_mask_lower.append(sat_pixels[:])
                #sat_pixels = None
            elif k==1 :
               median_mask_all_quad_B.append(median_mask[:])
               median_mask = None
                #saturation_mask_upper.append(sat_pixels[:])
                #sat_pixels = None
                    
                
   #*************************************************************************** 
    # The following function creates final mask based on 80% occurence criteria
    
    final_path = file_path+'\Outlier_median\save_masks'
    if not os.path.exists(final_path):
                os.makedirs(final_path)
    quad_name = 'QuadA_80%'
    title = 'Quad A outlier mask (80% occurence criteria)'
    final_mask_quad_A = create_final_mask(median_mask_all_quad_A, quad_name, title, final_path)
    quad_name = 'QuadB_80%'
    title = 'Quad B outlier mask (80% ocurence criteria)'
    final_mask_quad_B = create_final_mask(median_mask_all_quad_B, quad_name, title, final_path)
    final_outlier_mask = np.concatenate((final_mask_quad_A , final_mask_quad_B ), 
                                         axis=1)
    plt.figure()
    plt.imshow(final_outlier_mask, cmap='bwr', 
               interpolation='none', origin='lower')
    
    nx_quad, ny_quad = np.array(final_outlier_mask).shape
    final_outliers = np.array(np.where(np.reshape(final_outlier_mask,
                                          (nx_quad*ny_quad, 1))==0)).shape[1]     
    plt.title('Quad A and B Outlier Mask'+' (outliers = '+ str(final_outliers)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    
    
    plt.savefig(final_path+'/'+'final_mask.png', dpi=100, bbox_inches="tight")
    
   #*************************************************************************** 

   #*************************************************************************** 
    # if we want to 'OR" all the mask, use the following function
    
    quad_name = 'QuadA_ored'  
    title = 'Quad A outlier mask (Logical OR)'
    final_mask_lower_ored= create_ORed_mask(median_mask_all_quad_A, quad_name, title, final_path)
    quad_name = 'QuadB_ored'
    title = 'Quad B outlier mask (Logical OR)'
    final_mask_upper_ored = create_ORed_mask(median_mask_all_quad_B, quad_name, title, final_path)
    final_outlier_mask_ored = np.concatenate((final_mask_lower_ored,
                                              final_mask_upper_ored), axis=1)
    plt.figure()
    plt.imshow(np.invert(final_outlier_mask_ored), cmap='bwr', 
               interpolation='none', origin='lower')
    
    nx_quad, ny_quad = np.array(final_outlier_mask_ored).shape
    final_outliers_ored = np.array(np.where(np.reshape(final_outlier_mask_ored,
                                          (nx_quad*ny_quad, 1))==1)).shape[1]     
    plt.title('Quad A & B Outlier Mask'+ ' (outliers = '+ str(final_outliers_ored)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)    
    
    plt.savefig(final_path+'/'+'final_mask_ored.png', dpi=100, bbox_inches="tight")
    #**************************************************************************
    
    
    # TO DO : Add validation work 
if __name__ == "__main__":
    main()
