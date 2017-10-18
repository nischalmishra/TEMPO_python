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
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def create_image(image_data, title, figure_name):
    plt.figure()  
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='nipy_spectral', origin='lower')   
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    ax.grid(False)
    plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    #plt.show()    
    plt.close('all')

def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.06.28\saved_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask_spectrometer'
    all_int_files = [each for each in os.listdir(file_path) if each.endswith('.sav')]
    #data_path_all = [items for items in all_int_files if not items.endswith('_INT_0.dat.sav')]

    #names_to_check = ['SHORT', 'NOM', 'LONG']
    #names_to_check = ['SHORT', 'NOM','OP', 'LONG']
    median_mask_A = []
    median_mask_B = []    
    median_mask_C = []
    median_mask_D = []
    saturation_mask_A = []
    saturation_mask_B = []
    saturation_mask_C = []
    saturation_mask_D = []
    for data_path in all_int_files:
        print(data_path)
        data_file = os.path.join(file_path, data_path)
        IDL_variable = readsav(data_file)
        full_frame = IDL_variable.quads       
        active_quad_A = full_frame[0, 4:1028, 10:1034]
        active_quad_B = full_frame[1, 4:1028, 10:1034]
        active_quad_C = full_frame[2, 4:1028, 10:1034]
        active_quad_D = full_frame[3, 4:1028, 10:1034]      
        active_quads = [active_quad_A, active_quad_B, active_quad_C, active_quad_D]             
        quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']        
        # Now let's save the full_frame images
        int_time = 139.97
        collection_type = 'Nom Int.'
        frames = data_path.split('.')[0]
        
        for k in np.arange(0, 4):
            #mean_plot = r'Outlier_mean_tsoc'+'/'+ quads[k]
            
            median_plot = r'Outlier_median'+'/'+ quads[k]
            image_dir = 'Saved_Images'
            folder_name_hist = 'Hist_plot'
            folder_name_mask = 'Mask_plot'
            folder_name_sat = 'Saturation_plot'  
       
            save_median_hist = os.path.join(save_dir, median_plot, folder_name_hist)
            if not os.path.exists(save_median_hist):
                os.makedirs(save_median_hist)
                
            saved_images = os.path.join(save_dir, median_plot, image_dir)            
            if not os.path.exists(saved_images):
                os.makedirs(saved_images)
            title = quads[0]+'image (Spectrometer)'
            figure_name = saved_images+'/'+frames+'.png'
            create_image(active_quads[k], title, figure_name)
                
            save_median_mask = os.path.join(save_dir, median_plot, folder_name_mask)
            if not os.path.exists(save_median_mask):
                os.makedirs(save_median_mask)
             
            save_sat_mask = os.path.join(save_dir, median_plot, folder_name_sat)
            if not os.path.exists(save_sat_mask):
                os.makedirs(save_sat_mask)
            
            # Now let's call the outlier functions. First identify saturated
            #bits and then identify the 'hot and cold' outlier pixels
            
            title = 'Histogram of '+ quads[k]            
            
            sat_pixels = identify_saturation(active_quads[k],
                                                         save_sat_mask,
                                                         collection_type,
                                                         frames, int_time)            
            
            #outlier_filt_mean, outlier_mean = reject_outlier_mean(active_quads[k])
    
            outlier_filt_med, outlier_med = reject_outlier_median(active_quads[k])
    
            # Ok now lets' plot the histograms to understand what we see.
            # First the mean approach
            print('ok, let us plot the histograms now')           
            
            plot_hist_image(active_quads[k], title, collection_type, frames,
                            outlier_filt_med, outlier_med, int_time,
                            save_median_hist)
           
            if len(outlier_med) == 1:
                outlier_med = 0
            else:
                outlier_med = outlier_med[1] 
            title = 'Binary Outlier Mask of ' + quads[k]            
            median_mask = create_outlier_mask(active_quads[k], outlier_med,
                                              title, collection_type, frames,
                                              save_median_mask)
            active_quads[k] = None
            outlier_med = None      
            if k == 0:
                median_mask_A.append(median_mask[:])
                median_mask=None
                saturation_mask_A.append(sat_pixels[:])
                sat_pixels = None
            elif k==1 :
                median_mask_B.append(median_mask[:])
                median_mask=None
                saturation_mask_B.append(sat_pixels[:])
                sat_pixels = None
                   
            elif k==2 :
                median_mask_C.append(median_mask[:])
                median_mask=None
                saturation_mask_C.append(sat_pixels[:])
                sat_pixels = None 
                
            elif k==3 :
                median_mask_D.append(median_mask[:])
                median_mask=None
                saturation_mask_D.append(sat_pixels[:])
                sat_pixels = None
   #*************************************************************************** 
    # The following function creates final mask based on 80% occurence criteria
    
    final_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask_spectrometer\Outlier_median\save_masks_80%'
    if not os.path.exists(final_path):
                os.makedirs(final_path)
    quad_name = 'Quad A'
    title = 'Quad A outlier mask (50% criteria)'
    final_mask_A = create_final_mask(median_mask_A, quad_name, title, final_path)
    mask_name_A = final_path+'/'+'quad_A_outlier_mask.csv'
    np.savetxt(mask_name_A, np.array(final_mask_A), delimiter =",")
    
    quad_name = 'Quad B'
    title = 'Quad B outlier mask (50% criteria)'                      
    final_mask_B = create_final_mask(median_mask_B, quad_name, title, final_path)
    mask_name_B = final_path+'/'+'quad_B_outlier_mask.csv'
    np.savetxt(mask_name_B, np.array(final_mask_B), delimiter =",")
    
    quad_name = 'Quad C'
    title = 'Quad C outlier mask (50% criteria)'
    final_mask_C = create_final_mask(median_mask_C, quad_name, title, final_path)
    mask_name_C = final_path+'/'+'quad_C_outlier_mask.csv'
    np.savetxt(mask_name_C, np.array(final_mask_C), delimiter =",")
    
    quad_name = 'Quad D'
    title = 'Quad D outlier mask (50% criteria)'                      
    final_mask_D = create_final_mask(median_mask_D, quad_name, title, final_path)
    mask_name_D = final_path+'/'+'quad_D_outlier_mask.csv'
    np.savetxt(mask_name_D, np.array(final_mask_D), delimiter =",")
    
    lower_quads = np.concatenate((final_mask_A, final_mask_B), axis=1)
    upper_quads = np.concatenate((final_mask_D, final_mask_C), axis=1)    
    final_outlier_mask = np.concatenate((lower_quads, upper_quads), axis=0)
    
    mask_name = final_path+'/'+'final_outlier_mask_80%.csv'	
    np.savetxt(mask_name, np.array(final_outlier_mask), delimiter =",")
    
    plt.figure()
    final_outlier_mask = np.invert(final_outlier_mask)
    plt.imshow(final_outlier_mask, cmap='bwr', 
               interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    nx_quad, ny_quad = np.array(final_outlier_mask).shape
    final_outliers = np.array(np.where(np.reshape(final_outlier_mask,
                                          (nx_quad*ny_quad, 1))==0)).shape[1]     
    plt.title('TEMPO outlier Mask (50% Occurence Criteria)\n (Num. outliers = '+ str(final_outliers)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(final_path+'/'+'final_mask_50%.png', dpi=200, bbox_inches="tight")
    
   #*************************************************************************** 

   #*************************************************************************** 
    # if we want to 'OR" all the mask, use the following function
    
     
    final_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask_spectrometer\Outlier_median\save_masks_logical_OR'
    if not os.path.exists(final_path):
                os.makedirs(final_path)
    quad_name = 'Quad A'  
    title = 'Quads A outlier mask (Logical OR)'
    final_mask_A= create_ORed_mask(median_mask_A, quad_name, title, final_path)
    mask_name_A = final_path+'/'+'quad_A_outlier_mask_ored.csv'
    np.savetxt(mask_name_A, np.array(final_mask_A), delimiter =",")

    quad_name = 'Quad B'  
    title = 'Quads B outlier mask (Logical OR)'
    final_mask_B = create_ORed_mask(median_mask_B, quad_name, title, final_path)
    mask_name_B = final_path+'/'+'quad_B_outlier_mask_ored.csv'
    np.savetxt(mask_name_B, np.array(final_mask_B), delimiter =",")
    
    
    quad_name = 'Quad C'  
    title = 'Quads C outlier mask (Logical OR)'
    final_mask_C= create_ORed_mask(median_mask_C, quad_name, title, final_path)
    mask_name_C = final_path+'/'+'quad_C_outlier_mask_ored.csv'
    np.savetxt(mask_name_C, np.array(final_mask_C), delimiter =",")
    
    quad_name = 'Quad D'  
    title = 'Quads D outlier mask (Logical OR)'
    final_mask_D = create_ORed_mask(median_mask_D, quad_name, title, final_path)
    mask_name_D = final_path+'/'+'quad_D_outlier_mask_ored.csv'
    np.savetxt(mask_name_D, np.array(final_mask_D), delimiter =",")
    
    quad_name = 'upper_quads_ored'
    title = 'Quads C & D outlier mask (Logical OR)'
    
    lower_quads = np.concatenate((final_mask_A, final_mask_B), axis=1)
    upper_quads = np.concatenate((final_mask_D, final_mask_C), axis=1)    
    final_outlier_mask_ored = np.concatenate((lower_quads, upper_quads), axis=0)                       
    mask_name = final_path+'/'+'final_outlier_mask_ored.csv'	
    np.savetxt(mask_name, np.array(final_outlier_mask_ored), delimiter =",")
    
    plt.figure()
    plt.imshow(np.invert(final_outlier_mask_ored), cmap='bwr', 
               interpolation='none', origin='lower')
    nx_quad, ny_quad = np.array(final_outlier_mask_ored).shape
    final_outliers_ored = np.array(np.where(np.reshape(final_outlier_mask_ored,
                                          (nx_quad*ny_quad, 1))==1)).shape[1]     
    plt.title('TEMPO Outlier Mask (Logical OR)\n (Num. outliers = '+ str(final_outliers_ored)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(final_path+'/'+'final_mask_ored.png', dpi=200, bbox_inches="tight")
    #**************************************************************************
    
    
    # TO DO : Add validation work 
if __name__ == "__main__":
    main()
