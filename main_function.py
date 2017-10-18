# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

This code was written to process the tempo raw focal plane electronics (FPE)
Data. The hardcoded variables are standard to tempo instruments. Basically it's
taken from Dave Flittner's IDL script and pythonized. The sanity checks were
done by comparing the outputs from both the script.

"""

import os
import numpy as np

import matplotlib.pyplot as plt
from make_quads_from_raw_fpe import make_full_frame_from_raw_fpe,\
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

from outlier_detection import reject_outlier_median,\
                              reject_outlier_mean,\
                              create_outlier_mask ,\
                              perform_logical_and
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
    nx_quad = 1056 # For Tempo
    ny_quad = 1046 # For Tempo
    nlat = nx_quad*2
    nspec = ny_quad*2
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Background'
    data_path = [each for each in os.listdir(file_path) if each.endswith('.dat')]
    names_to_check = ['SHORT', 'NOM', 'LONG']
    
    intermediate_mask_mean = []
    intermediate_mask_median = []
    for data_path in data_path:      
    
    #data_path = 'SPI_20160912142518420_LIGHT_INTEGRATION_SWEEP_207_FT6_NOM_INT_99999.dat'
    # Now let's find the integtration type
    
        data_path_name_split = data_path.split('_')
        print(data_path_name_split)
        integration_type = [x for x in names_to_check if x in data_path_name_split]
        frames = data_path_name_split[-1]
        data_file = os.path.join(file_path, data_path)
        fileinfo = get_size(data_file)
        size_file = fileinfo.st_size
        #-------------Some printing to do Sanity checks----------------------------
        #print ("size of the file is ", size_file, "bytes")
        nframes = size_file/(nlat*nspec*2)
        #print("num of frames : ", nframes)
        #dt = np.dtype(')
        result = np.fromfile(data_file, dtype='>u2', count=-1) & int('3fff', 16)
        # Refer to 2450865 Ball document
        #print(result.size)
        #*Ok, let's do few steps before we move into data analysis
    
        # 1. Create full frame from raw FPE.
        # 2. Arrange full frame into various quads.
        full_frame = make_full_frame_from_raw_fpe(result, nx_quad, ny_quad)
        quads = full_frame_to_quads(full_frame, nx_quad, ny_quad)
        
        # Now let's save the full_frame images
        plot_dir = r'C:\Users\nmishra\Desktop\test'
        image_dir = 'Quad_images'
        save_quad_images= os.path.join(plot_dir,image_dir)
        if not os.path.exists(save_quad_images):os.makedirs(save_quad_images)        
        
        aligned_quads = 'aligned_quads'
        save_aligned_quads = os.path.join(plot_dir,aligned_quads)
        if not os.path.exists(save_aligned_quads):os.makedirs(save_aligned_quads)
        
        collection_type = integration_type[0]
        image_title = 'Full Frame Quad Image'
        plot_full_frame_image(full_frame, image_title, collection_type, frames,
                              save_quad_images)
        plot_each_quad(quads, image_title, collection_type, frames,
                       save_aligned_quads)
        
        quads = full_frame_to_quads(full_frame, nx_quad, ny_quad)
        active_pixel_quads = find_active_region_each_quads(quads)       
        title = 'Histogram of active pixels'  
               
        # Creste list of directories To save histograms before and after 
        #outlier rejections
        mean_plot = 'Outlier_mean'
        median_plot = 'Outlier_median'
        folder_name_hist = 'Hist_plot'
        folder_name_mask ='Mask_plot'
        
        
        save_mean_hist = os.path.join(plot_dir,mean_plot,folder_name_hist)
        if not os.path.exists(save_mean_hist ):os.makedirs(save_mean_hist)        
        
        save_median_hist = os.path.join(plot_dir,median_plot,folder_name_hist)
        if not os.path.exists(save_median_hist):os.makedirs(save_median_hist)
        
        save_mean_mask = os.path.join(plot_dir,mean_plot,folder_name_mask)
        if not os.path.exists(save_mean_mask):os.makedirs(save_mean_mask)        
        
        save_median_mask = os.path.join(plot_dir,median_plot,folder_name_mask)
        if not os.path.exists(save_median_mask):os.makedirs(save_median_mask)
        
        # Now let's call the outlier functions        
        #outlier_filt_mean, outlier_mean = reject_outlier_mean(active_pixel_quads)
        #outlier_filt_med, outlier_med = reject_outlier_median(active_pixel_quads)
         
        # Ok now lets' plot the histograms to understand what we see. 
        # First the mean approach
        #print('ok, let us plot the histograms now')
        #plot_hist_image(active_pixel_quads, title, collection_type, frames, 
                        #outlier_filt_mean, outlier_mean,save_mean_hist)
        # Secondly the median absolute difference approach
        #plot_hist_image(active_pixel_quads, title, collection_type, frames, 
                        #outlier_filt_med, outlier_med,save_median_hist)
       
    if __name__ == "__main__":
    main()
