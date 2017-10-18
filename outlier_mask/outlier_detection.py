# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 09:44:14 2017

@author: nmishra

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from copy import deepcopy


def identify_saturation(quads, plot_dir,collection_type, frames, int_time):
    """ This function identifies the saturated pixels in the quads. 
    Ideally saturation level is 2^N-1. But for additional threshold has been 
    added for analysis purpose
        
    """
    nx_quad = quads.shape[0]
    ny_quad = quads.shape[1]
    data = np.reshape(quads, (nx_quad*ny_quad, 1))
    #print(np.max(data[data<16000]))
    saturated_pixels = np.where(data>16000)[0]
    #print(saturated_pixels)
    
    mask = np.zeros((nx_quad * ny_quad, 1)).astype(int)
    mask[saturated_pixels] = 1
    mask = np.reshape(mask, (nx_quad, ny_quad))
    sat_mask = deepcopy(mask)
    
    plt.figure()
    plt.imshow(np.invert(mask), cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title('Saturation mask (Sat. Pixels = '+str(len(saturated_pixels))+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(plot_dir+'/'+ collection_type+'_'+\
                    frames+'.png', dpi=100, bbox_inches="tight")
    #plt.show()   
    plt.close('all')
    return  sat_mask


def reject_outlier_median(quads, sigma=4.):
    """ This function does a preliminary outlier classification of the data
    Median based filtering is more robust than mean based filter. The mean of
    a distribution is more biased by the outlier than the median. The algorithm
    steps are
    a) Find the relative distance of the data point to the median.
    b) Find the median distance of a.
    c) scale the relative distance by the median so that the threshold sigma is
        on a reasonable relative scale.
    """
    nx_quad = quads.shape[0]
    ny_quad = quads.shape[1]
    # All the print statement does intermediate sanity checks
    #print(ny_quad)
    #print(nx_quad)
    data = np.reshape(quads, (nx_quad*ny_quad, 1))
    #print('before= ', np.max(data))
    #data = data[data < (2**14-1)]
    #print('after= ', np.max(data))
    diff = abs(data - np.median(data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = data[measured_threshold < sigma]
    outlier_detectors = np.array(np.where([measured_threshold > sigma]))
    #print((np.array(outlier_detectors).shape)[1])
    #print('median->',  outlier_detectors.shape)

    return outlier_filtered_data, outlier_detectors

def reject_outlier_mean(quads, sigma=5.):
    """ Outlier classification based on mean.
    """
    nx_quad = quads.shape[0]
    ny_quad = quads.shape[1]
    data = (np.reshape(quads, (nx_quad*ny_quad, 1)))
    #data = data[data < (2**14-1)]
    #diff = np.abs(data - np.mean(data))
    outlier_filtered_data = data[abs(data - np.mean(data)) < sigma * np.std(data)]
    outlier_detectors = np.where(abs(data - np.mean(data)) > sigma * np.std(data))
    #print((np.array(outlier_detectors).shape))

    #print (data[diff > sigma * np.std(data)] )
    #print('mean->',np.where(diff > sigma * np.std(data)))
    #print('Mean Method outliers',np.array(outlier_detectors).shape)
    return outlier_filtered_data, outlier_detectors


def create_outlier_mask(quads, outlier_pixels, tit, collection_type, frames,
                        plot_dir):
    """Create a binary mask of the quads where 0 = outlier pixels
        and save the binary mask as png files and .mat files for future use.
        Location of mat files are created inside the function and this location
        is returned to the main function
    """
        
    num_outlier = 'outliers = '+ str(np.count_nonzero(outlier_pixels))
    
    title = tit + ' ('+ num_outlier +')'
    nx_quad = quads.shape[0]
    ny_quad = quads.shape[1]    
    #print(outlier_pixels)
    mask = np.zeros((nx_quad * ny_quad, 1)).astype(int)
    mask[outlier_pixels] = 1
    mask = np.reshape(mask, (nx_quad, ny_quad))
    out_mask = deepcopy(mask)
    #scipy.io.savemat(save_mask_mat_file+'/'+ collection_type+'_'+\
     #               frames+'.mat', mdict={'mask': out_mask})
    plt.figure()
    plt.imshow(np.invert(mask), cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title(title, fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(plot_dir+'/'+ collection_type+'_'+\
                    frames+'.png', dpi=400, bbox_inches="tight")
    plt.close('all')
    return  out_mask

def create_final_mask(outlier_mask, quad_name, title, final_path):
    """
    For each outlier mask, count how many times the pixel is bad. If the pixel
    appears as outlier more than 80% of the times, flag that as bad pixel in 
    the final mask
    """
    final_mask = [ ]    
    dim, nx_quad, ny_quad = np.array(outlier_mask).shape
    for i in range(0, nx_quad*ny_quad):
        quads = [a.reshape(nx_quad*ny_quad,1) for a in outlier_mask ]
        lst = [item[i] for item in quads]
       # print(len(lst))
        outliers = lst.count(1)/len(lst)
        #print(outliers)
        
        final_mask.append(outliers)
    for index, items in enumerate(final_mask):
        if not(items< 0.3):
            final_mask[index] = 0
        else:
            final_mask[index] = 1
    final_outliers_num = final_mask.count(0)
    #print(final_outliers_num)
   #cc            
    final_mask = np.reshape(final_mask,(nx_quad, ny_quad))
    plt.figure()
    plt.imshow(final_mask, cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title(title+' (outliers = '+ str(final_outliers_num)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)     
    plt.savefig(final_path+'/'+quad_name+'.png', dpi=200, bbox_inches="tight")
         
    return np.invert(final_mask)
    
def create_ORed_mask(outlier_mask, quad_name, title, final_path):
    
    final_mask = np.bitwise_or.reduce(outlier_mask[:])
    nx_quad, ny_quad = np.array( final_mask).shape
    final_outliers = np.array(np.where(np.reshape(final_mask,
                                          (nx_quad*ny_quad, 1))==1)).shape[1]
    plt.figure()
    plt.imshow(np.invert(final_mask), cmap='bwr', interpolation='none', origin='lower')
    #cbar = plt.colorbar(image)
    plt.title(title+'  (outliers = '+ str(final_outliers)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    
   
    plt.savefig(final_path+'/'+quad_name+'.png', dpi=400, bbox_inches="tight")
    return final_mask