# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:21:16 2017

@author: nmishra
This script is a kind of proof of concept of the outlier detection. 
The concept is for each dark current measurement, 
a. filter out the outlier mask using Median Absolute Deviation
b. Using Logical AND Operation, among all masks, to create a binary mask
c. Find Saturation Mask for each of the dark current measurements and Use
Logical OR operation between various saturation mask to find potential pixels
that tend to saturate and Create a single binary mask
d. Use Logical OR operation between mask from B and C to create a final outlier
mask

"""

import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns

def main():
# Note: The mat files are saved as dictionary. Use keyword 'outlier_mask_all'
# to read the content
    path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads\Outlier_median\save_masks'
    median_mask_lower = 'median_mask_lower.mat'
   # median_mask_upper = 'median_mask_upper.mat'
    sat_mask_lower = 'saturation_mask_lower.mat'
    #sat_mask_upper = 'saturation_mask_upper.mat'

    
    read_outlier_mask_lower = sio.loadmat(os.path.join(path,median_mask_lower))
    outlier_mask_lower = read_outlier_mask_lower['outlier_mask_lower'].astype(int)
    
    #read_outlier_mask_upper = sio.loadmat(os.path.join(path,median_mask_upper))
    #outlier_mask_upper = read_outlier_mask_upper['outlier_mask_upper'].astype(int)
    
    read_sat_mask_lower = sio.loadmat(os.path.join(path,sat_mask_lower))
    sat_mask_lower_quad = read_sat_mask_lower['sat_mask_lower'].astype(int)
    
    #read_sat_mask_upper = sio.loadmat(os.path.join(path,sat_mask_upper))
    #sat_mask_upper_quad = read_sat_mask_upper['sat_mask_upper'].astype(int)
    
    # perform the logical and operation between various dark currents in both upper
    # and lower mask
    
    final_mask_median_lower = np.bitwise_or.reduce(outlier_mask_lower[:])    
    #final_mask_median_upper = np.bitwise_or.reduce(outlier_mask_upper)    
    final_sat_mask_lower = np.bitwise_or.reduce(sat_mask_lower_quad[:])   
    #final_sat_mask_upper = np.bitwise_or.reduce(sat_mask_upper_quad)
    
    
    final_mask_lower = [final_mask_median_lower, final_sat_mask_lower]
    #final_mask_upper = [final_mask_median_upper, final_sat_mask_upper]
    
   
    
    plt.figure(frameon=False)
    plt.imshow((np.invert(np.bitwise_or.reduce(final_mask_lower[:]).astype(int))),
               cmap='bwr', interpolation='none', origin='lower')
    
    final_outliers_lower = np.array(np.where(np.reshape((np.bitwise_or.reduce(final_mask_lower[:])),
                                      (1024*2048, 1))==1)).shape[1]
    
    plt.title('Quad A & B final outlier mask  (outliers = '+ str(final_outliers_lower)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    sns.despine(left=False, bottom=True, right=True)
    plt.savefig(path +'/'+'final_mask_lower.png', dpi=100,
                bbox_inches="tight")
   
    
    plt.figure(frameon=False)
    plt.imshow((np.bitwise_or.reduce(final_mask_upper)),
               cmap='binary', interpolation='none', origin='lower')
    
    final_outliers_upper = np.array(np.where(np.reshape((np.bitwise_or.reduce(final_mask_upper)),
                                      (1024*2048, 1))==1)).shape[1]
    
    plt.title('Quad C & D final outlier mask  (outliers = '+ str(final_outliers_upper)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    sns.despine(left=False, bottom=True, right=True)
    plt.grid(False)
    plt.savefig(path +'/'+'final_mask_upper.png', dpi=100,
                bbox_inches="tight")
    
if __name__ == "__main__":
    main()


