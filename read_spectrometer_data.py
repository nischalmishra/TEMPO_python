# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

rewrite: This code was written to process the tempo spectrometer data.
The hardcoded variables are standard to tempo instruments. Basically it's
taken from Dave Flittner's IDL script and pythonized. The sanity checks were
done by comparing the outputs from both the script. This particular script
though reads the IDL .sav file from Dave's code.

"""

import os
import numpy as np
from scipy.io.idl import readsav
import scipy.io as sio
import matplotlib.pyplot as plt
import pandas as pd
#from make_quads_from_raw_fpe import make_quads_from_spectrometer_sav_file

#*****************************************************************************

def read_mat_file(filename):    
    mat_contents = sio.loadmat(filename)
    image_data = mat_contents['image']['img'][0, 0]
    return image_data 


def filter_outlier_median(quads):
    ndims, nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 3]
    return outlier_filtered_data

def perform_bias_subtraction (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([0*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[:, ::2]
    avg_bias_even = np.mean(even_detector_bias, axis=1)  
    odd_detector_bias = trailing_overclocks[:, 1::2]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]     
    odd_detector_active_quad = active_quad[:, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ny_quad, nx_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad


def read_linearity_look_up_table () :
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file= 'Linearity_look_up_table.csv'
    dframe = pd.read_csv(path+'/'+ linearity_file)
    return dframe
    

def perform_non_linearity_correction (active_quad, linearity_file, quad_name):
    
                   
                                                                        
        
    
    
    
    
    

def main():
    """
    Read in the saved IDL variables and makes quad images.
    """
    #To DO : Get feedback from B^2 on quads alignment
    
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\Spectrometer'
    data_path_all = [each for each in os.listdir(file_path)
                     if each.endswith('.sav')]
    for data_path in data_path_all:
        data_file = os.path.join(file_path, data_path)
        saved_variables = readsav(data_file)
        #print(saved_variables)
        #cc
        full_frame = saved_variables.thefpa
        
        
        plt.figure()
        plt.imshow(full_frame, cmap='nipy_spectral', origin='lower', aspect='auto')
        plt.grid(False)
        plt.colorbar()
        plt.show()

if __name__ == "__main__":
    main()
