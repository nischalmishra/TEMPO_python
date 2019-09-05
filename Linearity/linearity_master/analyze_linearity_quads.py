# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""

import os
import numpy as np
import pandas as pd
#from scipy.io.idl import readsav
from scipy.io.idl import readsav
import matplotlib.pyplot as plt

def read_outlier_mask():
    """ Read the outlier mask
    """
    outlier_mask = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_outlier_mask_2_sigma.csv', delimiter=',')
    quad_A = outlier_mask[0:1024, 0:1024]
    quad_B = outlier_mask[1024:, 0:1024]
    quad_C = outlier_mask[1024:, 1024:]
    quad_D = outlier_mask[0:1024:, 1024:]
    outlier_mask_final = [quad_A, quad_B, quad_C, quad_D]
    return outlier_mask_final


def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def filter_outlier_median(quads):
    if np.array(quads).ndim == 3:
        ndims, nx_quad, ny_quad = quads.shape
    else:
        ndims = 1
        nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads, (ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 4.]
    return outlier_filtered_data[outlier_filtered_data < 16383] # non-saturated_pixels

def perform_bias_subtraction(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    nx_quad, ny_quad = active_quad.shape    
    odd_detector_bias = trailing_overclocks[ :, ::2]
    odd_detector_bias = odd_detector_bias[:, 4:]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_bias = trailing_overclocks[:, 1::2]
    even_samples = even_detector_bias[:, 4:]
    rows, cols = even_samples.shape
    even_detector_bias = even_samples[:, 0:cols-1]
    avg_bias_even = np.mean(even_detector_bias, axis=1)
    odd_detector_active_quad = active_quad[:, ::2]
    even_detector_active_quad = active_quad[:, 1::2]
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]    
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None] 
    print(np.mean(avg_bias_odd), np.mean(avg_bias_even))
    
    return bias_subtracted_quad_odd, bias_subtracted_quad_even


def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2

   
    #outlier_mask = read_outlier_mask()
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'


    all_int_files = [each for each in os.listdir(file_path) \
                         if each.endswith('.dat.sav')]
    if 'Integration_Sweep' in file_path:
        saturated_collects = ['FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']
    elif 'Intensity_Sweep' in file_path:
        saturated_collects = ['162_OP_INT_118000.dat.sav', '164_OP_INT_118000.dat.sav'
                             ]
    nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files]
    
 
    save_dir = r'subset_linearity'
    plot_dir = os.path.join(file_path, save_dir)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir) 
    
 
    all_int_time = []
    all_med_quad_A_odd = []
    all_med_quad_B_odd = []
    all_med_quad_C_odd = []
    all_med_quad_D_odd = []
    all_med_quad_A_even = []
    all_med_quad_B_even = []
    all_med_quad_C_even = []
    all_med_quad_D_even = []
    
    all_std_quad_A_odd = []
    all_std_quad_B_odd = []
    all_std_quad_C_odd = []
    all_std_quad_D_odd = []
    all_std_quad_A_even = []
    all_std_quad_B_even = []
    all_std_quad_C_even = []
    all_std_quad_D_even = []
    dframe1 = pd.DataFrame()
    for data_files in nominal_int_files[10:]:
            data_path_name_split = data_files.split('_')
            print(data_files)
            if 'Intensity_Sweep' in file_path:
                int_time = data_path_name_split[0]
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                data_file = os.path.join(file_path, data_files)
                IDL_variable = readsav(data_file)
                all_full_frame = IDL_variable.q
                all_int_time.append(int_time)              
                for i in range(0, 4): # 4 quads                    

                    quad = all_full_frame[:, i, :, :]
                    active_quad = np.mean(quad[:, 2:1030, 10:1034], axis=0)                    
                    tsoc = np.mean(quad[:, 2:1030, 1034:1056], axis=0)
                    bias_subtracted_quad_odd, bias_subtracted_quad_even = perform_bias_subtraction(active_quad, tsoc)
                    bias_subtracted_quad_odd = bias_subtracted_quad_odd
                    bias_subtracted_quad_even = bias_subtracted_quad_even
                   
                    #unct_even = 10*np.std(filter_outlier_median(bias_subtracted_quad_even))/(np.median(filter_outlier_median(bias_subtracted_quad_even)))
                    #unct_odd = 10*np.std(filter_outlier_median(bias_subtracted_quad_odd))/(np.median(filter_outlier_median(bias_subtracted_quad_odd)))
                    
                    if i == 0:
                        all_med_quad_A_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                        all_med_quad_A_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                        #all_std_quad_A_even.append(unct_even)
                        #all_std_quad_A_odd.append(unct_odd)
                    elif i == 1:
                        all_med_quad_B_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                        all_med_quad_B_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                        #all_std_quad_B_even.append(unct_even)
                        #all_std_quad_B_odd.append(unct_odd)
    
                    elif i == 2:
                        all_med_quad_C_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                        all_med_quad_C_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                        #all_std_quad_C_even.append(unct_even)
                        #all_std_quad_C_odd.append(unct_odd)
    
                    else:
                        all_med_quad_D_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                        all_med_quad_D_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                        #all_std_quad_D_even.append(unct_even)
                        #all_std_quad_D_odd.append(unct_odd)
                cc
                    
    dframe1 = pd.DataFrame({'Int_time.' : all_int_time,
                            'Avg_Quad_A_odd' : all_med_quad_A_odd,
                            'Avg_Quad_A_even' : all_med_quad_A_even,
                            'Avg_Quad_B_odd' : all_med_quad_B_odd,
                            'Avg_Quad_B_even' : all_med_quad_B_even,
                            'Avg_Quad_C_odd' : all_med_quad_C_odd,
                            'Avg_Quad_C_even' : all_med_quad_C_even,
                            'Avg_Quad_D_odd' : all_med_quad_D_odd,
                            'Avg_Quad_D_even' : all_med_quad_D_even,
                             })

    dframe1.to_csv(plot_dir+'/'+'Linearity_Integration_Sweep_SAO.csv')
            #dframe1 = pd.DataFrame()
                
        #dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test\Unct_Quad_A_flooded.csv')
if __name__ == "__main__":
    main()
