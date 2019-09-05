# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""

import os
import numpy as np
import pandas as pd
#from scipy.io.idl import readsav
#import matplotlib.pyplot as plt
#from matplotlib.ticker import FormatStrFormatter
from scipy.io.idl import readsav
# pylint: disable=E1101




def perform_bias_subtraction(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    even_detector_bias = trailing_overclocks[:, :, 1::2]
    even_detector_bias = even_detector_bias[:, :, 4:-1]
    avg_bias_even = np.mean(even_detector_bias, axis=2)

    odd_detector_bias = trailing_overclocks[:, :, ::2]
    odd_detector_bias = odd_detector_bias[:, :, 4:-1]
    avg_bias_odd = np.mean(odd_detector_bias, axis=2)


    even_detector_active_quad = active_quad[:, :, 1::2]
    odd_detector_active_quad = active_quad[:, :, ::2]

    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, :, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, :, None]
    active_quad[:, :, 1::2] = bias_subtracted_quad_even
    active_quad[:, :, ::2] = bias_subtracted_quad_odd

    return active_quad


def perform_smear_subtraction(active_quad, int_time):
    """the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    """
    frame_transfer = 8.333
    smear_factor = (frame_transfer / (int_time+ frame_transfer))* np.mean(active_quad, axis=0)
    #print(active_quad.shape)
    #print(smear_factor.shape)
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad



def filter_outlier_median(quads):
    """ Apart from the fixed mask, there are times when the outlier needs to be
    run in order to get good statistics. This will be used in conjunction to
    the fixed mask """
    if np.array(quads).ndim == 3:
        ndims, nx_quad, ny_quad = quads.shape
    elif np.array(quads).ndim == 2:
        ndims = 1
        nx_quad, ny_quad = quads.shape
    else:
        nx_quad = 1
        ndims = 1
        ny_quad = len(quads)
    hist_data = np.reshape(quads, (ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 10.]
    return outlier_filtered_data

def main():
    """
    Tme main function
    """
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    save_dir_local_image = r'C:\Users\nmishra\Workspace\TEMPO\TEMPO_DC_FPA_TEMP_SENSTIVITY_NEW\CCD_based'
    temperature_files = [each for each in os.listdir(file_path) \
                        if each.endswith('_PT_Dark')]
    for k in range(1, len(temperature_files)):
        #temp = temperature_files[k][0:4]
        image_data_files = os.path.join(file_path, temperature_files[k],
                                        'Script_Data', 'saved_quads')
        all_int_files_image = [each for each in os.listdir(image_data_files) \
                               if not each.endswith('_118000.dat.sav')
                               if each.endswith('.dat.sav')]
        data_path_op_int = [each for each in os.listdir(image_data_files) \
                           if each.endswith('118000.dat.sav')]
        data_path_op_int = [data_path_op_int[4], data_path_op_int[5],
                            data_path_op_int[6]]
        all_int_files_image.append(data_path_op_int)
        save_dir = os.path.join(save_dir_local_image, temperature_files[k])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        all_int_time = []
        all_med_quad_A = []
        all_med_quad_B = []
        all_med_quad_C = []
        all_med_quad_D = []
#        saturated_collects = ['FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
#                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
#                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
#                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
#                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
#                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
#                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
#                              'FT6_LONG_INT_199999.dat.sav', 'FT6_LONG_INT_125047.dat.sav']

#        nominal_int_files = [items for items in all_int_files_image
#                             if not items.endswith(tuple(saturated_collects))
#                             if items in all_int_files_image]
        for data_files in data_path_op_int:
            data_path_name_split = data_files.split('_')
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
            int_time = int(int_time)/1000
            data_file = os.path.join(image_data_files, data_files)
            print(data_file)
            IDL_variable = readsav(data_file)
            all_full_frame = IDL_variable.q
            all_int_time.append(int_time)
            for i in range(0, 4): # 4 quads
                quad = all_full_frame[:, i, :, :]
                active_quad = quad[:, 2:1028, 10:1034]
                tsoc = quad[:, 2:1028, 1034:1056]
                tsoc = np.array(tsoc)
                bias_subtracted_quad = perform_bias_subtraction(np.array(active_quad),
                                                                np.array(tsoc))
                
                bias_subtracted_quad = perform_smear_subtraction(bias_subtracted_quad,
                                                                 int_time)
                bias_subtracted_quad = bias_subtracted_quad/int_time
                
                
                if i == 0:
                    all_med_quad_A.append(np.mean(filter_outlier_median
                                                  (bias_subtracted_quad)))
#                    print(all_med_quad_A)
#                    cc
                   
                elif i == 1:
                    all_med_quad_B.append(np.mean(filter_outlier_median
                                                  (bias_subtracted_quad)))
                elif i == 2:
                    all_med_quad_C.append(np.mean(filter_outlier_median
                                                  (bias_subtracted_quad)))
                else:
                    all_med_quad_D.append(np.mean(filter_outlier_median
                                                  (bias_subtracted_quad)))
        dframe1 = pd.DataFrame({'Int_time.' : all_int_time,
                                'Avg_Quad_A' : all_med_quad_A,
                                'Avg_Quad_B' : all_med_quad_B,
                                'Avg_Quad_C' : all_med_quad_C,
                                'Avg_Quad_D' : all_med_quad_D
                               })
        
        dframe1.to_csv(save_dir+'/'+temperature_files[k]+'_DC_all.csv')

if __name__ == "__main__":
    main()
