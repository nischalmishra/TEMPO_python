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
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2


    #outlier_mask = read_outlier_mask()
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Photon_transfer_analysis\Electronics_Gain\Electronics_Gain_Linearity_Data'

    image_data_files = os.path.join(file_path, 'saved_quads')
    all_int_files_image = [each for each in os.listdir(image_data_files) \
                         if each.endswith('.dat.sav')]
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    all_int_time = []
    all_med_quad_a_odd = []
    all_med_quad_b_odd = []
    all_med_quad_c_odd = []
    all_med_quad_d_odd = []
    all_med_quad_a_even = []
    all_med_quad_b_even = []
    all_med_quad_c_even = []
    all_med_quad_d_even = []

    all_std_quad_a_odd = []
    all_std_quad_b_odd = []
    all_std_quad_c_odd = []
    all_std_quad_d_odd = []
    all_std_quad_a_even = []
    all_std_quad_b_even = []
    all_std_quad_c_even = []
    all_std_quad_d_even = []

    saturated_collects = ['FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                          'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                          'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                          'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                          'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                          'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                          'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                          'FT6_LONG_INT_199999.dat.sav', 'FT6_LONG_INT_125047.dat.sav']
    nominal_int_files = [items for items in all_int_files_image
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files_image]
    for data_files in nominal_int_files:
        data_path_name_split = data_files.split('_')
        int_time = round(int(data_path_name_split[-1].split('.')[0]))
        int_time = int(int_time)/1000
        data_file = os.path.join(image_data_files, data_files)
        print(data_file)
        idl_variable = readsav(data_file)
        all_full_frame = idl_variable.q
        all_int_time.append(int_time)
        for i in range(0, 4): # 4 quads
            quad = all_full_frame[:, i, :, :]
            active_quad = quad[:, 2:1028, 10:1034]
            tsoc = quad[:, 2:1028, 1034:1056]
            tsoc = np.array(tsoc)
            bias_subtracted_quad = perform_bias_subtraction(np.array(active_quad),
                                                            np.array(tsoc))
            #bias_subtracted_quad = bias_subtracted_quad[:, 100:900, 99:899]
            bias_subtracted_quad_odd = bias_subtracted_quad[:, :, ::2]
            bias_subtracted_quad_even = bias_subtracted_quad[:, :, 1::2]
            variance_odd = np.var(bias_subtracted_quad_odd, axis=0)
            variance_odd = filter_outlier_median(variance_odd)
            variance_even = np.var(bias_subtracted_quad_even, axis=0)
            variance_even = filter_outlier_median(variance_even)
            # Now calculate the read noise
            tsoc = quad[:, 2:1028, 1034:1056]
            tsoc_odd = tsoc[:, :, ::2]
            tsoc_odd = tsoc_odd[:, :, 4:-1]
            variance_tsoc_odd = np.var(tsoc_odd, axis=0)

            tsoc_even = tsoc[:, :, 1::2]
            tsoc_even = tsoc_even[:, :, 4:-1]
            variance_tsoc_even = np.var(tsoc_even, axis=0)

            net_variance_odd = np.mean(variance_odd) - np.mean(variance_tsoc_odd)
            net_variance_even = np.mean(variance_even) - np.mean(variance_tsoc_even)


            if i == 0:
                all_med_quad_a_odd.append(np.mean(filter_outlier_median
                                                  (bias_subtracted_quad_odd)))
                all_med_quad_a_even.append(np.mean(filter_outlier_median
                                                   (bias_subtracted_quad_even)))
                all_std_quad_a_odd.append(net_variance_odd)
                all_std_quad_a_even.append(net_variance_even)

            elif i == 1:
                all_med_quad_b_odd.append(np.mean((bias_subtracted_quad_odd)))
                all_med_quad_b_even.append(np.mean((bias_subtracted_quad_even)))
                all_std_quad_b_odd.append(net_variance_odd)
                all_std_quad_b_even.append(net_variance_even)

            elif i == 2:
                all_med_quad_c_odd.append(np.mean((bias_subtracted_quad_odd)))
                all_med_quad_c_even.append(np.mean((bias_subtracted_quad_even)))
                all_std_quad_c_odd.append(net_variance_odd)
                all_std_quad_c_even.append(net_variance_even)

            else:
                all_med_quad_d_odd.append(np.mean((bias_subtracted_quad_odd)))
                all_med_quad_d_even.append(np.mean((bias_subtracted_quad_even)))
                all_std_quad_d_odd.append(net_variance_odd)
                all_std_quad_d_even.append(net_variance_even)

    dframe1 = pd.DataFrame({'Int_time.' : all_int_time,
                            'Avg_Quad_A_odd' : all_med_quad_a_odd,
                            'Avg_Quad_A_even' : all_med_quad_a_even,
                            'Avg_Quad_B_odd' : all_med_quad_b_odd,
                            'Avg_Quad_B_even' : all_med_quad_b_even,
                            'Avg_Quad_C_odd' : all_med_quad_c_odd,
                            'Avg_Quad_C_even' : all_med_quad_c_even,
                            'Avg_Quad_D_odd' : all_med_quad_d_odd,
                            'Avg_Quad_D_even' : all_med_quad_d_even,
                            'Var_quad_A_odd' :  all_std_quad_a_odd,
                            'Var_quad_A_even' :  all_std_quad_a_even,
                            'Var_quad_B_odd' :  all_std_quad_b_odd,
                            'Var_quad_B_even' :  all_std_quad_b_even,
                            'Var_quad_C_odd' :  all_std_quad_c_odd,
                            'Var_quad_C_even' :  all_std_quad_c_even,
                            'Var_quad_D_odd' :  all_std_quad_d_odd,
                            'Var_quad_D_even' :  all_std_quad_d_even
                           })
    dframe1.to_csv(save_dir+'/'+'Photon_transfer_data_all.csv')





if __name__ == "__main__":
    main()
