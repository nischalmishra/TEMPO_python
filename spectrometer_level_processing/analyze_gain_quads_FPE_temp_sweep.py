# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""

import os
import numpy as np
import pandas as pd
#from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.io.idl import readsav




def perform_bias_subtraction(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    ndims, nx_quad, ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[[0]*ndims]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, :, ::2]
    # remove outliers
    # First 4 hot lines in even and odd
    # last odd lne in odd
    even_detector_bias = even_detector_bias[:, :, 4:]
    avg_bias_even = np.mean(even_detector_bias, axis=2)
    odd_detector_bias = trailing_overclocks[:, :, 1::2]
    odd_samples = odd_detector_bias[:, :, 4:]
    ndims, rows, cols = odd_samples.shape
    odd_detector_bias = odd_samples[ :, :, 0:cols-1]
    avg_bias_odd = np.mean(odd_detector_bias, axis=2)
    even_detector_active_quad = active_quad[:, :, ::2]
    odd_detector_active_quad = active_quad[:, :, 1::2]
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, :, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, :, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ndims, nx_quad, ny_quad))
    bias_subtracted_quad[:, :, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, :, 1::2] = bias_subtracted_quad_odd
    #print(bias_subtracted_quad.shape)
    return bias_subtracted_quad, even_detector_bias, odd_detector_bias



def perform_bias_subtraction_ave(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    nx_quad, ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, ::2]
    # remove outliers
    # First 4 hot lines in even and odd
    # last odd lne in odd
    even_detector_bias = even_detector_bias[:, 4:]
   # print(even_detector_bias.shape)
    avg_bias_even = np.mean(even_detector_bias, axis=1)
    odd_detector_bias = trailing_overclocks[:, 1::2]
    odd_samples = odd_detector_bias[:, 4:]
    rows, cols = odd_samples.shape
    odd_detector_bias = odd_samples[:, 0:cols-1]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
#    print(odd_detector_bias.shape)
#    cc
    even_detector_active_quad = active_quad[:, ::2]
    odd_detector_active_quad = active_quad[:, 1::2]
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (nx_quad, ny_quad))
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad,

def calculate_dark_current(image, temp, i, int_time):
    """ Calculate the dark current based off the dark data
    It takes the filename of Light data and searches for the mathing integration
    time   in the dark data directory, finds the required dark data, subtracts off
    the offset and smear and computes the dark current
    """
    dark_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    tem_file = temp +r'_PT_Dark\Script_Data\saved_quads'
    data_file = os.path.join(dark_data_dir, tem_file)
    data_path_name_split = image.split('_')
    all_int_files = [each for each in os.listdir(data_file) \
                         if each.endswith('_'+data_path_name_split[-1])]
    #print(all_int_files)

    dark_data_file = os.path.join(data_file, all_int_files[0])
    IDL_variable = readsav(dark_data_file)
    all_full_frame = IDL_variable.q
    quad = all_full_frame[:, i, :, :]

    active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)
    tsoc = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
    dark_current_quad = perform_bias_subtraction_ave(active_quad, tsoc)
    return dark_current_quad


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

def plot_read_noise(bias_subtracted_quad, plot_save_dir, title, color, quads):
    """ PLot the read noise histogram
    """
    bias_subtracted_quad_odd = bias_subtracted_quad[:, :, 1::2]
    bias_subtracted_quad_even = bias_subtracted_quad[:, :, ::2]

    variance_odd = np.var(bias_subtracted_quad_odd, axis=0)
    variance_odd = filter_outlier_median(variance_odd)

    variance_even = np.var(bias_subtracted_quad_even, axis=0)
    variance_even = filter_outlier_median(variance_even)

    label1 = 'Mean Variance (Odd) = '+str(round(np.mean(variance_odd), 2)) +\
             '\n Std(Odd) = '+  str(round(np.std(variance_odd), 2))

    label2 = 'Mean Variance (Even) = '+str(round(np.mean(variance_even), 2)) +\
             '\n Std (Even) = '+  str(round(np.std(variance_even), 2))


    plt.figure(figsize=(7, 5))
    plt.hist(variance_odd, 50, normed=0, rwidth=1, facecolor=color, alpha=0.8,
             label=label1)
    plt.hist(variance_even, 50, normed=0, rwidth=1, facecolor='black', alpha=0.8,
             label=label2)
    plt.title(title, fontsize=14, fontweight="bold")
    plt.grid(True, linestyle=':')
    plt.xlabel('Temporal Noise (Variance['+r'$DN^{2}$'+'])', fontsize=14,
               fontweight="bold")
    plt.ylabel('Frequency', fontsize=14, fontweight="bold")
    plt.legend(loc='best')
    ax = plt.gca()
    #plt.xlim(0, 300)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #plt.show()
    plt.savefig(plot_save_dir, dpi=100)
    plt.close('all')


def plot_signal_level(bias_subtracted_quad, plot_save_dir, title, color, quads):
    """ PLot the  noise histogram
    """

    bias_subtracted_quad_odd = bias_subtracted_quad[:, :, 1::2]
    bias_subtracted_quad_even = bias_subtracted_quad[:, :, ::2]
    #mean_odd = np.mean(bias_subtracted_quad_odd, axis=0)
    mean_odd = bias_subtracted_quad_odd
    mean_odd = filter_outlier_median(mean_odd)
    #mean_even = np.mean(bias_subtracted_quad_even, axis=0)
    mean_even = bias_subtracted_quad_even
    mean_even = filter_outlier_median(mean_even)
    label1 = 'Mean (Odd) = '+str(round(np.mean(mean_odd), 2)) +\
             '\n Std(Odd) = '+  str(round(np.std(mean_odd), 2))

    label2 = 'Mean  (Even) = '+str(round(np.mean(mean_even), 2)) +\
             '\n Std (Even) = '+  str(round(np.std(mean_even), 2))

    plt.figure(figsize=(7, 5))
    plt.hist(filter_outlier_median(mean_odd), 100, rwidth=1, normed=0,
             facecolor=color, alpha=0.8, label=label1)
    plt.hist(filter_outlier_median(mean_even), 100, rwidth=1, normed=0,
             facecolor='black', alpha=0.8, label=label2)
    plt.title(title, fontsize=14, fontweight="bold")
    plt.grid(True, linestyle=':')
    plt.xlabel('Signal-Offset (DN)', fontsize=14, fontweight="bold")
    plt.ylabel('Frequency', fontsize=14, fontweight="bold")
    plt.legend(loc='best')
    ax = plt.gca()
    #plt.xlim(0, 300)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #plt.show()
#    cc
    plt.savefig(plot_save_dir, dpi=100)
    plt.close('all')



def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2


    #outlier_mask = read_outlier_mask()
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPE_Gain_vs_Temp'
    save_dir_local_image = r'C:\Users\nmishra\Workspace\TEMPO\TEMPO_FPE_TEMP_SENSTIVITY_LINEARITY'

    temperature_files = [each for each in os.listdir(file_path) \
                        if each.endswith('_PT_Light')]

    
    for k in range(0, len(temperature_files)):
        image_data_files = os.path.join(file_path, temperature_files[k],
                                        'Script_Data', 'saved_quads')
        
        data_path_all = [each for each in os.listdir(image_data_files)
                         if not each.endswith('_118000.dat.sav')
                         if each.endswith('.dat.sav')
                         ]
        
        data_path_op_int = [each for each in os.listdir(image_data_files) \
							if each.endswith('118000.dat.sav')]
        
        #data_path_all.append(data_path_op_int)
        temp = temperature_files[k][0:4]
        save_dir = os.path.join(save_dir_local_image, temperature_files[k])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
#        print(save_dir)
#        cc
        all_int_time = []
        all_med_quad_A_odd = []
        all_med_quad_B_odd = []
        all_med_quad_C_odd = []
        all_med_quad_D_odd = []
        all_med_quad_A_even = []
        all_med_quad_B_even = []
        all_med_quad_C_even = []
        all_med_quad_D_even = []

        all_var_quad_A_odd = []
        all_var_quad_B_odd = []
        all_var_quad_C_odd = []
        all_var_quad_D_odd = []
        all_var_quad_A_even = []
        all_var_quad_B_even = []
        all_var_quad_C_even = []
        all_var_quad_D_even = []
        saturated_collects = ['FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav', 'FT6_LONG_INT_125047.dat.sav']

        nominal_int_files = [items for items in data_path_all
                         if not items.endswith(tuple(saturated_collects))
                         if items in data_path_all]
        for data_files in nominal_int_files:

            data_path_name_split = data_files.split('_')
            #print(data_path_name_split)
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
            int_time = int(int_time)/1000
            #print('integ. time= ', int_time)
            data_file = os.path.join(image_data_files, data_files)
            print(data_file)
            #print(int_time)
            IDL_variable = readsav(data_file)
            all_full_frame = IDL_variable.q
            all_int_time.append(int_time)
            quads = ['QuadA', 'QuadB', ' QuadC', 'QuadD']
            color = ['blue', 'red', 'green', 'magenta']
            for i in range(0, 4): # 4 quads

                quad = all_full_frame[:, i, :, :]
                active_quad = quad[:, 2:1028, 10:1034]
                tsoc = quad[:, 2:1028, 1034:1056]
                tsoc = np.array(tsoc)
                bias_subtracted_quad, odd_detector_bias, even_detector_bias = perform_bias_subtraction(np.array(active_quad),
                                                                np.array(tsoc))
                
                #dark_current_quad = calculate_dark_current(data_files, temp, i, int_time)
                #bias_subtracted_quad = bias_subtracted_quad - dark_current_quad
                bias_subtracted_quad = bias_subtracted_quad[:, :, :]
                bias_subtracted_quad_odd = bias_subtracted_quad[:, :, ::2]
                bias_subtracted_quad_even = bias_subtracted_quad[:, :, 1::2]                
                
                # for read noise, calculate variance of tsoc
                variance_odd = np.var(bias_subtracted_quad_odd, axis=0)
                variance_odd = filter_outlier_median(variance_odd)
                variance_even = np.var(bias_subtracted_quad_even, axis=0)
                variance_even = filter_outlier_median(variance_even)
                unct_odd = np.mean(variance_odd)
                unct_even = np.mean(variance_even)
                

                if i == 0:
                    all_med_quad_A_odd.append(np.mean(filter_outlier_median
                                                      (bias_subtracted_quad_odd)))
                    all_med_quad_A_even.append(np.mean(filter_outlier_median
                                                       (bias_subtracted_quad_even)))
                    all_var_quad_A_odd.append(unct_odd)
                    all_var_quad_A_even.append(unct_even)

                elif i == 1:
                    all_med_quad_B_odd.append(np.mean((bias_subtracted_quad_odd)))
                    all_med_quad_B_even.append(np.mean((bias_subtracted_quad_even)))
                    all_var_quad_B_odd.append(unct_odd)
                    all_var_quad_B_even.append(unct_even)
                elif i == 2:
                    all_med_quad_C_odd.append(np.mean((bias_subtracted_quad_odd)))
                    all_med_quad_C_even.append(np.mean((bias_subtracted_quad_even)))
                    all_var_quad_C_odd.append(unct_odd)
                    all_var_quad_C_even.append(unct_even)
                else:
                    all_med_quad_D_odd.append(np.mean((bias_subtracted_quad_odd)))
                    all_med_quad_D_even.append(np.mean((bias_subtracted_quad_even)))
                    all_var_quad_D_odd.append(unct_odd)
                    all_var_quad_D_even.append(unct_even)
        dframe1 = pd.DataFrame({'Int_time.' : all_int_time,
                                'Avg_Quad_A_odd' : all_med_quad_A_odd,
                                'Avg_Quad_A_even' : all_med_quad_A_even,
                                'Avg_Quad_B_odd' : all_med_quad_B_odd,
                                'Avg_Quad_B_even' : all_med_quad_B_even,
                                'Avg_Quad_C_odd' : all_med_quad_C_odd,
                                'Avg_Quad_C_even' : all_med_quad_C_even,
                                'Avg_Quad_D_odd' : all_med_quad_D_odd,
                                'Avg_Quad_D_even' : all_med_quad_D_even,
                                'Var_quad_A_odd' : all_var_quad_A_odd,
                                'Var_quad_A_even' : all_var_quad_A_even,
                                'Var_quad_B_odd' : all_var_quad_B_odd,
                                'Var_quad_B_even' : all_var_quad_B_even,
                                'Var_quad_C_odd' : all_var_quad_C_odd,
                                'Var_quad_C_even' : all_var_quad_C_even,
                                'Var_quad_D_odd' : all_var_quad_D_odd,
                                'Var_quad_D_even' : all_var_quad_D_even
                               })
        #csv_save_dir = save_dir+'/' + temperature_files[k]
        dframe1.to_csv(save_dir+'/'+temperature_files[k]+'_read_noise_FPE_Temp.csv')

if __name__ == "__main__":
    main()
