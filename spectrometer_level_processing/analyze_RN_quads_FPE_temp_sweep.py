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
import random



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





def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2


    #outlier_mask = read_outlier_mask()
    
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    save_dir_local_image = r'C:\Users\nmishra\Workspace\TEMPO\Photon_transfer_analysis\Electronics_Gain\Electronics_Gain_Vs_FPA_Temp'
    

    temperature_files = [each for each in os.listdir(file_path) \
                        if each.endswith('_PT_Dark')]
    for k in range(0, len(temperature_files)):
        temp = temperature_files[k][0:4]
        image_data_files = os.path.join(file_path, temperature_files[k],
                                        'Script_Data', 'saved_quads')        

        op_int_files = [each for each in os.listdir(image_data_files) \
                         if each.endswith('_118000.dat.sav')]
        all_int_files_image = [each for each in os.listdir(image_data_files) \
                         if not each.endswith('_118000.dat.sav')]     
        
        op_int_files_random = random.choice(op_int_files)       
        all_int_files_image = all_int_files_image + [op_int_files_random]   


        save_dir = os.path.join(save_dir_local_image, temperature_files[k])

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
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
        all_int_time = []
        all_var_quad_A_odd = []
        all_var_quad_B_odd = []
        all_var_quad_C_odd = []
        all_var_quad_D_odd = []
        all_var_quad_A_even = []
        all_var_quad_B_even = []
        all_var_quad_C_even = []
        all_var_quad_D_even = []
        
        for data_files in nominal_int_files:
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
                tsoc = quad[:, 2:1028, 1034:1056]                
                tsoc_odd = tsoc[:, :, ::2]
                tsoc_odd = tsoc_odd[:,:, 4:-1]               
                variance_tsoc_odd = np.var(tsoc_odd, axis=0)                
                tsoc_even = tsoc[:, :, 1::2]
                tsoc_even = tsoc_even[:,:, 4:-1]
                variance_tsoc_even = np.var(tsoc_even, axis=0)               
                # for read noise, calculate variance of tsoc
                if i == 0:                    
                    all_var_quad_A_odd.append(np.mean(variance_tsoc_odd))
                    all_var_quad_A_even.append(np.mean(variance_tsoc_even))
                elif i == 1:
                    all_var_quad_B_odd.append(np.mean(variance_tsoc_odd))
                    all_var_quad_B_even.append(np.mean(variance_tsoc_even))
                elif i == 2:
                    all_var_quad_C_odd.append(np.mean(variance_tsoc_odd))
                    all_var_quad_C_even.append(np.mean(variance_tsoc_even))
                else:
                    all_var_quad_D_odd.append(np.mean(variance_tsoc_odd))
                    all_var_quad_D_even.append(np.mean(variance_tsoc_even))
                    
        dframe1 = pd.DataFrame({'Int_time.' : all_int_time,                               
                            'Var_quad_A_odd' : all_var_quad_A_odd,
                            'Var_quad_A_even' : all_var_quad_A_even,
                            'Var_quad_B_odd' : all_var_quad_B_odd,
                            'Var_quad_B_even' : all_var_quad_B_even,
                            'Var_quad_C_odd' : all_var_quad_C_odd,
                            'Var_quad_C_even' : all_var_quad_C_even,
                            'Var_quad_D_odd' : all_var_quad_D_odd,
                            'Var_quad_D_even' : all_var_quad_D_even
                               })
       
        dframe1.to_csv(save_dir+'/'+temperature_files[k]+'_read_noise_FPE_Temp.csv')


if __name__ == "__main__":
    main()
