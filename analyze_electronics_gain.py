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
from scipy.io.idl import readsav  
def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def filter_outlier_median(quads):
    if np.array(quads).ndim ==3:
        ndims, nx_quad, ny_quad = quads.shape
    else:
        ndims=1
        nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 3.]
    #print(outlier_filtered_data)
    return outlier_filtered_data

   
def perform_bias_subtraction (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    ndims, nx_quad,ny_quad = active_quad.shape
    
    bias_subtracted_quad = np.array([[[0]*ndims]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[:, :, ::2]
    avg_bias_even = np.mean(even_detector_bias, axis=2)  
    odd_detector_bias = trailing_overclocks[:, :, 1::2]
    avg_bias_odd = np.mean(odd_detector_bias, axis=2)
    even_detector_active_quad = active_quad[:, :, ::2]     
    odd_detector_active_quad = active_quad[:, :, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, :,None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, :, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ndims, ny_quad, nx_quad))    
    bias_subtracted_quad[:,:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:,:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad



def perform_bias_subtraction_ave (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, ::2]
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
    
def plot_few_tsocs(serial_overclocks,plot_dir, quad, int_time):
    # let's take the mean tsoc for 100 frames
    
    mean_tsoc = np.mean(serial_overclocks, axis=0)
    plt.plot(mean_tsoc, '.')
    plt.xlim(0, 1100)
    #plt.ylim(650,  900)
    plt.title('Trailing Overclock profile, '+ quad+ ', Integ.time = '+ str(int_time) +' microsec', 
              fontsize=12, fontweight='bold')
    plt.xlabel('Spectral pixel indices (#)', fontsize=12, fontweight='bold')
    plt.ylabel('Overclock signal level (DN)', fontsize=12, fontweight='bold')
  
    plt.savefig(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\Overclocks\plot_sample'+'/'+ quad+'.png',dpi=100,bbox_inches="tight")
    
    

def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    
    all_int_time = [ ]
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
    
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    nominal_int_files = [each for each in os.listdir(file_path) \
                         if each.endswith('.dat.sav')]
    
    save_dir = r'Trailing_overclocks'
    plot_dir = os.path.join(file_path, save_dir)
    if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
    for data_files in nominal_int_files:           
            
            data_path_name_split = data_files.split('_')  
            print(data_files)
            #integration_type = [x for x in names_to_check if x in data_path_name_split]
            #collection_type = integration_type[0]
            #frames_int = data_path_name_split[-1]            
            #for integ_sweep  
           
            if 'Intensity_Sweep' in file_path:
                 int_time = data_path_name_split[0]
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                    
            data_file = os.path.join(file_path, data_files)  
            IDL_variable = readsav(data_file)            
            all_full_frame = IDL_variable.q  
            all_int_time.append(int_time)
            
            # Let us calculate mean first
            #quads = [all_full_frame[:, 0, :, :], all_full_frame[:, 1, :, :], 
                    #all_full_frame[:, 2, :, :], all_full_frame[:, 3, :, :]]
           
            for i in range(0, 4): # 4 quads
                #print(i)
                quad = all_full_frame[:, i, :, :]
                active_quad_all = quad[:, 4:1028, 10:1034]
                active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)
                bias = np.mean(quad[:, 4:1028, 1034:1056], axis=0)             
                avg_bias_even = np.mean(bias[:,::2], axis=1)  
                avg_bias_odd = np.mean(bias[:, 1::2], axis=1)
                even_detector_active_quad = active_quad[:, ::2]     
                odd_detector_active_quad = active_quad[:, 1::2]    
                bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
                bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
                # subtract smear too
#                tFT = 8*10**(3)
#                ti = int_time
#                smear_factor_even = (tFT / (ti+ tFT))* np.mean(bias_subtracted_quad_even, axis=0)
#                bias_subtracted_quad_even = bias_subtracted_quad_even-smear_factor_even[None, :] 
                #print(smear_factor_even)    
   
              
#                smear_factor_odd= (tFT / (ti+ tFT))* np.mean(bias_subtracted_quad_odd, axis=0)
#                #print(smear_factor_odd)
#                bias_subtracted_quad_odd= bias_subtracted_quad_odd-smear_factor_odd[None, :]   
#                

                # To calculate the std follow the steps below
                active_quad_100th_frame = quad[50, 4:1028, 10:1034]
                active_quad_101th_frame = quad[51, 4:1028, 10:1034]
                
                # Lets estimate noise for even pixels
                bias_100th_frame = quad[50, 4:1028, 1034:1056]
                active_quad_even_100th_frame = active_quad_100th_frame[:, ::2]
                avg_bias_even_100th_frame = np.mean(bias_100th_frame [:, ::2], axis=1)
                bias_sub_100th_frame_even = active_quad_even_100th_frame - avg_bias_even_100th_frame[:, None]
                
                
                bias_101th_frame = quad[51, 4:1028, 1034:1056]
                active_quad_even_101th_frame = active_quad_101th_frame[:,::2]
                avg_bias_even_101th_frame = np.mean(bias_101th_frame [:,::2], axis=1)
                bias_sub_101th_frame_even = active_quad_even_101th_frame - avg_bias_even_101th_frame[:, None]
                diff_two_frames_even = 1000 + bias_sub_101th_frame_even -  bias_sub_100th_frame_even
                sigma_s_even = np.std(diff_two_frames_even)/2
              
                # Lets estimate noise for odd pixels
                active_quad_odd_100th_frame = active_quad_100th_frame[:,1::2]
                avg_bias_odd_100th_frame = np.mean(bias_100th_frame[:,1::2], axis=1)
                bias_sub_100th_frame_odd = active_quad_odd_100th_frame - avg_bias_odd_100th_frame[:, None]
                
                active_quad_odd_101th_frame = active_quad_101th_frame[:,1::2]
                avg_bias_odd_101th_frame = np.mean(bias_101th_frame[:,1::2], axis=1)
                bias_sub_101th_frame_odd = active_quad_odd_101th_frame - avg_bias_odd_101th_frame[:, None]
                
                diff_two_frames_odd = 1000 + bias_sub_101th_frame_odd -  bias_sub_100th_frame_odd
                sigma_s_odd = np.std(diff_two_frames_odd)/2
                
                
                if i == 0: 
                    
                    all_med_quad_A_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                    all_med_quad_A_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                    all_std_quad_A_even.append(sigma_s_even)
                    all_std_quad_A_odd.append(sigma_s_odd)
                
                elif i == 1: 
                    all_med_quad_B_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                    all_med_quad_B_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))              
                    all_std_quad_B_even.append(sigma_s_even)
                    all_std_quad_B_odd.append(sigma_s_odd)
                  
                    
                    
                elif i == 2: 
                    all_med_quad_C_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                    all_med_quad_C_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                    all_std_quad_C_even.append(sigma_s_even)
                    all_std_quad_C_odd.append(sigma_s_odd)
                    
                    
                else: 
                   all_med_quad_D_even.append(np.mean(filter_outlier_median(bias_subtracted_quad_even)))
                   all_med_quad_D_odd.append(np.mean(filter_outlier_median(bias_subtracted_quad_odd)))
                   all_std_quad_D_even.append(sigma_s_even)
                   all_std_quad_D_odd.append(sigma_s_odd)
                  
                    
                active_quad = bias = avg_bias_even = avg_bias_odd = None                
                even_detector_active_quad = odd_detector_active_quad = None
                bias_subtracted_quad_even = bias_subtracted_quad_odd = None
        
    
    dframe1 = pd.DataFrame(
                    {'Int_time.' : all_int_time,
                     'Avg_Quad_A_odd' : all_med_quad_A_odd,
                     'Avg_Quad_A_even' : all_med_quad_A_even,
                     'Avg_Quad_B_odd' : all_med_quad_B_odd,
                     'Avg_Quad_B_even' : all_med_quad_B_even,                       
                     'Avg_Quad_C_odd' : all_med_quad_C_odd,
                     'Avg_Quad_C_even' : all_med_quad_C_even,                                       
                     'Avg_Quad_D_odd' : all_med_quad_D_odd,
                     'Avg_Quad_D_even' : all_med_quad_D_even,
                     'Var_Quad_A_odd ': all_std_quad_A_odd,
                     'Var_Quad_A_even': all_std_quad_A_even,
                     'Var_Quad_B_odd ': all_std_quad_B_odd,
                     'Var_Quad_B_even': all_std_quad_B_even,
                     'Var_Quad_C_odd ': all_std_quad_C_odd,
                     'Var_Quad_C_even': all_std_quad_C_even,
                     'Var_Quad_D_odd ': all_std_quad_D_odd,
                     'Var_Quad_D_even': all_std_quad_D_even,
                     
                     })
                
    
    
    dframe1.to_csv(file_path+'/'+'FPS_gain_analysis_ping_pong.csv')
        #dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test\Unct_Quad_A_flooded.csv')
if __name__ == "__main__":
    main()