# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""


import os
import numpy as np
import pandas as pd
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from make_quads_from_raw_fpe import make_quads_from_IDL_sav_file,\
                                    make_full_frame_from_raw_fpe,\
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

from outlier_detection import identify_saturation,\
                              reject_outlier_median,\
                              reject_outlier_mean,\
                              create_outlier_mask,\
                              create_final_mask
def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def filter_outlier_median(quads):
    ndims, nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 4]
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

def plot_hist(data, radiation_test, plot_dir, figure_name, quad, temp, int_time):
    
    nrows = 1 # to create 2*2 plot for each quads in each direction
    ncols = 1
    #figure_name = 'hist_'+str(radiatio)
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols)
    mean = np.mean(data)
    sigma = np.std(data)
    med = np.median(data)
    max_val = np.max(data)
    min_val = np.min(data)
    label = 'Mean = '+ str(round(mean, 1)) + \
            '\n Median = '+ str(round(med, 2)) + \
            '\n Std. = '+ str(round(sigma, 2))+ \
            '\n Max = '+ str(round(max_val, 2)) + \
            '\n Min = '+ str(round(min_val, 1)) + \
             '\n Temp = '+ str(round(temp, 2))+'K' + \
            '\n Integ. Time = '+ str(int_time)+ ' microsecs'                

    title = 'Histograms of Dark Current (' + str(radiation_test) +')'
    #print(title)
    sns.set_context("talk")
    with sns.axes_style("darkgrid"):
        
        ax.hist(data, 50, normed=0, facecolor='blue', alpha=1.75, label=label)
        ax.tick_params(axis='x', pad=10)
        ax.grid(True, linestyle=':')
        legend = ax.legend(loc='best', ncol=3, shadow=True,
                           prop={'size':12}, numpoints=1)
        legend.get_frame().set_edgecolor('r')
        legend.get_frame().set_linewidth(2.0)
        ax.set_ylabel('Frequency (# of pixels)', fontsize=15,
                      fontweight="bold")
        #ax.set_xlim(10000, 14000)
        ax.set_xlabel('Counts (DNs)', fontsize=14, fontweight="bold")        
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        fig.savefig(plot_dir+'/'+'hist'+quad+ figure_name+'.png', dpi=100)
        plt.close('all')
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo'
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    
     
    
    
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark_Current_Radiation_testing'
    temp_folder = os.listdir(file_path)
    
    
    for files in range(0, len(temp_folder)):        
        saved_quads = os.path.join( file_path, temp_folder[files])        
        
        nominal_int_files = [each for each in os.listdir(saved_quads) \
                             if each.endswith('.fits')]
      
        all_int_time = [ ]
        all_temp = [ ]
        all_med_quad_A = [ ]
        all_med_quad_B = [ ]
        
        all_var_quad_A = [ ]
        all_var_quad_B = [ ]
        
        dframe1 = pd.DataFrame()
        save_dir = r'dark_current_analysis'
        plot_dir = os.path.join(saved_quads, save_dir)
        
        if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
              
        for data_files in nominal_int_files:  
            data_path_name_split = data_files.split('_')           
            #integration_type = [x for x in names_to_check if x in data_path_name_split]
            #collection_type = integration_type[0]
            #frames_int = data_path_name_split[-1]            
            temp = round(int(data_path_name_split[-2].split('-')[1][0:3])) #for integ_sweep  
            temp = temp-273.15
            int_time = round(int(data_path_name_split[-3].split('-')[1]))
            all_temp.append(temp)
            all_int_time.append(int_time)
           
            data_file = os.path.join(saved_quads, data_files)  
            hdulist = fits.open(data_file)
            full_frame = hdulist[0].data
                            
            ndimes, nx_quad, ny_quad = full_frame.shape             
            quad_A = full_frame[0:ndimes:2, :, :]
            
            active_quad_A = np.array(quad_A[:, 4:1028, 10:1034])            
            bias_A =  bias_A = quad_A[:, 4:1028, 1034:1056] 
            bias_subtracted_quad_A = perform_bias_subtraction (active_quad_A, bias_A)      
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_A)
            var_quad_A = np.var(outlier_filtered_data)     
            median_quad_A = np.mean(outlier_filtered_data)            
            all_med_quad_A.append(median_quad_A)
            all_var_quad_A.append(var_quad_A)
            quad = 'Quad_A'
            plot_hist(outlier_filtered_data, temp_folder[files], plot_dir, data_files,quad, temp, int_time)
                        
            # For quad B
            quad_B = full_frame[1:ndimes:2, :, :]
            active_quad_B = np.array(quad_B[:, 4:1028, 10:1034])            
            bias_B = quad_B[:, 4:1028, 1034:1056] 
            bias_subtracted_quad_B = perform_bias_subtraction (active_quad_B, bias_B)      
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_B)
            var_quad_B = np.var(outlier_filtered_data)     
            median_quad_B = np.mean(outlier_filtered_data)            
            all_med_quad_B.append(median_quad_B)
            all_var_quad_B.append(var_quad_B)
            quad = 'Quad_B'
            plot_hist(outlier_filtered_data, temp_folder[files], plot_dir, data_files,quad, temp, int_time)
            
           
            
            

       
        print('A:', all_med_quad_A)
        print('B:', all_med_quad_B)
        
         
        dframe1 = pd.DataFrame(
                    {'Int_time.' : all_int_time,
                     'Temp.' : all_temp,
                     'Avg. Quad A' : all_med_quad_A,
                     'Avg. Quad B' : all_med_quad_B,                     
                     'Avg. Var_Quad_A': all_var_quad_A,
                     'Avg. Var_Quad_B': all_var_quad_B
                     
                     
                     })
                
        dframe1.to_csv(plot_dir +'/'+'Median_DN_vs_Var.csv')
        #dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test\Unct_Quad_A_flooded.csv')
        
if __name__ == "__main__":
    main()