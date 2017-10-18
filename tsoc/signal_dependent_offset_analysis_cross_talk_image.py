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

from outlier_detection import identify_saturation,\
                              reject_outlier_median,\
                              reject_outlier_mean,\
                              create_outlier_mask,\
                              create_final_mask,\
                              create_ORed_mask 
                              
from analytical_functions import plot_full_frame_image,\
                                 plot_each_quad,\
                                 plot_hist_image,\
                                 plot_hist_each_quad
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
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    #print(outlier_filtered_data)
    return outlier_filtered_data


def plot_hist(even_samples, odd_samples, title, figure_name) :
    
    if np.array(even_samples).ndim ==3:         
        num_dims, nx_quad, ny_quad = even_samples.shape  
    
    elif np.array(even_samples).ndim ==2:      
        
        nx_quad, ny_quad = even_samples.shape
        num_dims = 1
    else:
        nx_quad= 1        
        ny_quad = len(even_samples)
        num_dims=1
    
    
    if np.array(odd_samples).ndim ==3:         
        num_dims1, nx_quad1, ny_quad1 = odd_samples.shape  
    
    elif np.array(odd_samples).ndim ==2:      
        
        nx_quad1, ny_quad1 = odd_samples.shape
        num_dims1 = 1
    else:
        nx_quad1= 1        
        ny_quad1 = len(odd_samples)
        num_dims1=1
        
    mean_diff = np.mean(odd_samples)- np.mean(even_samples) 
    text1 = 'Mean Diff (Odd-Even) = ' +  str(round(mean_diff, 2)) +'DN'
    
    
    plt.figure(figsize=(8, 5))
    plt.hist(np.reshape(even_samples, (num_dims*nx_quad* ny_quad, 1)), facecolor='red', label='Even Lines Samples')
    plt.hist(np.reshape(odd_samples, (num_dims1*nx_quad1* ny_quad1, 1)),facecolor='blue', label='Odd Lines Samples')
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.xlim(800, 835)
    plt.ylim(0, 1000)
    plt.text(845, 300, text1)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' Signal Counts (DN)  ', fontsize=12,
              fontweight="bold")
    plt.title(title)
    #plt.show()
    #cc
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")   
    plt.close('all')
        
   
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
    
def plot_few_tsocs(even_samples_avg, odd_samples_avg, figure_name, title):
    # let's take the mean tsoc for 100 frames
    nrows = 2 
    ncols = 1    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    ax[0].plot(even_samples_avg, '.', label='Even Lines')
   
    ax[0].set_title(title+' (Even Lines)', fontsize=12, fontweight='bold')
    #ax[0].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')    
    ax[0].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')
    ax[0].set_ylim(700, 1300)
    # Now for the odd
    ax[1].plot(odd_samples_avg, '.', label='Odd Lines')   
    ax[1].set_title(title+' (Odd Lines)', fontsize=12, fontweight='bold')
    ax[1].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    ax[1].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')    
    ax[1].set_ylim(700, 1300)   
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    plt.show()
    
    
    
    

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
    
    
    all_avg_tsoc_quad_A_odd = []
    all_avg_tsoc_quad_B_odd = [] 
    all_avg_tsoc_quad_C_odd = []
    all_avg_tsoc_quad_D_odd = [] 
    all_avg_tsoc_quad_A_even = []
    all_avg_tsoc_quad_B_even = [] 
    all_avg_tsoc_quad_C_even = []
    all_avg_tsoc_quad_D_even = [] 
    
    all_std_offset_quad_A_odd = []
    all_std_offset_quad_B_odd = []
    all_std_offset_quad_C_odd = []
    all_std_offset_quad_D_odd = []
    all_std_offset_quad_A_even = []
    all_std_offset_quad_B_even = []
    all_std_offset_quad_C_even = []
    all_std_offset_quad_D_even = []
    
    
    
    
    dframe1 = pd.DataFrame()    
    
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Crosstalk\Channel_B\Script_Data\saved_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Signal_dependent_offset\Cross_Talk_Image\Channel_B'
    all_int_files = [each for each in os.listdir(file_path1) \
                         if each.endswith('.dat.sav')]
#    saturated_collects = [ 'FT6_LONG_INT_130018.dat.sav',
#                              'FT6_LONG_INT_134990.dat.sav',
#                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
#                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
#                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
#                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
#                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
#                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
#                              'FT6_LONG_INT_199999.dat.sav']
    
    for i in range(0, 4):
        
        for data_files in all_int_files:
            data_path_name_split = data_files.split('_')              
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)            
            print(data_files)
                         
            if 'Crosstalk' in file_path1:
                int_time = data_path_name_split[5]
                string1 = 'Input='
                string2 = 'Input = '
            
            print(int_time)
            
            all_int_time.append(int_time)
            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data

            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q            
            all_int_time.append(int_time)
            
            
            all_full_frame = IDL_variable.q
            quad = all_full_frame[:, i, :, :]
            tsoc_all = quad[:, 4:1028, 1034:1056]
            quad_dir = quads[i]
            
            #----------------------------------------------------------------#
            # Ok, let's plot the histogram of saved quads
            
            all_frames_hist = 'all_frames_hist'
            save_dir_image = os.path.join(save_dir, quad_dir, all_frames_hist)
            if not os.path.exists(save_dir_image):
               os.makedirs(save_dir_image)
           
            # separate out even and odd lines
            
            even_samples_all = tsoc_all[:, :, ::2]
            odd_samples_all = tsoc_all[:, :, 1::2]           
            
            even_samples_avg = np.mean(even_samples_all, axis=0)            
            odd_samples_avg = np.mean(odd_samples_all, axis=0)
            
            title = 'Histogram of Serial Overclocks (All 100 Frames)\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name = save_dir_image + '/'+ string1 + str(int_time) + '_image.png'  
            #plot_hist(even_samples_all, odd_samples_all, title, figure_name) 
          
            avg_frames_hist = 'avg_frames_hist'
            save_dir_image = os.path.join(save_dir, quad_dir, avg_frames_hist)
            if not os.path.exists(save_dir_image):
               os.makedirs(save_dir_image)
            title = 'Histogram of Serial Overclocks (Avg. of 100 Frames)\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'                                               
            #plot_hist(even_samples_avg, odd_samples_avg,  title, figure_name) 
            
            
            final_two_lines = 'final_two_lines'
            save_dir_image = os.path.join(save_dir, quad_dir, final_two_lines)
            if not os.path.exists(save_dir_image):
                 os.makedirs(save_dir_image)
            even_samples_used = np.mean(even_samples_avg, axis=1)
            odd_samples_used = np.mean(odd_samples_avg, axis=1)
            
            title = 'Histogram of Serial Overclocks (Avg  for even and odd lines)\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'                                               
            plot_hist(even_samples_used, odd_samples_used,  title, figure_name)
            
            tsoc_profile = 'tsoc_plot'
            save_tsoc_profile  = os.path.join(save_dir, quad_dir, tsoc_profile)
            if not os.path.exists(save_tsoc_profile):
                 os.makedirs(save_tsoc_profile)
            figure_name= save_tsoc_profile + '/'+ string1 + str(int_time) + '_image.png'
            title = 'Profile of Serial Overclocks '+ quads[i]+', ' + string2 + str(int_time)
            plot_few_tsocs(even_samples_avg, odd_samples_avg, figure_name, title)
           
            # Now, let us run the outlier_filter for ping and pong using MAD technique
            odd_even_lines = [even_samples_avg, odd_samples_avg]
            odd_even_lines_name = ['Even Lines','Odd Lines']
            for k in np.arange(0, 2):
                # Lets make directory of 
                #k=1
                median_plot = r'Outlier_median_tsoc'+'/'+ odd_even_lines_name[k]
                folder_name_hist = 'Hist_plot'
                folder_name_mask = 'Mask_plot'
                folder_name_sat = 'Saturation_plot'
                save_median_hist = os.path.join(save_dir,quad_dir, median_plot, folder_name_hist)
                if not os.path.exists(save_median_hist):
                    os.makedirs(save_median_hist)
                save_median_mask = os.path.join(save_dir, quad_dir,median_plot, folder_name_mask)
                if not os.path.exists(save_median_mask):
                    os.makedirs(save_median_mask)             
                save_sat_mask = os.path.join(save_dir, quad_dir,median_plot, folder_name_sat)
                if not os.path.exists(save_sat_mask):
                    os.makedirs(save_sat_mask)
               
                #title = 'Histogram of Serial Overclocks \n '+odd_even_lines_name[k]+', '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
                
                figure_name = save_sat_mask + '/'+ string1 + str(int_time) + '_image.png'  
                sat_pixels = identify_saturation(odd_even_lines[k], int_time, figure_name)    
                outlier_filt_med, outlier_med = reject_outlier_median(odd_even_lines[k])
                if len(outlier_filt_med)==1:
                    outlier_med = 0
                else:
                    outlier_med = outlier_med[1]
                    
                title = 'Binary Outlier Mask of Trailing Overclocks (' + quads[i]+', '+ odd_even_lines_name[k]+', '+str(int_time)+' micro secs)'
                 
                figure_name = save_median_mask + '/'+ string1 + str(int_time) + '_image.png'  
                median_mask = create_outlier_mask(odd_even_lines[k], outlier_med,
                                              title, figure_name)
                
                #title = 'Histogram of Trailing Overclocks (' + quads[i]+', '+ odd_even_lines_name[k]+')'
                #title1 = 'outliers = '+ str(outlier_med.shape[0])
                #figure_name = save_median_hist + '/'+ string1 + str(int_time) + '_image.png'  
                #plot_hist_image(odd_even_lines[k],title,title1, outlier_filt_med, outlier_med, figure_name)
                
    cc
    dframe1 = pd.DataFrame(
                    {'Int_time.' : all_int_time,
                     'Avg_Quad_A_odd' : all_med_quad_A_odd,
                     'Avg_Quad_A_even' : all_med_quad_A_even,
                     'Avg_Offset_Quad_A_odd' : all_avg_tsoc_quad_A_odd,
                     'Avg_Offset_Quad_A_even':all_avg_tsoc_quad_A_even,
                     'Avg_unct_quad_A_odd' :all_std_offset_quad_A_odd, 
                     'Avg_unct_quad_A_even' :all_std_offset_quad_A_even,
                     
                     'Avg_Quad_B_odd' : all_med_quad_B_odd,
                     'Avg_Quad_B_even' : all_med_quad_B_even,
                     'Avg_Offset_Quad_B_odd' : all_avg_tsoc_quad_B_odd,
                     'Avg_Offset_Quad_B_even':all_avg_tsoc_quad_B_even,
                     'Avg_unct_quad_B_odd' :all_std_offset_quad_B_odd, 
                     'Avg_unct_quad_B_even' :all_std_offset_quad_B_even,
                                          
                     'Avg_Quad_C_odd' : all_med_quad_C_odd,
                     'Avg_Quad_C_even' : all_med_quad_C_even,
                     'Avg_Offset_Quad_C_odd' : all_avg_tsoc_quad_C_odd,
                     'Avg_Offset_Quad_C_even':all_avg_tsoc_quad_C_even,
                     'Avg_unct_quad_C_odd' :all_std_offset_quad_C_odd, 
                     'Avg_unct_quad_C_even' :all_std_offset_quad_C_even,
                     
                     'Avg_Quad_D_odd' : all_med_quad_D_odd,
                     'Avg_Quad_D_even' : all_med_quad_D_even,
                     'Avg_Offset_Quad_D_odd' : all_avg_tsoc_quad_D_odd,
                     'Avg_Offset_Quad_D_even':all_avg_tsoc_quad_D_even,
                     'Avg_unct_quad_D_odd' :all_std_offset_quad_D_odd, 
                     'Avg_unct_quad_D_even' :all_std_offset_quad_D_even,
                     })
                
    dframe1.to_csv(file_path+'/'+'Signal_dependent_offset.csv')
        #dframe2.to_csv(r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test\Unct_Quad_A_flooded.csv')
if __name__ == "__main__":
    main()