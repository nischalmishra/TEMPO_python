# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 08:16:04 2017

@author: nmishra
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""

import os
import numpy as np
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
    elif np.array(quads).ndim ==2:      
        ndims=1
        nx_quad, ny_quad = quads.shape
    else:
        nx_quad= 1
        ndims=1
        ny_quad = len(quads)
        
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    #print(outlier_filtered_data)
    return outlier_filtered_data

def perform_bias_subtraction_ave (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, ::2]
    # remove outliers 
    # First 4 hot lines in even and odd
    # last odd lne in odd
    even_detector_bias = even_detector_bias[:, 4:]
    avg_bias_even = np.mean(even_detector_bias, axis=1)  
    odd_detector_bias = trailing_overclocks[:, 1::2]
    odd_samples = odd_detector_bias[:, 4:]
    rows, cols = odd_samples.shape   
    odd_detector_bias = odd_samples[:, 0:cols-1] 
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]     
    odd_detector_active_quad = active_quad[:, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (nx_quad, ny_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad
    
 
def plot_few_storage_lines(storage_lines, figure_name, title):
    #print(storage_lines.shape)
    line1 = storage_lines[0,:].T
    line2 = storage_lines[1,:].T
    plt.plot(filter_outlier_median(line1),'.', color='blue', label='Line 1')
    plt.plot(filter_outlier_median(line2),'.', color='red', label='Line2')
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.grid(True, linestyle=':')
    plt.title(title, fontsize=12, fontweight= 'bold')
    plt.xlabel('Spatial Pixel Indices(#)', fontsize=12, fontweight= 'bold')
    #plt.ylim(700, 1000)
    plt.ylabel('Raw Signal (DN)', fontsize=12, fontweight='bold')        
    plt.savefig(figure_name, dpi=500, bbox_inches="tight")
    plt.show()
    cc
    plt.close('all')
    
def plot_few_tsoc_lines(tsoc_lines, figure_name, title):
    #print(storage_lines.shape)
    line1 = tsoc_lines[0,:].T
    line2 = tsoc_lines[1,:].T
    plt.plot((line1), color='blue', label='Line 1')
    plt.plot((line2), color='red', label='Line2')
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.grid(True, linestyle=':')
    plt.title(title, fontsize=12, fontweight= 'bold')
    plt.xlabel('Spatial Pixel Indices(#)', fontsize=12, fontweight= 'bold')
    #plt.ylim(700, 1000)
    plt.ylabel('Raw Signal (DN)', fontsize=12, fontweight='bold')        
    plt.savefig(figure_name, dpi=500, bbox_inches="tight")
    plt.close('all')
    
def create_hist(image, title, figure_name) : 
    # hisrogram of line 2
    image = image[1,:]
    even_pixels = filter_outlier_median(image[::2])
    odd_pixels = filter_outlier_median(image[1::2])        
    label1 = 'Mean (Even Pixels) = '+ str(round(np.mean(even_pixels), 2))+ 'DN' 
    label2 = 'Mean (Odd Pixels) = '+ str(round(np.mean(odd_pixels), 2))+  'DN'           
    plt.figure(figsize=(8, 5))
    plt.hist(even_pixels, 40, facecolor='green', label=label1)
    plt.hist(odd_pixels, 40, facecolor='orange', label=label2)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    #plt.xlim(700, 1000)
    #plt.ylim(0, 300)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' Raw Signal (DN) ', fontsize=12,
              fontweight="bold")
    plt.title(title,fontsize=12,
              fontweight="bold")    
    plt.savefig(figure_name, dpi=500, bbox_inches="tight")
    plt.show()    
    plt.close('all')
    
    
    
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    save_file_path = r'C:\Users\nmishra\Workspace\TEMPO\FPA_Dark_Current_Stability'    
    temp_files = os.listdir(file_path) 
    for files in range(1, len(temp_files)):              
        save_dir =  os.path.join(save_file_path,temp_files[files])         
        if not os.path.exists(save_dir):
               os.makedirs(save_dir)
        saved_data_files = os.path.join(file_path, temp_files[files])         
        
        all_int_files = [each for each in os.listdir(saved_data_files) \
                     if each.endswith('.dat.sav')]          
        for data_files in all_int_files:
            print(data_files)            
            data_file = os.path.join(saved_data_files, data_files)           
            IDL_variable = readsav(data_file)
            data_path_name_split = data_files.split('_')
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
           # temp = int(data_path_name_split[-4][0:3])
            #temp_all.append(temp)
            string1 = 'Integ_time_'
            string2 = 'Int.time = '                    
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q           
            for i in range(0, 4):
                quad_full_frame = all_full_frame[:, i, :, :]
                avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
                tsoc_storage = avg_quad[0:2, 1034:1056]
                storage_quads = avg_quad[0:2, 10:1034]
                # let us examine the storage region
                storage_region_plot = 'Storage_region_profile'
                save_dir_image = os.path.join(save_dir, quads[i],storage_region_plot)
                if not os.path.exists(save_dir_image):
                    os.makedirs(save_dir_image)
                title = 'Storage Region Profile, '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
                figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
                plot_few_storage_lines(storage_quads, figure_name, title)
                
                # Let us examine the tsocs of the storage region
                storage_region_tsoc_plot = 'Storage_region_tsoc'
                save_dir_image = os.path.join(save_dir, quads[i], storage_region_tsoc_plot)
                if not os.path.exists(save_dir_image):
                    os.makedirs(save_dir_image)
                title = 'Trailing Serial Overclock Profile (Storage Region),\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
                figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
                plot_few_tsoc_lines(tsoc_storage, figure_name, title)
            
            
            
#                hist_dir = 'Hist_storage_region'
#                diff_hist_dir = os.path.join(save_dir, quads[i], hist_dir)
#                if not os.path.exists(diff_hist_dir):
#                    os.makedirs(diff_hist_dir)
#                figure_name = diff_hist_dir + '/'+ string1 + str(int_time) + '_hist_storage_region.png'
#                title = 'Histogram of Storage Region (Second Line) \n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
#                create_hist(storage_quads, title, figure_name)  
      
if __name__ == "__main__":
    main()
    