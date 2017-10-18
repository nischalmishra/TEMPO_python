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
    
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads'
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
            int_time = round(int(data_path_name_split[-1].split('.')[0])) #for integ_sweep             
            data_file = os.path.join(file_path, data_files)  
            IDL_variable = readsav(data_file)            
            all_full_frame = IDL_variable.q  
            all_int_time.append(int_time)
            
            quad_A = all_full_frame[:, 0, :, :]
            active_quad_A = np.mean(quad_A[:, 4:1028, 10:1034], axis=0)
            bias_A = np.mean(quad_A[:, 4:1028, 1034:1056], axis=0)
            avg_bias_A_even  = np.mean(filter_outlier_median(bias_A[:, ::2])) 
            avg_bias_A_odd  = np.mean(filter_outlier_median(bias_A[:, 1::2])) 
            avg_quad_A_even = np.mean(filter_outlier_median(active_quad_A[:, ::2]))  
            avg_quad_A_odd = np.mean(filter_outlier_median(active_quad_A[:, 1::2]))
            #print(np.mean(filter_outlier_median(bias_subtracted_quad_A)))
            #print(avg_quad_A_odd-avg_bias_A_odd)
            #print(avg_quad_A_even-avg_bias_A_even)
            bias_subtracted_quad = (np.mean(filter_outlier_median(perform_bias_subtraction(quad_A[:, 4:1028, 10:1034], quad_A[:, 4:1028, 1034:1056]))))
            all_med_quad_A_odd.append(bias_subtracted_quad)
            all_avg_tsoc_quad_A_odd.append(avg_bias_A_odd)
            all_std_offset_quad_A_odd.append(np.std(filter_outlier_median(bias_A[:, 1::2])))
            all_med_quad_A_even.append(bias_subtracted_quad)
            all_avg_tsoc_quad_A_even.append(avg_bias_A_even)
            all_std_offset_quad_A_even.append(np.std(filter_outlier_median(bias_A[:, ::2])))
            bias_subtracted_quad = None   
            quad = 'Quad A'
            #plot_few_tsocs(quad_A[:, 4:1028, 1034:1056], plot_dir, quad, int_time)
            
            
            # For quad B
            quad_B = all_full_frame[:, 1, :, :]
            active_quad_B = np.mean(quad_B[:, 4:1028, 10:1034], axis=0) 
            bias_B = np.mean(quad_B[:, 4:1028, 1034:1056], axis=0) 
            avg_bias_B_even  = np.mean(filter_outlier_median(bias_B[:, ::2])) 
            avg_bias_B_odd  = np.mean(filter_outlier_median(bias_B[:, 1::2]))   
            avg_quad_B_even = np.mean(filter_outlier_median(active_quad_B[:, ::2]))  
            avg_quad_B_odd = np.mean(filter_outlier_median(active_quad_B[:, 1::2]))
            bias_subtracted_quad = (np.mean(filter_outlier_median(perform_bias_subtraction(quad_B[:, 4:1028, 10:1034], quad_B[:, 4:1028, 1034:1056]))))
            all_med_quad_B_odd.append(bias_subtracted_quad)
            all_avg_tsoc_quad_B_odd.append(avg_bias_B_odd)
            all_std_offset_quad_B_odd.append(np.std(filter_outlier_median(bias_B[:, 1::2])))
            all_med_quad_B_even.append(bias_subtracted_quad)
            all_avg_tsoc_quad_B_even.append(avg_bias_B_even)
            all_std_offset_quad_B_even.append(np.std(filter_outlier_median(bias_B[:, ::2])))  
            quad = 'Quad B'
           # plot_few_tsocs(quad_B[:, 4:1028, 1034:1056], plot_dir, quad, int_time)
            bias_subtracted_quad = None
            # For quad C
            quad_C = all_full_frame[:,2,:,:]
            active_quad_C = np.mean(quad_C[:, 4:1028, 10:1034], axis=0) 
            bias_C = np.mean(quad_C[:, 4:1028, 1034:1056], axis=0) 
            avg_bias_C_even  = np.mean(filter_outlier_median(bias_C[:, ::2])) 
            avg_bias_C_odd  = np.mean(filter_outlier_median(bias_C[:, 1::2]))   
            avg_quad_C_even = np.mean(filter_outlier_median(active_quad_C[:, ::2]))  
            avg_quad_C_odd = np.mean(filter_outlier_median(active_quad_C[:, 1::2]))
            bias_subtracted_quad = (np.mean(filter_outlier_median(perform_bias_subtraction(quad_C[:, 4:1028, 10:1034], quad_C[:, 4:1028, 1034:1056]))))
            
            
            all_med_quad_C_odd.append(bias_subtracted_quad)
            all_avg_tsoc_quad_C_odd.append(avg_bias_C_odd)
            all_std_offset_quad_C_odd.append(np.std(filter_outlier_median(bias_C[:, 1::2])))
            all_med_quad_C_even.append(bias_subtracted_quad)
            all_avg_tsoc_quad_C_even.append(avg_bias_C_even)
            all_std_offset_quad_C_even.append(np.std(filter_outlier_median(bias_C[:, ::2])))                  
            quad = 'Quad C'
            #plot_few_tsocs(quad_C[:, 4:1028, 1034:1056], plot_dir, quad, int_time)
            
            quad_D = all_full_frame[:,3,:,:]
            active_quad_D = np.mean(quad_D[:, 4:1028, 10:1034], axis=0) 
            bias_D = np.mean(quad_D[:, 4:1028, 1034:1056], axis=0) 
            avg_bias_D_even  = np.mean(filter_outlier_median(bias_D[:, ::2])) 
            avg_bias_D_odd  = np.mean(filter_outlier_median(bias_D[:, 1::2]))   
            avg_quad_D_even = np.mean(filter_outlier_median(active_quad_D[:, ::2]))  
            avg_quad_D_odd = np.mean(filter_outlier_median(active_quad_D[:, 1::2]))
            bias_subtracted_quad = (np.mean(filter_outlier_median(perform_bias_subtraction(quad_D[:, 4:1028, 10:1034], quad_D[:, 4:1028, 1034:1056]))))

            all_med_quad_D_odd.append( bias_subtracted_quad)
            all_avg_tsoc_quad_D_odd.append(avg_bias_D_odd)
            all_std_offset_quad_D_odd.append(np.std(filter_outlier_median(bias_D[:, 1::2])))
            all_med_quad_D_even.append( bias_subtracted_quad)
            all_avg_tsoc_quad_D_even.append(avg_bias_D_even)
            all_std_offset_quad_D_even.append(np.std(filter_outlier_median(bias_D[:, ::2])))                
            quad = 'Quad D'
            #plot_few_tsocs(quad_D[:, 4:1028, 1034:1056], plot_dir, quad, int_time)
            bias_subtracted_quad = None
            # let's save quad images too
             
#            quad_a = np.mean(quad_A, axis=0)            
#            quad_b = np.fliplr(np.mean(quad_B, axis=0))
#            quad_c = np.rot90(np.mean(quad_C, axis=0), 2)
#            #print(quad_c.shape)
#            quad_d = np.flipud(np.mean(quad_D, axis=0))            
#            lower_quads = np.concatenate((quad_a, quad_b), axis=1)
#            upper_quads = np.concatenate((quad_d, quad_c), axis=1)
#            quads = np.concatenate((lower_quads, upper_quads), axis=0)
#            plt.figure()
#            image = plt.imshow(quads, cmap='nipy_spectral', origin='lower', interpolation='none')
#            plt.grid(b=False)
#            cbar = plt.colorbar(image)
#            title = 'Full Frame Quad Image' + ' (Int.time = '+ str(int_time)+')' 
#            plt.title(title, fontsize=14)
#            plt.xlabel('# of spatial pixels', fontsize=12)
#            plt.ylabel('# of spectral pixels', fontsize=12)
#            plt.grid(False)   
#            plt.savefig(plot_dir+'/'+ data_files+'.png',dpi=100,bbox_inches="tight")
#            plt.close('all')
#       
#    print('A:', all_med_quad_A)
#    print('B:', all_med_quad_B)
#    print('C:',  all_med_quad_C)
#    print('D:', all_med_quad_D)
         
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