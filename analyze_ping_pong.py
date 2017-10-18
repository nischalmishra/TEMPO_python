# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 14:08:35 2017

@author: nmishra
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def read_outlier_mask():
    outlier_mask= np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_outlier_mask_2_sigma.csv', delimiter=',')
    quad_A = outlier_mask[0:1024, 0:1024]
    #print(quad_A[20])
    quad_B = outlier_mask[1024:, 0:1024]
    quad_C = outlier_mask[1024:, 1024:]
    quad_D = outlier_mask[0:1024:, 1024:]
    outlier_mask_final = [quad_A, quad_B, quad_C, quad_D]
    return outlier_mask_final
    
    """
    This function reads the outlier_mask
    """
    

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
    outlier_filtered_data = hist_data[measured_threshold < 6.]
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
    
def perform_bias_subtraction_ave_sto (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
   
    
    bias_subtracted_quad = np.zeros((1,1024))
    even_detector_bias = trailing_overclocks[:, ::2]
    even_detector_bias = even_detector_bias[:, 1:]
    avg_bias_even = np.mean(even_detector_bias)  
   #print(np.mean(avg_bias_even))
    odd_detector_bias = trailing_overclocks[:, 1::2] 
    odd_detector_bias = odd_detector_bias[:, 1:10 ]
    avg_bias_odd = np.mean(odd_detector_bias)
#    plt.plot(np.mean(even_detector_bias, axis=0).T,'.', color='blue')
#    plt.plot(np.mean(odd_detector_bias, axis=0).T,'.', color='black')
#    plt.show()
#    cc
#    
    #print(np.mean(avg_bias_odd))   
    even_detector_active_quad = active_quad[::2]     
    odd_detector_active_quad = active_quad[1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd     
    bias_subtracted_quad[:, ::2] = np.array(bias_subtracted_quad_even)
    bias_subtracted_quad[:, 1::2] = np.array(bias_subtracted_quad_odd)
    #print(avg_bias_even, avg_bias_odd, np.mean(bias_subtracted_quad_even), np.mean(bias_subtracted_quad_odd))
    return bias_subtracted_quad       
   
def perform_smear_subtraction(active_quad, int_time):
    # the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    tFT = 8*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)
    #print(smear_factor.shape)
    #cc
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad 
    

def create_image(image_data, title, figure_name, spot):
    plt.figure()
    ax = plt.gca()
    if spot==2:
    
        image = ax.imshow(image_data[720:860, 720:860], cmap='nipy_spectral', origin='lower')
    elif spot==1:
        image = ax.imshow(image_data[185:325, 170:310], cmap='nipy_spectral', origin='lower')
    
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)
    plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    #plt.show()    
    plt.close('all')   
    
    
def create_hist(image,outlier_mask, title, figure_name, COLOR) : 
    
    if np.array(image).ndim ==2:      
        
        nx_quad, ny_quad = image.shape
    else:
        nx_quad= 1        
        ny_quad = len(image)
        #print(ny_quad)
        #cc
    a = np.ma.array(image, mask = outlier_mask)
    label = 'Mean = '+ str(round(image.mean(), 2))             
    plt.figure(figsize=(8, 5))
    plt.hist(np.reshape(image, (nx_quad* ny_quad, 1)),50, facecolor=COLOR, label=label)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    #plt.xlim(-10, 10)
    #plt.ylim(0, 40000)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' Signal-Offset (DN)  ', fontsize=12,
              fontweight="bold")
    plt.title(title)     
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.show()
    cc
    plt.close('all')
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Ping_Pong'     
    outlier_mask = read_outlier_mask()     
    all_int_files = [each for each in os.listdir(file_path) \
                     if each.endswith('.dat.sav')]
    
    for files in range(0, 4):       
  
       
        for data_files in all_int_files:
            data_file = os.path.join(file_path, data_files) 
            print(data_files)
            IDL_variable = readsav(data_file)
            data_path_name_split = data_files.split('_')
            int_time = int(data_path_name_split[-1].split('.')[0])
            string1 = 'Integ_time_'
            string2 = 'Int.time = '  
            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            color = ['Blue','Green','Red','Orange']
            all_full_frame = IDL_variable.q          
            for i in range(0, 4):
                quad_full_frame = all_full_frame[:, i, :, :]
                avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
                active_quad = avg_quad[4:1028, 10:1034]               
                tsoc = avg_quad[4:1028, 1034:1056] 
                
                bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc)
                bias_subtracted_quad = np.array(bias_subtracted_quad)              
                
                outlier_mask = outlier_mask[i]
                outlier_mask = (outlier_mask==1).astype(int)
                outlier_masked = ma.masked_values(bias_subtracted_quad, outlier_mask)                             
                print((outlier_masked[~outlier_masked.mask]))
                #print(bias_subtracted_quad.view(ma.MaskedArray))
                #print(np.mean(bias_subtracted_quad))
                
             
                hist_save = 'saved_hist'
                save_hist_plot = os.path.join(save_dir, quads[i], hist_save)
                if not os.path.exists( save_hist_plot):
                    os.makedirs( save_hist_plot)
                figure_name= save_hist_plot + '/'+ string1 + str(int_time) + '_image.png'
                title ='Histogram of '+ quads[0] + ', ' + string2 + str(int_time)+ r" $\mu$" +'secs'                                             
                create_hist(bias_subtracted_quad, outlier_mask[i],title, figure_name, color[i])
                
                plt.plot(bias_subtracted_quad[400,:],'b')
                plt.show()
                cc

if __name__ == "__main__":
    main()