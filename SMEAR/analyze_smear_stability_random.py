# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""
# This script uses the PRNU map derived from integration time sweep data to some
# of the images acquired during intesnity sweep.

from random import randint
import os
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo


def get_oultier_mask():
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask'
    file_name = 'final_outlier_mask.csv'
    outlier_mask = genfromtxt(file_path+'/'+file_name, delimiter=',')
    rows, cols = outlier_mask.shape
    mask_A = outlier_mask[0: rows/2, 0:cols/2]
    mask_B = outlier_mask[rows/2:rows, 0:cols/2]
    mask_C = outlier_mask[rows/2:rows, cols/2:cols]
    mask_D = outlier_mask[0:rows/2, cols/2:cols]
    outlier_mask = [mask_A, mask_B, mask_C, mask_D]
    return outlier_mask
    

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
    
def perform_linearity(active_quad, quad_name):
    
    # perform linearity based on look up table. 
        
    quad = quad_name.split(' ')[1]
    active_quad[active_quad<0] = 0
    
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file= 'Linearity_look_up_table_final_DN.csv'
    linearity_file = pd.read_csv(path+'/'+ linearity_file)
    # for even quad
    lut_even =  linearity_file[['Input Signal (DN)', quad+str(1)]].dropna()   
    lut_even['DN'] = lut_even.pop('Input Signal (DN)').astype(int)
    lut_even['Value'] = lut_even.pop(quad+str(1)).astype(float)
    lut_even = lut_even.set_index('DN') 
    even_detector_active_quad = active_quad[:, ::2]
    nx_half_quad, ny_half_quad = even_detector_active_quad.shape
    even_detector_active_quad = np.reshape(np.array(even_detector_active_quad), (nx_half_quad*ny_half_quad, 1))
    df = pd.DataFrame(data=even_detector_active_quad)
    df.columns = ['DN']
    datin_even = df['DN'].astype(int)
    linearized_quad_even = datin_even.apply(lambda x:lut_even.ix[x])
    even_pixels_quad = np.reshape(linearized_quad_even.values,(nx_half_quad, ny_half_quad))    
    df = None
    
    lut_odd =  linearity_file[['Input Signal (DN)', quad+str(2)]].dropna()   
    lut_odd['DN'] = lut_odd.pop('Input Signal (DN)').astype(int)
    lut_odd['Value'] = lut_odd.pop(quad+str(2)).astype(float)
    lut_odd = lut_odd.set_index('DN')    
    odd_detector_active_quad = active_quad[:, 1::2]
    
    nx_half_quad, ny_half_quad = odd_detector_active_quad.shape
    odd_detector_active_quad = np.reshape(np.array(odd_detector_active_quad), (nx_half_quad*ny_half_quad, 1))
    df = pd.DataFrame(data=odd_detector_active_quad)
    df.columns = ['DN']
    datin_odd = df['DN'].astype(int) 
    
    linearized_quad_odd = datin_odd.apply(lambda x: lut_odd.ix[x])
    odd_pixels_quad = np.reshape(linearized_quad_odd.values,(nx_half_quad, ny_half_quad))
    # now let's put quad in place
    
    nx_quad,ny_quad = active_quad.shape
    linearized_quad = np.array([[0]*ny_quad]*nx_quad)
    linearized_quad = np.reshape(linearized_quad, (nx_quad, ny_quad))
    linearized_quad[:, ::2] = even_pixels_quad 
    linearized_quad[:, 1::2] = odd_pixels_quad 
    return linearized_quad

def perform_smear_subtraction(active_quad, int_time):
    # the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis. 
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO 
    # tft = 8ms 
    tFT = 8*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)   
      
    smear_subtracted_quad = active_quad - smear_factor[None, :]
#    print(np.shape(smear_factor[None, :]))
#    print(np.shape(active_quad))
#    cc
    #column_average = np.mean(active_quad, axis=0)
    return smear_subtracted_quad, smear_factor

def smear_calculation_overclocks(overclocks, trailing_overclocks, active_quad):
    # Ok perform smear subtraction using parallel overclocks
    # seperate out the even and odd
    overclocks= overclocks[2:, :]
    nx_quad, ny_quad = overclocks.shape
    even_detector_smear_quad = overclocks[:, ::2]       
    odd_detector_smear_quad = overclocks[:, 1::2]
    even_detector_bias = trailing_overclocks[ :, ::2]
    odd_detector_bias = trailing_overclocks[:, 1::2]    
    smear_quad = np.array([[0]*ny_quad]*nx_quad)
    even_smear = even_detector_smear_quad- np.mean(filter_outlier_median(even_detector_bias))    
    odd_smear = odd_detector_smear_quad- np.mean(filter_outlier_median(odd_detector_bias))
    smear_quad = np.reshape(smear_quad,(nx_quad, ny_quad))
    smear_quad[:, ::2] = even_smear
    smear_quad[:, 1::2] = odd_smear     
    avg_smear = np.mean(smear_quad, axis=0)
    smear_subtracted_quad = active_quad - avg_smear[None, :] 
    
    return smear_subtracted_quad, avg_smear
       

def plot_smear_signal(smear, int_time_all, quad, color):
    plt.figure()
    for xe, ye in zip(int_time_all, smear):
      plt.scatter([xe] * len(ye), ye, marker='.')  
    
    title = 'Estimate SMEAR Vs. Integration time, '+ quad
    plt.title(title ,fontsize=12,fontweight="bold")
    r" $\mu$" +'secs'
    plt.xlabel('Integration time ('+r" $\mu$" +'secs)')
    plt.ylabel('Smear Signal Estimate (DN)')
    plt.grid(True, linestyle=':')
    plt.ylim(700, 940)
    plt.show()
    #plt.xlim(800, 16000)    
    #plt.savefig(figure_name,dpi=100,bbox_inches="tight")    
    plt.close('all')

def main():
    """
    Tme main function
    """
  
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    if 'Integration_Sweep' in file_path1:
    #if 'Bar_target' in file_path1:
        saturated_collects = [ 'FT6_LONG_INT_130018.dat.sav','FT6_SHORT_INT_0.dat.sav',
                              'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']

    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')] 
    nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files] # for cross talk image
    
  # # nominal_int_files = 'SPI_20160912100114585_ACLMP-0X2E-5_DCLMP-0X2E-5_ASMPL-0X8-7_DSMPL-0X8-7_FT6_SHORT_INT_15009.dat.sav'                     
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Smear\Using_integ_sweep_data'
    smear_A = []
    smear_B = []
    smear_C = []
    smear_D = []
    int_time_all = [ ]
    
    
    for i in range(0, 4): # for the 4 quads
                
        for data_files in nominal_int_files:
            print(data_files)
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)          
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
               
            int_time_all.append(int_time)
            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data

                   
            all_full_frame = IDL_variable.q 
            print(all_full_frame.shape)
            cc
            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)       
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)     
 
            #------perform bias subtraction using trailing overclocks----------
            active_quad_A = perform_bias_subtraction_ave(quad_A[4:1028, 10:1034], tsoc_A)
            
            # perform non-lineairty correction using look up table
            #linearized_quad_A = perform_linearity(active_quad_A, quads[i])
            
            # Skip lineairty if you are using a bright homegeneous photon transfer data
            linearized_quad_A = active_quad_A     
          
            # Now the method of active quad            
            active_quad_A_corr, smear_signal_active_quad = perform_smear_subtraction(linearized_quad_A, int_time)
           
            if i==0:
                smear_A.append(smear_signal_active_quad)
            elif i==1:
                smear_B.append(smear_signal_active_quad)
            elif i==2:
                 smear_C.append(smear_signal_active_quad)
            elif i==3:
                 smear_D.append(smear_signal_active_quad)
         
        
        if not os.path.exists(save_dir):
               os.makedirs(save_dir)
        csv_name_A = save_dir + '/' + '_SMEAR_A.csv'
        csv_name_B = save_dir + '/' + '_SMEAR_B.csv'
        csv_name_C = save_dir + '/' + '_SMEAR_C.csv'
        csv_name_D = save_dir + '/' + '_SMEAR_D.csv'           
        int_name_all  = save_dir + '/'+ 'all_int_time.csv'
      
        np.savetxt(csv_name_A, np.array(smear_A), delimiter=",", fmt='%1.2f')
        np.savetxt(csv_name_B, np.array(smear_B), delimiter=",", fmt='%1.2f')
        np.savetxt(csv_name_C, np.array(smear_C), delimiter=",", fmt='%1.2f')
        np.savetxt(csv_name_D, np.array(smear_D), delimiter=",", fmt='%1.2f')
        np.savetxt(int_name_all, np.array( int_time_all), delimiter=",", fmt='%1.2f')    
       
#        quad = ['Quad A' ,'Quad B', 'Quad C' , 'Quad D']
#        colors = ['blue','green','red','orange']
#        plot_smear_signal(smear_A, int_time_all, quad[0], colors[0])
#        plot_smear_signal(smear_B, int_time_all, quad[1], colors[1]) 
#        plot_smear_signal(smear_C, int_time_all, quad[2], colors[2])  
#        plot_smear_signal(smear_D, int_time_all, quad[3], colors[3])         
if __name__ == "__main__":
    main()
    