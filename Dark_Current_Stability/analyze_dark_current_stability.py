# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""


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
   
    
def create_image(image_data, title, figure_name):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='bwr', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)
    #plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    plt.show()    
    plt.close('all')   
    
    
def create_hist(image, title, figure_name, COLOR) : 
    
    if np.array(image).ndim ==2:      
        
        nx_quad, ny_quad = image.shape
    else:
        nx_quad= 1        
        ny_quad = len(image)
        #print(ny_quad)
        #cc
    
    label = 'Mean = '+ str(round(np.mean(image), 2))             
    plt.figure(figsize=(8, 5))
    plt.hist(np.reshape(image, (nx_quad* ny_quad, 1)),10, facecolor=COLOR, label=label)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    #plt.xlim(-10, 10)
    #plt.ylim(0, 40000)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' Dark current (DN)  ', fontsize=12,
              fontweight="bold")
    plt.title(title)     
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")    
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
        all_med_quad_A_odd = [ ]
        all_med_quad_B_odd = [ ]
        all_med_quad_C_odd = [ ]
        all_med_quad_D_odd = [ ]
        all_med_quad_A_even = [ ]
        all_med_quad_B_even = [ ]
        all_med_quad_C_even = [ ]
        all_med_quad_D_even = [ ]
        dframe1 = pd.DataFrame()
        all_med_sto_A_odd = [ ]
        all_med_sto_B_odd = [ ]
        all_med_sto_C_odd = [ ]
        all_med_sto_D_odd = [ ]
        all_med_sto_A_even = [ ]
        all_med_sto_B_even = [ ]
        all_med_sto_C_even = [ ]
        all_med_sto_D_even = [ ]
        dframe2 = pd.DataFrame()
        int_time_all = []
        temp_all = []        
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
            int_time_all.append(int_time)            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q           
            for i in range(0, 4):
                quad_full_frame = all_full_frame[:, i, :, :]
                avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
                active_quad = avg_quad[4:1028, 10:1034]
                tsoc= avg_quad[4:1028, 1034:1056]
                tsoc_storage = avg_quad[0:2, 1034:1056]
                storage_quad = avg_quad[1, 10:1034]
                title='Quad A Image'
                figure_name = 'test.png'
                create_image(active_quad, title, figure_name)
                
                # let us examine the 
                
                #------perform bias subtraction using trailing overclocks and save the dark current image----------
                bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc)
                bias_subtracted_quad_storage = perform_bias_subtraction_ave_sto(storage_quad, tsoc_storage)
                                
                dark_current_rate = bias_subtracted_quad
                dark_current_rate_storage = bias_subtracted_quad_storage
                
                #dark_current_rate = dark_current_rate/np.max(dark_current_rate)
                
                if i == 0:
                    all_med_quad_A_even.append(np.median((dark_current_rate[:, ::2])))
                    all_med_quad_A_odd.append(np.median((dark_current_rate[:, 1::2])))
                    all_med_sto_A_even.append(np.median((dark_current_rate_storage[:, ::2])))
                    all_med_sto_A_odd.append(np.median((dark_current_rate_storage[:, 1::2])))
                    
                elif i == 1:
                    all_med_quad_B_even.append(np.median((dark_current_rate[:, ::2])))
                    all_med_quad_B_odd.append(np.median((dark_current_rate[:, 1::2])))
                    all_med_sto_B_even.append(np.median((dark_current_rate_storage[:, ::2])))
                    all_med_sto_B_odd.append(np.median((dark_current_rate_storage[:, 1::2])))
                elif i == 2:
                    all_med_quad_C_even.append(np.median((dark_current_rate[:, ::2])))
                    all_med_quad_C_odd.append(np.median((dark_current_rate[:, 1::2])))
                    all_med_sto_C_even.append(np.median((dark_current_rate_storage[:, ::2])))
                    all_med_sto_C_odd.append(np.median((dark_current_rate_storage[:, 1::2])))
                elif i == 3:                    
                    all_med_quad_D_even.append(np.median((dark_current_rate[:, ::2])))
                    all_med_quad_D_odd.append(np.median((dark_current_rate[:, 1::2])))
                    all_med_sto_D_even.append(np.median((dark_current_rate_storage[:, ::2])))
                    all_med_sto_D_odd.append(np.median((dark_current_rate_storage[:, 1::2])))
                  
                quad_save = 'Dark_Current_Image'
                save_dir_image = os.path.join(save_dir, quads[i],quad_save)
                if not os.path.exists(save_dir_image):
                    os.makedirs(save_dir_image)                 
                title = 'Dark Current Rate Image, '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
                figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'                  
                create_image(dark_current_rate, title, figure_name)
                
                # Hsotgram of dark current in active region
                hist_dir = 'Hist_Dark_Current_Active'
                diff_hist_dir = os.path.join(save_dir, quads[i], hist_dir)
                if not os.path.exists(diff_hist_dir):
                    os.makedirs(diff_hist_dir)
                COLOR = 'blue'
                figure_name = diff_hist_dir + '/'+ string1 + str(int_time) + '_hist_dark_current.png'
                title = 'Histogram of dark current (Active Region) \n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
                if int_time>0:
                    dark_current_rate = filter_outlier_median(dark_current_rate)
                create_hist(dark_current_rate, title, figure_name, COLOR)  
                 
                # Histogram of dark current in storage region
                hist_dir = 'Hist_Dark_Current_Storage'
                diff_hist_dir = os.path.join(save_dir, quads[i], hist_dir)
                if not os.path.exists(diff_hist_dir):
                    os.makedirs(diff_hist_dir)
                COLOR= 'red'
                figure_name = diff_hist_dir + '/'+ string1 + str(int_time) + '_hist_dark_current.png'
                title = 'Histogram of dark current (Storage) \n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
                if int_time>0:
                    dark_current_rate_storage = filter_outlier_median(dark_current_rate_storage)
                create_hist(filter_outlier_median(dark_current_rate_storage), title, figure_name, COLOR)
                
        
        dframe1 = pd.DataFrame(
                    {'Int_time.' : int_time_all,
                     'Quad_A_Even' : all_med_quad_A_even,
                     'Quad_A_Odd' : all_med_quad_A_odd,
                     'Quad_B_Even' : all_med_quad_B_even,
                     'Quad_B_Odd' : all_med_quad_B_odd,
                     'Quad_C_Even' : all_med_quad_C_even,
                     'Quad_C_Odd' : all_med_quad_C_odd,
                     'Quad_D_Even' : all_med_quad_D_even,
                     'Quad_D_Odd' : all_med_quad_D_odd
                     
                     })
    
    
        dframe2 = pd.DataFrame(
                    {'Int_time.' : int_time_all,
                     'sto_A_Even' : all_med_sto_A_even,
                     'sto_A_Odd' : all_med_sto_A_odd,
                     'sto_B_Even' : all_med_sto_B_even,
                     'sto_B_Odd' : all_med_sto_B_odd,
                     'sto_C_Even' : all_med_sto_C_even,
                     'sto_C_Odd' : all_med_sto_C_odd,
                     'sto_D_Even' : all_med_sto_D_even,
                     'sto_D_Odd' : all_med_sto_D_odd
                     
                     })
    
        #dframe1.to_csv(save_dir +'/'+'Active_region_Dark_Current_vs_Temp.csv')
        #dframe2.to_csv(save_dir +'/'+'Storage_region_Dark_Current_vs_Temp.csv')
        dframe1= None
        dframe2= None
        all_med_quad_A_even =  all_med_quad_A_odd= None
        all_med_quad_B_even =  all_med_quad_B_odd= None
        all_med_quad_C_even =  all_med_quad_C_odd= None
        all_med_quad_D_even =  all_med_quad_D_odd= None
        
        all_med_sto_A_even =  all_med_sto_A_odd= None
        all_med_sto_B_even =  all_med_sto_B_odd= None
        all_med_sto_C_even =  all_med_sto_C_odd= None
        all_med_sto_D_even =  all_med_sto_D_odd= None
         
         
        
         
          
if __name__ == "__main__":
    main()
    