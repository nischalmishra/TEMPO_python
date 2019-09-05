# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def read_outlier_mask():
    outlier_mask= np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask\final_outlier_mask_2_sigma.csv', delimiter=',')
    quad_A = outlier_mask[0:1024, 0:1024]
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

def extract_tsoc(raw_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads_odd = []
    all_quads_even = []
    num_quads = raw_quads.shape[0]
    
    #num_quads=2
    for quads in range(0, num_quads):        
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]        
        odd_detector_bias = trailing_overclocks[ :, ::2]
        
                
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        odd_detector_bias = odd_detector_bias[:, 4:]  
        
        even_detector_bias = trailing_overclocks[:, 1::2]
        even_detector_bias= even_detector_bias[:, 4:-1]        
        all_quads_even.append(np.mean(even_detector_bias, axis=1))        
        all_quads_odd.append(np.mean(odd_detector_bias, axis=1))
       
        
    
    return np.array(all_quads_even), np.array(all_quads_odd)


    
 
def plot_row_avg(row_avg, title, figure_name, COLOR,  xlabel):
    # let's take the mean tsoc for 100 frames
    
    nrows = 1
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(7,5))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    
    ax.plot(row_avg, '.', color=COLOR)
    ax.grid(True, linestyle=':')
    ax.set_title(title, fontsize=12, fontweight='bold')    
    ax.set_ylabel('Signal - Offset (DN)', fontsize=12, fontweight='bold')   
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    #ax.set_ylim(ylim[0], ylim[1])
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight") 
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
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Crosstalk'
    save_file_path = r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test_After_TVAC' 
    temp_files = os.listdir(file_path) 
    
    for files in range(0, 4):   
        dframe1 = []
        dframe2 = []
        rows_max_A = [ ]
        tsoc_A_even = [ ] 
        tsoc_A_odd = [ ] 
        rows_max_B = [ ]
        cols_max_B = [ ]
        rows_max_C = [ ]
        cols_max_C = [ ]
        rows_max_D = [ ]
        cols_max_D = [ ]
        tsoc_B_even = [ ] 
        tsoc_B_odd = [ ]
        tsoc_C_even = [ ] 
        tsoc_C_odd = [ ]
        tsoc_D_even = [ ] 
        tsoc_D_odd = [ ]
        
        save_dir =  os.path.join(save_file_path, temp_files[files])             
        if not os.path.exists(save_dir):
               os.makedirs(save_dir)
        saved_data_files = os.path.join(file_path, temp_files[files],'Script_Data','saved_quads') 
        
        all_int_files = [each for each in os.listdir(saved_data_files) \
                         if each.endswith('dat.sav')] 
             
        for data_files in all_int_files:
            data_file = os.path.join(saved_data_files, data_files) 
            print(data_file)
            
            IDL_variable = readsav(data_file)
            data_path_name_split = data_files.split('_')
            int_time = int(data_path_name_split[-1].split('.')[0])        
            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            color = ['Blue','Green','Red','Orange']
            all_full_frame = IDL_variable.q
            all_full_frame = np.mean(all_full_frame, axis=0)
            even_tsoc, odd_tsoc = extract_tsoc(all_full_frame)
                    
            ylim1= [0, 16000]
            ylim2 = [-6, 6]
            for i in range(0, 4):
                quad_full_frame = all_full_frame[i, :, :]             
                active_quad = quad_full_frame[2:1028, 10:1034]               
             
                cross_talk_array = active_quad
                nx1, ny1 = cross_talk_array.shape
                # let's reshape the array to 1-D so we can work with single loop
                #cross_talk_array = np.reshape(cross_talk_array, (nx1*ny1, 1))
                
                
              
                if(temp_files[files] in("Channel_A", "Channel_C")) :
                    
                    if len(data_path_name_split)>9:
                        spot = 2
                        input_signal = (data_path_name_split[-5])
                        quad_illuminated = data_path_name_split[-6]
                    else:
                        spot = 1
                        input_signal = (data_path_name_split[-4])
                        quad_illuminated = data_path_name_split[-5]
           
                    # subtract off the dark current
                    
                    cross_talk_array = cross_talk_array 
                    row_average = np.mean(cross_talk_array, axis=1)                
                    column_average = np.mean(cross_talk_array, axis=0)
                    
                    string1 = quad_illuminated[0:4]+' '+quad_illuminated[4]+ ' Illuminated'    
                    string2 = 'Input Signal = '+ input_signal
                    string3 = 'spot'+ str(spot)                    
                    title1 = quads[i]+' Image\n ' + string1+' @'+string3+', '+string2
                    title2 = quads[i]+' Row Average Profile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                    title3 = quads[i]+' Column Average Profile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                    title4 = quads[i]+' Image Profile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                elif(temp_files[files] in("Channel_D")) :
                    
                    if len(data_path_name_split)>9:
                        spot = 2
                        input_signal = (data_path_name_split[-5])
                        quad_illuminated = data_path_name_split[-6]
                    else:
                        spot = 1                        
                        quad_illuminated = data_path_name_split[-4]
                        input_signal ='Not Given'
           
                    string1 = quad_illuminated[0:4]+' '+quad_illuminated[4]+ ' Illuminated'    
                    string2 = 'Input Signal = '+ input_signal
                    string3 = 'spot'+ str(spot)
#                    print(string1)
#                    print(string2)
#                    print(string3)
                    
                    title1 = quads[i]+' Image\n ' + '('+ string1+' @'+string3+', '+string2+')'
                    title2 = quads[i]+' Row Average Profile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                    title3 = quads[i]+' Column Average Profile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                    title4 = quads[i]+' ImageProfile \n ' +'('+ string1+' @'+string3+', '+string2 +')'
                
                else:
                    
                    if len(data_path_name_split)>8:
                        spot = 2
                        quad_illuminated = data_path_name_split[-5]
                    else:
                        spot=1
                        quad_illuminated = data_path_name_split[-4]
                    
                    string1 = quad_illuminated[0:4]+' '+quad_illuminated[4]+ ' Illuminated' 
                    string3 = 'Spot'+ str(spot)                            
                    title1 = quads[i]+' Image\n ' + string1+' @'+string3
                    title2 = quads[i]+' Row Average Profile \n ' + string1+' @'+string3
                    title3 = quads[i]+' Column Average Profile \n ' + string1+' @'+string3 
                    title4 = quads[i]+' Image Profile \n ' + string1+' @'+string3
                if quad_illuminated.lower() == quads[i].replace(" ","").lower():
                    ylim = ylim1
                    #print(ylim)
                else:
                    cross_talk_array[(cross_talk_array>1500)] = np.mean(cross_talk_array)
                    ylim = ylim2
                
                if spot == 1:
                    #rows, cols = np.reshape()
                    if i == 0: 
                        
                        rows_max_A.append(cross_talk_array[185:325, 170:310])
                        even_tsoc = even_tsoc[0, 185:325]
                        odd_tsoc = odd_tsoc[0, 185:325]
                        tsoc_A_even.append(even_tsoc)                        
                        tsoc_A_odd.append(odd_tsoc)
                                
                        
                    elif i == 1: 
                        
                        rows_max_B.append(cross_talk_array[185:325, 170:310])
                        even_tsoc = even_tsoc[1, 185:325]
                        odd_tsoc = odd_tsoc[1, 185:325]
                        tsoc_B_even.append(even_tsoc)
                        tsoc_B_odd.append(odd_tsoc)
                        
                    elif i == 2:                        
                        rows_max_C.append(cross_talk_array[185:325, 170:310])
                        even_tsoc = even_tsoc[2, 185:325]
                        odd_tsoc = odd_tsoc[2, 185:325]
                        tsoc_C_even.append(even_tsoc)
                        tsoc_C_odd.append(odd_tsoc)
                       
                    elif i == 3:                        
                        rows_max_D.append(cross_talk_array[185:325, 170:310])
                        even_tsoc = even_tsoc[3, 185:325]
                        odd_tsoc = odd_tsoc[3, 185:325]
                        tsoc_D_even.append(even_tsoc)
                        tsoc_D_odd.append(odd_tsoc)
                          
                
                elif spot==2:
                    if i==0:
                        rows_max_A.append(cross_talk_array[720:860, 720:860])
                        even_tsoc = even_tsoc[0, 720:860]
                        odd_tsoc = odd_tsoc[0, 720:860]
                        tsoc_A_even.append(even_tsoc)
                        tsoc_A_odd.append(odd_tsoc)
                       
                    elif i==1:
                        rows_max_B.append(cross_talk_array[720:860, 720:860])
                        even_tsoc = even_tsoc[1, 720:860]
                        odd_tsoc = odd_tsoc[1, 720:860]
                        tsoc_B_even.append(even_tsoc)
                        tsoc_B_odd.append(odd_tsoc)
                       
                    elif i==2:
                        rows_max_C.append(cross_talk_array[720:860, 720:860])
                        even_tsoc = even_tsoc[2, 720:860]
                        odd_tsoc = odd_tsoc[2, 720]
                        tsoc_C_even.append(even_tsoc)
                        tsoc_C_odd.append(odd_tsoc)
                        
                    elif i==3:
                        rows_max_D.append(cross_talk_array[720:860, 720:860])
                        even_tsoc = even_tsoc[3, 720:860]
                        odd_tsoc = odd_tsoc[2, 720:860]
                        tsoc_D_even.append(even_tsoc)
                        tsoc_D_odd.append(odd_tsoc)
                      
                dframe1['quadA'] = tsoc_A_even
                dframe1['quadB'] = tsoc_B_even
                dframe1['quadC'] = tsoc_C_even
                dframe1['quadD'] = tsoc_D_even               
                quad_save = 'Cross_Talk_tsoc_Ghost'
                save_dir_tsoc = os.path.join(save_dir, quad_save)
                if not os.path.exists(save_dir_tsoc):
                    os.makedirs(save_dir_tsoc)             
              

        ndims, row_s,col_s = np.array(rows_max_A).shape        
        row_s_tsoc, col_s_tsoc = np.array(tsoc_A_odd).shape 
        #print(np.array(tsoc_A_odd).shape)
        #cc
        
        rows_max_A = np.reshape(np.array(rows_max_A), (ndims*row_s*col_s, 1))
        rows_max_B = np.reshape(np.array(rows_max_B), (ndims* row_s*col_s,1 ))
        rows_max_C = np.reshape(np.array(rows_max_C), (ndims*row_s*col_s, 1))
        rows_max_D = np.reshape(np.array(rows_max_D), (ndims*row_s*col_s, 1))
        
        tsoc_A_odd = np.reshape(np.array(tsoc_A_odd), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_B_odd = np.reshape(np.array(tsoc_B_odd), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_C_odd = np.reshape(np.array(tsoc_C_odd), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_D_odd = np.reshape(np.array(tsoc_D_odd), (row_s_tsoc*col_s_tsoc, 1))
        
        tsoc_A_even = np.reshape(np.array(tsoc_A_even), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_B_even = np.reshape(np.array(tsoc_B_even), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_C_even = np.reshape(np.array(tsoc_C_even), (row_s_tsoc*col_s_tsoc, 1))
        tsoc_D_even = np.reshape(np.array(tsoc_D_even), (row_s_tsoc*col_s_tsoc, 1))
        
        
        csv_name_A = save_dir+'/'+temp_files[files]+'_cross_talk_A.csv'        
        csv_name_B = save_dir+'/'+temp_files[files]+'_cross_talk_B.csv'
        csv_name_C = save_dir+'/'+temp_files[files]+'_cross_talk_C.csv'
        csv_name_D = save_dir+'/'+temp_files[files]+'_cross_talk_D.csv'
        
        csv_name_tsoc_A_odd = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_A_odd.csv'        
        csv_name_tsoc_B_odd = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_B_odd.csv'
        csv_name_tsoc_C_odd = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_C_odd.csv'
        csv_name_tsoc_D_odd = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_D_odd.csv'
        
                                                     
        csv_name_tsoc_A_even = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_A_even.csv'        
        csv_name_tsoc_B_even = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_B_even.csv'
        csv_name_tsoc_C_even = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_C_even.csv'
        csv_name_tsoc_D_even = save_dir+'/'+temp_files[files]+'_cross_talk_tsoc_D_even.csv'
        
        np.savetxt(csv_name_A, np.asarray(rows_max_A), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_B, np.asarray(rows_max_B), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_C, np.asarray(rows_max_C), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_D, np.asarray(rows_max_D), delimiter=',', fmt='%1.2f')
        
        np.savetxt(csv_name_tsoc_A_odd, np.asarray(tsoc_A_odd), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_B_odd, np.asarray(tsoc_B_odd), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_C_odd, np.asarray(tsoc_C_odd), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_D_odd, np.asarray(tsoc_D_odd), delimiter=',', fmt='%1.2f')
        
        np.savetxt(csv_name_tsoc_A_even, np.asarray(tsoc_A_even), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_B_even, np.asarray(tsoc_B_even), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_C_even, np.asarray(tsoc_C_even), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_tsoc_D_even, np.asarray(tsoc_D_even), delimiter=',', fmt='%1.2f')
        
        
       
        #cc
if __name__ == "__main__":
    main()
    