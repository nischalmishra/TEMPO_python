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
    tFT = 8.3333*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)
    #print(smear_factor.shape)
    #cc
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad 
    
def perform_Dark_removal(data_file, i):
    # calculate dark current
     IDL_variable = readsav(data_file)
     all_full_frame = IDL_variable.q
     quad_full_frame = all_full_frame[:, i , :, :]
     avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
     active_quad = avg_quad[4:1028, 10:1034]
     tsoc = avg_quad[4:1028, 1034:1056] 
     dark_current = perform_bias_subtraction_ave(active_quad, tsoc)
     return dark_current
    
    
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
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show()
    plt.close('all')
    
 
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
    file_path_dark = r'F:\TEMPO\Data\GroundTest\FPS\Crosstalk'
    save_file_path = r'C:\Users\nmishra\Workspace\TEMPO\Cross_Talk_Test' 
    outlier_mask = read_outlier_mask()   
    temp_files = os.listdir(file_path) 
    
    for files in range(0, 4):   
        dframe1 = []
        dframe2 = []
        rows_max_A = [ ]
        cols_max_A = [ ] 
        rows_max_B = [ ]
        cols_max_B = [ ]
        rows_max_C = [ ]
        cols_max_C = [ ]
        rows_max_D = [ ]
        cols_max_D = [ ]
        save_dir =  os.path.join(save_file_path, temp_files[files])             
        if not os.path.exists(save_dir):
               os.makedirs(save_dir)
        saved_data_files = os.path.join(file_path, temp_files[files],'Script_Data','saved_quads') 
        saved_dark_files = os.path.join(file_path_dark, temp_files[files],'Script_Data','saved_quads','Dark')
        
        all_int_files = [each for each in os.listdir(saved_data_files) \
                         if each.endswith('dat.sav')] 
        # all_dark_files = [each for each in os.listdir(saved_dark_files) \
                          # if each.endswith('dat.sav')] 
       
        for data_files in all_int_files:
            data_file = os.path.join(saved_data_files, data_files) 
            print(data_file)
            
            IDL_variable = readsav(data_file)
            data_path_name_split = data_files.split('_')
            int_time = int(data_path_name_split[-1].split('.')[0])        
            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            color = ['Blue','Green','Red','Orange']
            all_full_frame = IDL_variable.q
        
            ylim1= [0, 16000]
            ylim2 = [-6, 6]
            for i in range(0, 4):
                quad_full_frame = all_full_frame[:, i, :, :]
                avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
                active_quad = avg_quad[4:1028, 10:1034]
                tsoc = avg_quad[4:1028, 1034:1056]               
                
                #------perform bias subtraction using trailing overclocks and save the dark current image----------
                #bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc) 
                # mask out the outliers
               
                
                cross_talk_array = avg_quad
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
           
                    
                    if spot==1:
                        
                        dark_data_file= os.path.join(saved_dark_files, all_dark_files[0])
                   
                    elif spot==2:
                        dark_data_file= os.path.join(saved_dark_files, all_dark_files[1])
                        
                   
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
                    print(string1)
                    print(string2)
                    print(string3)
                    
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
                    smear_subtracted_quad[(smear_subtracted_quad>1500)] = np.mean(smear_subtracted_quad)
                    ylim = ylim2
                
                if spot == 1:
                    #rows, cols = np.reshape()
                    if i == 0: 
                        
                        rows_max_A.append(cross_talk_array[185:325, 170:310])
                        
                    elif i == 1:                        
                        rows_max_B.append(cross_talk_array[185:325, 170:310])
                        
                    elif i == 2:                        
                        rows_max_C.append(cross_talk_array[185:325, 170:310])
                       
                    elif i == 3:                        
                        rows_max_D.append(cross_talk_array[185:325, 170:310])
                          
                
                elif spot==2:
                    if i==0:
                        rows_max_A.append(cross_talk_array[720:860, 720:860])
                       
                    elif i==1:
                        rows_max_B.append(cross_talk_array[720:860, 720:860])
                       
                    elif i==2:
                        rows_max_C.append(cross_talk_array[720:860, 720:860])
                        
                    elif i==3:
                        rows_max_D.append(cross_talk_array[720:860, 720:860])
                      
                
                quad_save = 'Cross_Talk_Image_Ghost'
                save_dir_image = os.path.join(save_dir, quads[i], quad_save)
                if not os.path.exists(save_dir_image):
                    os.makedirs(save_dir_image)             
                figure_name = save_dir_image + '/'+ data_files + '.png'                  
                create_image(cross_talk_array, title1, figure_name, spot)
                
                #save_plot = 'plot_row_average'
                save_plot = 'plot_all_data'
                save_dir_plot = os.path.join(save_dir, quads[i], save_plot)
                if not os.path.exists(save_dir_plot):
                    os.makedirs(save_dir_plot)                 
                
                figure_name = save_dir_plot + '/'+ data_files + '.png' 
                xlabel = 'Pixel Indices (#)'
                plot_row_avg(cross_talk_array, title4, figure_name, color[i], xlabel)
                
                save_plot = 'plot_column_average'
                save_dir_plot = os.path.join(save_dir, quads[i], save_plot)
                if not os.path.exists(save_dir_plot):
                    os.makedirs(save_dir_plot)                 
                
                figure_name = save_dir_plot + '/'+ data_files + '.png' 
                xlabel = 'Spatial Pixel Indices (#)'
                #plot_row_avg(column_average, title3, figure_name, color[i], xlabel)
                
                save_plot = 'plot_row_average'
                save_dir_plot = os.path.join(save_dir, quads[i], save_plot)
                if not os.path.exists(save_dir_plot):
                    os.makedirs(save_dir_plot)                 
                
                figure_name = save_dir_plot + '/'+ data_files + '.png' 
                xlabel = 'Spectral Pixel Indices (#)'
                #plot_row_avg(row_average, title2, figure_name, color[i], xlabel)
        #cc        
#        dframe1 = pd.DataFrame(
#                  {'Quad_A_rows' : rows_max_A,
#                   'Quad_B_rows' : rows_max_B,
#                   'Quad_C_rows' : rows_max_C,
#                   'Quad_D_rows': rows_max_D,
#                   })
#        dframe2 = pd.DataFrame(
#                  {'Quad_A_cols' : cols_max_A,
#                   'Quad_B_cols' : cols_max_B,
#                   'Quad_C_cols' : cols_max_C,
#                   'Quad_D_cols': cols_max_D,
#                   })
        ndims, row_s,col_s = np.array(rows_max_A).shape 
        rows_max_A = np.reshape(np.array(rows_max_A), (ndims*row_s*col_s, 1))
        rows_max_B = np.reshape(np.array(rows_max_B), (ndims* row_s*col_s,1 ))
        rows_max_C = np.reshape(np.array(rows_max_C), (ndims*row_s*col_s, 1))
        rows_max_D = np.reshape(np.array(rows_max_D), (ndims*row_s*col_s, 1))
        
        csv_name_A = save_dir+'/'+temp_files[files]+'_cross_talk_A.csv'
        csv_name_B = save_dir+'/'+temp_files[files]+'_cross_talk_B.csv'
        csv_name_C = save_dir+'/'+temp_files[files]+'_cross_talk_C.csv'
        csv_name_D = save_dir+'/'+temp_files[files]+'_cross_talk_D.csv'
        
        np.savetxt(csv_name_A, np.asarray(rows_max_A), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_B, np.asarray(rows_max_B), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_C, np.asarray(rows_max_C), delimiter=',', fmt='%1.2f')
        np.savetxt(csv_name_D, np.asarray(rows_max_D), delimiter=',', fmt='%1.2f')
        
        #csv_name_cols = save_dir+'/'+temp_files[files]+'_cols_mean.csv'
        #dframe1.to_csv(csv_name_rows, header=True, columns=['Quad_A_rows','Quad_B_rows','Quad_C_rows','Quad_D_rows'])
        #dframe2.to_csv(csv_name_cols, header=True, columns=['Quad_A_cols','Quad_B_cols','Quad_C_cols','Quad_D_cols'])
       
        #cc
if __name__ == "__main__":
    main()
    