# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""


import os
import numpy as np
import seaborn as sns
import pandas as pd
#from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from scipy.io.idl import readsav 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

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
    outlier_filtered_data = hist_data[measured_threshold < 3]
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
    
def plot_hist(data, plot_dir, figure_name, quad, int_time, color):
    save_dir = r'Hist_TSOC'+'/'+ quad
    plot_dir_hist = os.path.join(plot_dir, save_dir)
    if not os.path.exists(plot_dir_hist):
                os.makedirs(plot_dir_hist)
    if np.array(data).ndim > 1:
        ndims, nx_quad, ny_quad = data.shape
        data = np.reshape(data, (nx_quad*ny_quad*ndims, 1))
    
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
             '\n Integ. Time = '+ str(int_time)+ ' microsecs'                

    title = 'Histograms of TSOC (' +quad +')'
    #print(title)
    sns.set_context("talk")
    with sns.axes_style("darkgrid"):
        
        ax.hist(data, 23, normed=0, facecolor=color, alpha=1.75, label=label)
        ax.tick_params(axis='x', pad=10)
        ax.grid(True, linestyle=':')
        ax = plt.gca()
        legend = ax.legend(loc='best', ncol=3, shadow=True,
                           prop={'size':12}, numpoints=1)
        legend.get_frame().set_edgecolor('r')
        legend.get_frame().set_linewidth(2.0)
        ax.set_ylabel('Frequency (# of pixels)', fontsize=15,
                      fontweight="bold")
        ax.set_ylim(ymin=0)
        #ax.set_xlim(xmin=min(data))
        ax.set_xlabel('Counts (DNs)', fontsize=14, fontweight="bold")        
        ax.set_title(title, fontsize=14, fontweight="bold")
        plt.ticklabel_format(useOffset=False)
        fig.savefig(plot_dir_hist+'/'+'Hist_'+quad+'_'+ figure_name+'.png', dpi=100)
        plt.close('all')

def plot_dark_current_image(bias_subtracted_quad_A, plot_dir, data_files, quad, int_time, color):
    dark_current_image = np.mean(bias_subtracted_quad_A, axis=0)
    #dark_current_image = bias_subtracted_quad_A[0,:,:]
    #outliers= np.array(np.where([dark_current_image > 35]))
    #dark_current_image[outliers]=np.mean(dark_current_image)
    threshold=30
    threshold1 = -5
    plt.figure()
    image = plt.imshow(dark_current_image, cmap=plt.cm.hot, vmin = threshold1, vmax= threshold, origin='lower', interpolation='none')
    cbar = plt.colorbar(image, extend='max')
    cbar.cmap.set_over('green')
    title = 'Dark Current Image of '+quad + ' (Int.time = '+ str(int_time)+')' 
    plt.title(title, fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    #plt.show()
    save_dir = r'Dark_Current_image'+'/'+ quad
    plot_dir_dark_cur = os.path.join(plot_dir, save_dir)  
    if not os.path.exists(plot_dir_dark_cur):
                os.makedirs(plot_dir_dark_cur)
    plt.savefig(plot_dir_dark_cur+'/'+ data_files+'.png',dpi=100,bbox_inches="tight")
    plt.close('all')


def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    
    all_int_time = [ ]
    all_med_quad_A = [ ]
    all_med_quad_B = [ ]
    all_med_quad_C = [ ]
    all_med_quad_D = [ ]
    all_var_quad_A = [ ]
    all_var_quad_B = [ ]
    all_var_quad_C = [ ]
    all_var_quad_D = [ ]
    dframe1 = pd.DataFrame() 
    avg_quad_all = [ ]
    VA_setting = [ ]
    
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Intensity_Sweep\Light\saved_quads'
    nominal_int_files = [each for each in os.listdir(file_path) \
                         if each.endswith('.dat.sav')]
    save_dir = r'quad_images'
    plot_dir = os.path.join(file_path, save_dir)
    if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
    for data_files in nominal_int_files:           
            
            dframe1 = pd.DataFrame()
            data_path_name_split = data_files.split('_')  
            print(data_path_name_split)            
            int_time = round(int(data_path_name_split[-1].split('.')[0])) #for integ_sweep            
            VA =  data_path_name_split[0]
            #print(VA_setting)
            VA_setting.append(VA)
            data_file = os.path.join(file_path, data_files)

            IDL_variable = readsav(data_file)            
            all_full_frame = IDL_variable.q  
            all_int_time.append(int_time)
            
            quad_A = all_full_frame[:, 0, :, :]
            active_quad_A = quad_A[:, 4:1028, 10:1034]
            bias_A = quad_A[:, 4:1028, 1034:1056]           
            bias_subtracted_quad_A = perform_bias_subtraction (active_quad_A, bias_A)      
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_A)                    
            all_med_quad_A.append(np.mean(outlier_filtered_data))
            #print(bias_subtracted_quad_A.shape)
            diff_two_frames =  bias_subtracted_quad_A[1,:,:] - bias_subtracted_quad_A[0,:,:]
            sigma_s = np.std(diff_two_frames)/2
            all_var_quad_A.append(sigma_s)
            quad = 'Quad A'
            color='blue'
            plot_hist(bias_A, plot_dir, data_files,quad, int_time, color)
            plot_dark_current_image(bias_subtracted_quad_A, file_path, data_files, quad, int_time, color)
            
            
            
            outlier_filtered_data = None
            diff_two_frames = None
            sigma_s = None
 
            # For quad B
            quad_B = all_full_frame[:, 1, :, :]
            active_quad_B = quad_B[:, 4:1028, 10:1034]
            bias_B = quad_B[:, 4:1028, 1034:1056]
            bias_subtracted_quad_B = perform_bias_subtraction (active_quad_B, bias_B)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_B) 
            #outlier_filtered_data = filter_outlier_median(active_quad_B)                             
            #median_quad_B = np.mean(bias_subtracted_quad_B)           
            #all_med_quad_B.append(median_quad_B)
            all_med_quad_B.append(np.mean(outlier_filtered_data))
            diff_two_frames =  bias_subtracted_quad_B[1,:,:]- bias_subtracted_quad_B[0,:,:]
            #diff_two_frames = filter_outlier_median(diff_two_frames)
            sigma_s = np.std(diff_two_frames)/2
            all_var_quad_B.append(sigma_s)
            quad = 'Quad B'
            color='green'
            plot_hist(bias_B,plot_dir, data_files,quad, int_time, color)
            outlier_filtered_data = None
            diff_two_frames = None
            sigma_s = None
            
#            
#            # For quad C
            quad_C = all_full_frame[:,2,:,:]
            active_quad_C = quad_C[:, 4:1028, 10:1034]
            bias_C = quad_C[:, 4:1028, 1034:1056]
            bias_subtracted_quad_C = perform_bias_subtraction (active_quad_C, bias_C)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_C)
            #outlier_filtered_data = filter_outlier_median(active_quad_C)                                
            #median_quad_C = np.mean(bias_subtracted_quad_C)           
            #all_med_quad_C.append(median_quad_C)
            all_med_quad_C.append(np.mean(outlier_filtered_data))
            diff_two_frames =  bias_subtracted_quad_C[1,:,:] - bias_subtracted_quad_C[0,:,:]
            #diff_two_frames = filter_outlier_median(diff_two_frames)
            sigma_s = np.std(diff_two_frames)/2
            all_var_quad_C.append(sigma_s)
            quad = 'Quad C'
            color='red'
            plot_hist(bias_C,plot_dir, data_files,quad, int_time, color)
            outlier_filtered_data = None
            diff_two_frames = None
            sigma_s = None
            
            quad_D = all_full_frame[:,3,:,:]
            active_quad_D = quad_D[:, 4:1028, 10:1034]
            bias_D = quad_D[:, 4:1028, 1034:1056]
            bias_subtracted_quad_D = perform_bias_subtraction (active_quad_D, bias_D)
            outlier_filtered_data = filter_outlier_median(bias_subtracted_quad_D) 
            #outlier_filtered_data = filter_outlier_median(active_quad_D)                             
            #median_quad_D = np.mean(bias_subtracted_quad_D)           
            #all_med_quad_D.append(median_quad_D)
            all_med_quad_D.append(np.mean(outlier_filtered_data))
            diff_two_frames =  bias_subtracted_quad_D[1,:,:] - bias_subtracted_quad_D[0,:,:]
            #diff_two_frames = filter_outlier_median(diff_two_frames)
            sigma_s = np.std(diff_two_frames)/2
            all_var_quad_D.append(sigma_s)
            quad = 'Quad D'
            color = 'magenta'
            plot_hist(bias_D,plot_dir, data_files,quad, int_time, color)
            outlier_filtered_data = None
            diff_two_frames = None
            sigma_s = None
            data = [np.mean(active_quad_A), np.mean(active_quad_B), np.mean(active_quad_C), np.mean(active_quad_D)]
            avg_quad_all.append(np.mean(data))
            
            
             #let's save quad images too
             
#            quad_a = np.mean(bias_A, axis=0)            
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
       
#    print('A:', all_med_quad_A)
#    print('B:', all_med_quad_B)
#    print('C:',  all_med_quad_C)
#    print('D:', all_med_quad_D)
#         
    dframe1 = pd.DataFrame(
                    {'Int_time.' : ['%.2f' % elem for elem in all_int_time],
                     'VA_Setting': VA_setting,
                     'Avg_Active_Area': ['%.2f' % elem for elem in avg_quad_all],
                     'Avg. Quad A' : ['%.2f' % elem for elem in all_med_quad_A],
                     'Avg. Quad B' : ['%.2f' % elem for elem in all_med_quad_B],
                     'Avg. Quad C' : ['%.2f' % elem for elem in all_med_quad_C],
                     'Avg. Quad D' : ['%.2f' % elem for elem in all_med_quad_D],
                     'Avg. Var_Quad_A': ['%.2f' % elem for elem in all_var_quad_A],
                     'Avg. Var_Quad_B': ['%.2f' % elem for elem in all_var_quad_B],
                     'Avg. Var_Quad_C': ['%.2f' % elem for elem in all_var_quad_C],
                     'Avg. Var_Quad_D': ['%.2f' % elem for elem in all_var_quad_D]                     
                     })
                
    dframe1.to_csv(plot_dir+'/'+'DN_Vs_Variance_Intensity_sweep'+'.csv')
if __name__ == "__main__":
    main()