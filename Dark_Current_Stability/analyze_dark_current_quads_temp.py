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
from mpl_toolkits.axes_grid1 import make_axes_locatable


def perform_bias_subtraction(active_quad, trailing_overclocks):

    """This function calculates the offset with trailing overclocks and subtracts
    the offset from the active area
    """
    spec_pix, spat_pix = active_quad.shape
    even_detector_bias = trailing_overclocks[ :, ::2]
    odd_detector_bias = trailing_overclocks[:, 1::2]
    bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
    # Toss out the outliers

    avg_bias_even = np.mean(even_detector_bias, axis=1)
    odd_detector_bias = trailing_overclocks[:, 1::2]
    # Toss out the outliers
    odd_detector_bias = odd_detector_bias[:, 0:10]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]
    odd_detector_active_quad = active_quad[:, 1::2]
    #Subtract off the offset from active region
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (spec_pix, spat_pix))
    # Reform the quad
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad

def perform_smear_subtraction(active_quad, int_time):
    """the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    """
    frame_transfer = 8
    smear_factor = (frame_transfer / (int_time+ frame_transfer))* np.mean(active_quad, axis=0)
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad

def filter_outlier_median(quads):
    """ Apart from the fixed mask, there are times when the outlier needs to be
    run in order to get good statistics. This will be used in conjunction to
    the fixed mask """
    if np.array(quads).ndim == 3:
        ndims, nx_quad, ny_quad = quads.shape
    elif np.array(quads).ndim == 2:
        ndims = 1
        nx_quad, ny_quad = quads.shape
    else:
        nx_quad = 1
        ndims = 1
        ny_quad = len(quads)
    hist_data = np.reshape(quads, (ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 6.]
    return outlier_filtered_data

def create_image(image_data, title, figure_name):
    """ Create image as a sanity check step
    """
    plt.figure()
    axes_ = plt.gca()
    image_data[image_data>0.03] = np.median(image_data)
    image_data[image_data<0.018] = np.median(image_data)
    image = axes_.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(axes_)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
    cc
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')

def create_hist(dark_current, title, figure_name, COLOR):
    """ Create histogram of adark current
    """

    dark_current = filter_outlier_median(dark_current)
   
#    rows, cols = storage.shape
#    storage = np.reshape(storage, (rows*cols, 1))
    uncertainty_dc = 100*np.std(dark_current)/np.mean(dark_current)
    

    label1 = 'Uncertainty (Image) = '+ str(round(uncertainty_dc, 2))+'%'
   
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
    fig.subplots_adjust(left= 0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    ax.hist(dark_current, bins=10,facecolor=COLOR, label=label1)
    ax.grid(True, linestyle=':')
    legend = ax.legend(loc='best', ncol=1,edgecolor=None, shadow=True,
                        prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    ax.set_ylabel('Frequency (# of pixels)', fontsize=12, fontweight="bold")
    ax.set_xlabel('Dark Current Rates (DN/ms)', fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.close('all')
    #ax.set_ylim(0, 90000)   
    plt.show()
    #cc
    

def plot_few_tsocs(tsoc_image, tsoc_storage,  title, figure_name):
    # let's take the mean tsoc for 100 frames
    even_samples_avg_image = np.mean(tsoc_image[:, ::2], axis=1)
    odd_samples_avg_image = tsoc_image[:, 1::2]
    #print(len(odd_samples_avg))
    odd_samples_avg_image = np.mean(odd_samples_avg_image[:, 0:10], axis=1)
    even_samples_avg_sto = np.mean(tsoc_storage[:, ::2], axis=1)
    odd_samples_avg_sto = tsoc_storage[:, 1::2]
    odd_samples_avg_sto = np.mean(odd_samples_avg_sto[:, 0:10], axis=1)
    nrows = 2
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 8))
    fig.subplots_adjust(left= 0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    #even_samples_avg = even_samples_avg
    ax[0].plot(even_samples_avg_image, 'b.', label='Image Region Even ')
    ax[0].plot(odd_samples_avg_image, 'r.', label='Image Region Odd')
    legend = ax[0].legend(loc='best', ncol=1, shadow=True,
                          prop={'size':12}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    ax[0].set_title(title+' (Image Region)', fontsize=12, fontweight='bold')
    ax[0].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')
    ax[0].grid(True, linestyle=':')
    ax[1].plot(even_samples_avg_sto, 'g.', label='Storage Region Even')
    ax[1].plot(odd_samples_avg_sto, 'm.', label='Storage Region Odd')
    ax[1].set_title(title+' (Storage Region)', fontsize=12, fontweight='bold')
    ax[1].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    ax[1].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')
    ax[1].grid(True, linestyle=':')
    legend = ax[1].legend(loc='best', ncol=1, shadow=True,
                          prop={'size':12}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.grid(True, linestyle=':')
    plt.show()
    plt.close('all')

    # let make the plot continuous
    range1 = np.arange(0, 1028)
    range2 = np.arange(1028, 2056)
    plt.plot(range1, even_samples_avg_image[::-1], 'b.', range2, even_samples_avg_sto[::-1],'g.')
    plt.plot(range1, odd_samples_avg_image[::-1], 'r.', range2, odd_samples_avg_sto[::-1],'m.')
    plt.grid(True, linestyle=':')
    plt.title(title, fontsize=12, fontweight='bold')
    plt.xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    plt.ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')


    plt.show()
    cc


def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2


    #outlier_mask = read_outlier_mask()
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    save_dir_local_image = r'C:\Users\nmishra\Workspace\TEMPO\Storage_region_analysis\Image_Sto_Comparisons\Image'
    save_dir_local_sto = r'C:\Users\nmishra\Workspace\TEMPO\Storage_region_analysis\Image_Sto_Comparisons\Storage'
     
    save_dir_local = [save_dir_local_image, save_dir_local_sto]

    temperature_files = [each for each in os.listdir(file_path) \
                        if each.endswith('_Darks')]

    for k in range(3, len(temperature_files)):


        image_dark_files = os.path.join(file_path, temperature_files[k],
                                        'Script_Data', 'saved_quads')

        sto_dark_files = os.path.join(file_path, temperature_files[k],
                                        'Dark_Imaging', 'saved_quads')

        all_dark_files = [image_dark_files, sto_dark_files]

        for num_files in range(1, len(all_dark_files)):
            all_int_time = []
            all_med_quad_A = []
            all_med_quad_B = []
            all_med_quad_C = []
            all_med_quad_D = []
            all_std_quad_A = []
            all_std_quad_B = []
            all_std_quad_C = []
            all_std_quad_D = []
            all_tsoc_image = []
            all_dark_current_image = []
            all_tsoc_storage = []
            all_dark_current_storage = []

            dframe1 = pd.DataFrame()
            all_int_files_image = [each for each in os.listdir(all_dark_files[num_files]) \
                                 if each.endswith('.dat.sav')]

            nominal_int_files = [items for items in all_int_files_image ]
            save_dir = os.path.join(save_dir_local[num_files], temperature_files[k])
            if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
            #print(save_dir)
            for data_files in nominal_int_files[2:]:
                   
                    data_path_name_split = data_files.split('_')
                    if num_files == 0:
                        int_time = round(int(data_path_name_split[-1].split('.')[0]))
                        int_time = int(int_time)/1000
                    elif num_files == 1:
                        int_time =  data_path_name_split[4][0:3]
                       # print(int_time[0])
                        if 's' in int_time:
                            int_time = int(int_time[0])*1000

                        else:
                            int_time = int(int_time)
                    #print('integ. time= ', int_time)
                    data_file = os.path.join(all_dark_files[num_files], data_files)
                    print(data_file)
                    print(int_time)
                    IDL_variable = readsav(data_file)
                    all_full_frame = IDL_variable.q
                    all_int_time.append(int_time)
                    quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
                    for i in range(0, 4): # 4 quads
                        quad_name = quads[i]          
                        quad = all_full_frame[:, i, :, :]
                        active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)
                        tsoc = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
                        bias_subtracted_quad = perform_bias_subtraction(active_quad, tsoc)
                        smear_subtracted_quad =  perform_smear_subtraction(bias_subtracted_quad, int(int_time))
                        #create_image(smear_subtracted_quad/int_time, title='a', figure_name='b')
                        if num_files == 0:
                            title = 'Image Region Dark Current Histogram\n' + str(quad_name) + ', Int. Time = '+ str(int_time) + ' ms, '+ 'Temp = ' + str(temperature_files[k][0:4])
                            COLOR ='blue'                    
                        elif num_files == 1:
                            title = 'Storage Region Dark Current Histogram\n' + str(quad_name) + ', Int. Time = '+ str(int_time) + ' ms, '+ 'Temp = ' + str(temperature_files[k][0:4])
                            COLOR='red'
                        figure_name = save_dir +'/' + quad_name.replace(" ","") + '_' + str(int(int_time)) + 'ms_hist_dark_current' 
                        
                        
                        
                        
                        # calculate the dark current rates
                        if int(int_time)==0:
                            smear_subtracted_quad = smear_subtracted_quad
                        else:
                            smear_subtracted_quad = smear_subtracted_quad/int(int_time)
                        
                       
                        create_hist(smear_subtracted_quad, title, figure_name, COLOR)
                        cc
                        unct = 10*np.std(filter_outlier_median(smear_subtracted_quad))/np.median(filter_outlier_median(smear_subtracted_quad))
                        if i == 0:
                            all_med_quad_A.append(np.mean(filter_outlier_median(smear_subtracted_quad)))
                            all_std_quad_A.append(unct)
                                           
                            
                        elif i == 1:
                            all_med_quad_B.append(np.mean(filter_outlier_median(smear_subtracted_quad)))
                            all_std_quad_B.append(unct)

                        elif i == 2:
                            all_med_quad_C.append(np.mean(filter_outlier_median(smear_subtracted_quad)))
                            all_std_quad_C.append(unct)

                        else:
                            all_med_quad_D.append(np.mean(filter_outlier_median(smear_subtracted_quad)))
                            all_std_quad_D.append(unct)

                        
#                        if num_files == 0: 
#                             all_tsoc_image.append(tsoc)
#                             all_dark_current_image.append(smear_subtracted_quad)
#                            
#                        elif num_files==1:
#                             all_tsoc_storage.append(tsoc)
#                             all_dark_current_storage.append(smear_subtracted_quad)
#                        
#                        active_quad = None
#                        bias_subtracted_quad = None
#                        smear_subtracted_quad = None
#
#            
#            plot_few_tsocs(all_tsoc_image, all_tsoc_storage, )
            
            dframe1 = pd.DataFrame(
                            {'Int_time.' : all_int_time,
                             'Avg_Quad_A' : all_med_quad_A,
                             'Avg_Quad_B' : all_med_quad_B,
                             'Avg_Quad_C' : all_med_quad_C,
                             'Avg_Quad_D' : all_med_quad_D,
                             'Var_Quad_A': all_std_quad_A,
                             'Var_Quad_B': all_std_quad_B,
                             'Var_Quad_C': all_std_quad_C,
                             'Var_Quad_D': all_std_quad_D,


                             })

            #dframe1.to_csv(save_dir+'/'+temperature_files[k]+'_Dark_current_rate.csv')


if __name__ == "__main__":
    main()