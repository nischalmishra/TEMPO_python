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



def fill_with_tsoc(active_quad_sto_filled, tsoc_image):
    print(active_quad_sto_filled.shape)
    spec_pix, spat_pix = active_quad_sto_filled.shape
              
    filled_quad = np.array([[0]*spec_pix]*spat_pix)
    filled_quad = np.reshape(filled_quad,(spec_pix, spat_pix))
    active_odd = active_quad_sto_filled[:, 1::2]
    active_even = active_quad_sto_filled[:, 0::2]
    
    tsoc_odd = tsoc_image[:, 1::2]
    tsoc_odd = tsoc_odd[:,0:10]
    tsoc_odd = np.mean(tsoc_odd)
    tsoc_even = np.mean(tsoc_image[:, ::2])
    
    filled_quad[:, ::2] = active_even+tsoc_even
    filled_quad[:, 1::2] = active_odd+ tsoc_odd   
    return filled_quad
    


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
    
    image_data[image_data>2000] =np.median(image_data)
#    image_data[image_data<0.018] = np.median(image_data)
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
    plt.plot(range1, even_samples_avg_image[::-1], 'r.', range2, even_samples_avg_sto[::-1],'b.')
    plt.plot(range1, odd_samples_avg_image[::-1], 'r.', range2, odd_samples_avg_sto[::-1],'b.')
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
    save_dir_local= r'C:\Users\nmishra\Workspace\TEMPO\Storage_region_analysis\Image_Sto_Comparisons_FPE_Temp\Image_storage_validation'

    temperature_files = [each for each in os.listdir(file_path) \
                        if each.endswith('_Darks')]

    for k in range(3, len(temperature_files)):


        image_dark_files = os.path.join(file_path, temperature_files[k],
                                        'Script_Data', 'saved_quads')

        sto_dark_files = os.path.join(file_path, temperature_files[k],
                                        'Dark_Imaging', 'saved_quads')

        dframe1 = pd.DataFrame()
        all_int_files_image = [each for each in os.listdir(image_dark_files) \
                             if each.endswith('.dat.sav')]
        all_int_files_sto = [each for each in os.listdir(sto_dark_files) \
                             if each.endswith('.dat.sav')]
        
        nominal_int_files_image = [items for items in all_int_files_image ]
        nominal_int_files_sto = [items for items in all_int_files_sto]
        
        save_dir = os.path.join(save_dir_local, temperature_files[k])
        if not os.path.exists(save_dir):
                os.makedirs(save_dir)
        
        for data_files_image, data_files_sto in zip(nominal_int_files_image[3:], nominal_int_files_sto[3:]):
                print(data_files_image)
                print(data_files_sto)
                
               
                data_path_name_split_image = data_files_image.split('_')
                int_time_image = round(int(data_path_name_split_image[-1].split('.')[0]))                    
                int_time_image = int(int_time_image)/1000
                  
               
                data_path_name_split_sto= data_files_sto.split('_')
                int_time_sto = data_path_name_split_sto[4][0:3]
                
                
                if 's' in int_time_sto:
                        int_time_sto = int(int_time_sto[0])*1000

                else:
                        int_time_sto = int(int_time_sto)
                
                data_file_image = os.path.join(image_dark_files, data_files_image)
                data_file_sto = os.path.join(sto_dark_files, data_files_sto)
                
                print(data_file_image)
                print(data_file_sto)
                
                IDL_variable_image = readsav(data_file_image)
                all_full_frame_image = IDL_variable_image.q
                print(all_full_frame_image.shape)
               
                IDL_variable_image= None
                
                IDL_variable_sto = readsav(data_file_sto)
                all_full_frame_sto = IDL_variable_sto.q
                print(all_full_frame_sto.shape)
                
               
                quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
                for i in range(0, 4): # 4 quads
                    quad_name = quads[i] 
                    print(int_time_image)
                    
                    quad_image_sto = np.mean(all_full_frame_sto[:, i, :, :], axis=0)
                                        
                    active_quad_sto = quad_image_sto[0:1028, 10:1034] # avg of 100 frames
                    tsoc_sto = quad_image_sto[0:1028, 1034:1056]           
                    bias_subtracted_quad_sto = perform_bias_subtraction(active_quad_sto, tsoc_sto)
                    
                    
                    print(bias_subtracted_quad_sto.shape)        
                    active_quad_sto_filled = np.sum(bias_subtracted_quad_sto[255:255+512, :], axis=0)
                    active_quad_sto_filled = np.reshape(active_quad_sto_filled, (1024, 1))
                    active_quad_sto_filled = fill_with_tsoc(np.array(active_quad_sto_filled), tsoc_sto)
                    #plt.plot(active_quad_sto_filled,'r')
                    #plt.show()
#                    cc
                    
                    
                    
                    
#                    
#                    create_image(bias_subtracted_quad_sto, title='a', figure_name='b')
#                    
#                    plt.plot(active_quad_sto[300, :], 'r--',label='300')
#                    plt.plot(active_quad_sto[350, :], 'b--',label='350')
#                    plt.plot(active_quad_sto[400, :], 'g--',label='400')
#                    plt.plot(active_quad_sto[500, :], 'm--',label='500')
#                    #plt.plot(active_quad_sto_filled,'k', label='Summed Storage Region (512 rows)')
#                    plt.title('Storage Region Profile from Dark Data, Quad B, -20C, int.time=2s' )
#                    plt.xlabel('Spatial Pixel Indices (#)')
#                    plt.ylabel('Raw DN')
#                    plt.legend(loc='best')
#                    plt.show()
#                    cc

                                   
                    quad_image_active = np.mean(all_full_frame_image[:, i, :, :], axis=0)                   
                    active_quad_image = quad_image_active[0:1028, 10:1034] # avg of 100 frames
                    tsoc_image = quad_image_active[0:1028, 1034:1056]
                    #active_quad_image =fill_with_tsoc(active_quad_image, tsoc_image) 
                    
                    
                    plt.plot(active_quad_sto_filled,'r', label='Estimated Storage summation row from SR ')
                    plt.plot(active_quad_image[0, :], 'k', label='Storage summation row from IR (First Line)')
                    plt.title('Storage Summation, quad B, -20C @int. time = 2S')
                    plt.legend(loc='best')
                    plt.show()
                    
                    #title = 'Serial overclock profile of Dark Data, Quad B, -20C, int.time=2s'
                   # plot_few_tsocs(tsoc_image, tsoc_sto,  title= title, figure_name='b')
                    
                    #print(active_quad_sto_filled.shape)
                    active_quad_sto_filled = fill_with_tsoc(bias_subtracted_quad_sto, tsoc_image)
                    #active_quad_sto_filled = np.ravel(active_quad_sto_filled.T)
#                    
                    #active_quad_sto_filled = fill_with_tsoc(active_quad_sto_filled, tsoc_sto)
                    
                    plt.plot(active_quad_sto_filled, label='Estimated Storage summation row from SR Image')
                    plt.plot(active_quad_image[0,:], label='Storage summation row from IR Image)')
                    #plt.plot(active_quad_image[1,:], label='Second Line (IR)')
                    #plt.plot(active_quad_sto_filled.T, label='Estimated Storage summation row')
                    plt.legend(loc='best')
                    plt.xlabel('Spatial Pixel Indices (#)')
                    plt.ylabel('DN')
                    plt.show()
                    cc
                    
                    plt.title('Comparisons of 2 Storage Summation Rows in Imaging Region\n  Quad B, temp=-20C, int. time = 118 ms')
                    plt.ylabel('Raw DN')
                    plt.xlabel('Spatial Pixel Indices (#)')
                    #plt.plot(np.mean(tsoc_image, axis=1), 'k')
                    #plt.plot(active_quad_image[1,:],'g')
                    plt.legend(loc='best')
                    plt.show()
                    cc
                    
                    
                    
#                    bias_subtracted_quad_active = perform_bias_subtraction(active_quad_image,
#                                                                      tsoc_image)
                    
                    #smear_subtracted_quad_active = perform_smear_subtraction(bias_subtracted_quad_active, int(int_time_image))
                   # smear_subtracted_quad_active  = active_quad_image
                   
                    
                    
                    
                    #smear_subtracted_quad_sto = perform_smear_subtraction(bias_subtracted_quad_sto, int(int_time_image))
                    #smear_subtracted_quad_sto  = active_quad_sto
                   
                    title = 'Image Region, QuadB, -18C, int.time = 118ms'
                    #create_image(active_quad_image, title, figure_name='b')
                    
                    
                    # to fill the storage region, skip 256 rows, sum 512 rows
                    #smear_subtracted_quad_sto = smear_subtracted_quad_sto[254:, :]
                  
                    
                    active_quad_sto_filled= np.sum(bias_subtracted_quad_sto[255:255+512, 10:1034], axis=0)/int(int_time_image)
                    #print(active_quad_sto_fill.mean())
#                    print(np.mean(tsoc_sto, axis=1).shape)
#                    cc
                    
                   
                    storage_line_active_quad = active_quad_sto_filled+np.mean(tsoc_sto[4:128, 0 ])+800
                    #plt.plot(active_quad_sto_filled[10:1034])
#                    plt.plot(storage_line_active_quad[10:1034],'r')
#                    plt.plot(active_quad_image[0, 10:1034],'k', label='1 (Image Region)')
#                    plt.show()
#                    cc
                    
                    
                   # print(np.mean(active_quad_sto[255:255+512, :]))
                   
                    
                    plt.plot(active_quad_image[0, :],'r', label='1 (Image Region)')
                    plt.plot(active_quad_image[1, :],'b', label='2 (Image Region)')
                    plt.plot(active_quad_sto_filled,'k', label='Summed Storage Region')
                    #plt.plot( active_quad_image[2, :],'g', label='3')
                    #plt.plot( active_quad_image[3, :],'y' ,label='4')
                    #plt.plot( active_quad_image[4, :], 'm',label='5')
                    plt.legend(loc='best')
                    #plt.title('Image Region Profile from Dark Data, Quad B, -20C, int.time=118ms' )
                    plt.show()
                    
                    
                    
                    plt.plot(active_quad_sto[300, :], 'r',label='300')
                    plt.plot(active_quad_sto[350, :], 'b',label='350')
                    plt.plot(active_quad_sto[400, :], 'g',label='400')
                    plt.plot(active_quad_sto[500, :], 'm',label='500')
                    plt.plot(active_quad_sto_filled,'k', label='Summed Storage Region (512 rows)')
                    plt.title('Storage Region Profile from Dark Data, Quad B, -20C, int.time=118ms' )
                    
                    plt.legend(loc='best')
                
                    #plt.plot(smear_subtracted_quad_sto, 'r')
                   
                    plt.show()
                    cc
                    
                    
                    


if __name__ == "__main__":
    main()