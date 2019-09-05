# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 07:24:40 2017
This script reads the fits file and analyzes the storage region behavior
@author: nmishra

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import matplotlib as mpl
#mpl.style.use('classic')
#mpl.rcParams['figure.facecolor'] = '0.75'
#mpl.rcParams['grid.color'] = 'k'
#mpl.rcParams['grid.linestyle'] = ':'
#mpl.rcParams['grid.linewidth'] = 0.5

from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

def make_quads_from_fits_file(data_file_name, data_file_full_path):
    """ Function to read the fits file
    """
    read_out_time = 0.095601715
    data_path_name_split = data_file_name.split('_')
    num_lines = data_path_name_split[4].split('-')[1]
    
    int_time = float(read_out_time)*int(num_lines)    
     
    print(int_time)
    cc
    temp = data_path_name_split[-2].split('-')[1]
    #print(temp)
    hdulist = fits.open(data_file_full_path)
    full_frame = hdulist[0].data
    # This data contains two quads, arrange in alternated fashion.
    ndims, rwos, cols = full_frame.shape
#    print(full_frame.shape)
#    cc
    quad_A = full_frame[0:ndims:2, :, :]
    quad_B = full_frame[1:ndims:2, :, :]
    quad_A = np.mean(quad_A, axis=0)
    quad_B = np.mean(quad_B, axis=0)
    full_frame = [quad_A, quad_B]
    return full_frame, int_time, temp


def split_frames_active_storage(full_frame):
    """ Function to return active and storage region from full_frame
    """
    quad_A = full_frame[0]
    quad_B = full_frame[1]
    storage_region_A = quad_A[0:1046, :]
    storage_region_B = quad_B[0:1046, :]
    image_region_A = quad_A[1046:, :]
    image_region_B = quad_B[1046:, :]
    # storage region has parallel 16 parallel overclocks and 2 storage region
    #to begin with. Let's
    #toss those out
    storage_region_A = storage_region_A[0:1028, :]
    storage_region_B = storage_region_B[0:1028, :]
    image_region = [image_region_A, image_region_B]
    storage_region = [storage_region_A, storage_region_B]
    return image_region, storage_region

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
    outlier_filtered_data = hist_data[measured_threshold < 10.]
    return outlier_filtered_data

def create_hist(active,storage, title, figure_name, xlabel):
    """ Create histogram of active and storage region after bias removal
    """
    
    active = filter_outlier_median(active)
    storage = filter_outlier_median(storage)
#    rows, cols = storage.shape
#    storage = np.reshape(storage, (rows*cols, 1))
    uncertainty_image = 100*np.std(active)/np.mean(active)
    uncertainty_storage = 100*np.std(storage)/np.mean(storage)
    
    label1 = 'Uncertainty (Image) = '+ str(round(uncertainty_image, 2))+'%'     
    label2 = 'Uncertainty (Storage) = '+ str(round(uncertainty_storage, 2))+'%'
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))
    fig.subplots_adjust(left= 0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)    
    ax[0].hist(active, bins=74, facecolor='blue', label=label1) 
    ax[0].grid(True, linestyle=':')
    legend = ax[0].legend(loc='best', ncol=1,edgecolor=None, shadow=True,
                        prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)    
    ax[0].set_ylabel('Frequency (# of pixels)', fontsize=12, fontweight="bold")
    ax[0].set_xlabel(xlabel, fontsize=11, fontweight="bold")
    ax[0].set_title(title, fontsize=12, fontweight="bold")
    ax[0].set_ylim(0, 90000)
    
    ax[1].hist(storage, bins=80, facecolor='red', label=label2) 
    ax[1].grid(True, linestyle=':')
    legend = ax[1].legend(loc='best', ncol=1,edgecolor=None, shadow=True,
                        prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)    
    ax[1].set_ylabel('Frequency (# of pixels)', fontsize=12, fontweight="bold")
    ax[1].set_xlabel(xlabel, fontsize=11, fontweight="bold")
    ax[1].set_ylim(0, 90000)
    #ax[1].set_title(title+ '(Storage Region)')
    plt.show()
    cc
    
    
    
    
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
#    plt.show()
#    cc
    plt.close('all')


def create_image(image_data, title, figure_name):
    """ Create image as a sanity check step
    """
    #plt.figure()
    fig, ax = plt.subplots(figsize=(10,10))
    image_data[image_data>0.09] = np.mean(image_data)
    image_data[image_data<0] = 0
    #image_data[image_data<0.04] = np.mean(image_data)
              
    image = ax.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cbar = fig.colorbar(image, extend='max')
    #cbar.cmap.set_over('green')
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
    cc
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')


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
    The main function
    """
    file_path1 = r'C:\Users\nmishra\Workspace\TEMPO\Storage_region_analysis\CCD_FullFrameDark'
    ccd_name = ['FT6BS_017_FF_Darks', 'FT6BS_014_FF_Darks']
 
    for ccd in ccd_name:
        
        file_path2 = os.path.join(file_path1, ccd)
        #print(file_path2)

        save_dir = file_path2 +'_analysis'
        all_int_files = [each for each in os.listdir(file_path2) \
                         if each.endswith('.fits')]
        dframe1 = pd.DataFrame()
        active_region_all = []
        storage_region_all = []
        quad_A_dark_current_image = []
        quad_A_dark_current_image_unct = []        
        quad_A_dark_current_storage = []
        quad_A_dark_current_storage_unct = []   
        quad_B_dark_current_image = []
        quad_B_dark_current_image_unct = []   
        quad_B_dark_current_storage = []
        quad_B_dark_current_storage_unct = []
        all_int = []
        all_temp = []
        for data_files in all_int_files[7:]:            
            data_file_path = os.path.join(file_path2, data_files)
            print(data_file_path)
            full_frame,  int_time, temp = make_quads_from_fits_file(data_files,
                                                                   data_file_path)
            image_region, storage_region = split_frames_active_storage(full_frame)
            string1 = 'Int. time = ' +str(round(int_time,2))+'ms, Temp = ' +str(temp)
            #print(string1)
            all_int.append(int_time)
            all_temp.append(temp[:-1])
          
            # Note that image region and storage region contains quad A and quad B
            if 'FT6BS_017' in file_path2:
                quads = ['Quad A', 'Quad B']
            elif 'FT6BS_014' in file_path2:
                quads = [' Quad C', 'Quad D']
            
            for i in range(0, 2):
                quad_image = np.array(image_region)[i, :, :]
                quad_image_active = quad_image[:, 10:1034]
                quad_image_tsoc = quad_image[:, 1034:]
              
                quad_storage = np.array(storage_region)[i, :, :]
                quad_storage_active = quad_storage[:, 10:1034]
                quad_storage_tsoc = quad_storage[:, 1034:]
                # let's compare the raw image
                hist_dir = 'compare_hist_active_region_before_offset'
                save_hist_plot = os.path.join(save_dir, hist_dir, quads[i])
                if not os.path.exists(save_hist_plot):
                       os.makedirs(save_hist_plot)
                figure_name =  save_hist_plot + '/'+ string1+ '_tsoc.png'
                title = 'Histogram of active regions\n ' + quads[i]+', ' +string1
                xlabel = 'Raw Signal (DN)'
#                create_hist(quad_image_active, quad_storage_active,
#                            title, figure_name, xlabel)
    
                # Ok, let's compare the trailing overclocks for the two region
                tsoc_comp = 'compare_overclocks'
                save_tsoc_plot = os.path.join(save_dir, tsoc_comp, quads[i])
                if not os.path.exists(save_tsoc_plot):
                       os.makedirs(save_tsoc_plot)
                figure_name =  save_tsoc_plot + '/'+ string1+ '_tsoc.png'
                title = 'Comparison of  trailing serial overclocks\n '+ quads[i]+', ' +string1
                #plot_few_tsocs(quad_image_tsoc, quad_storage_tsoc, title, figure_name)
    
                # Ok, let's perform the bias subtraction
                # Also , append the quads so that we can create images at the end
                bias_subtracted_quad_active = perform_bias_subtraction(quad_image_active,
                                                                      quad_image_tsoc)
                smear_subtracted_quad_active = perform_smear_subtraction(bias_subtracted_quad_active, int(int_time))
                #smear_subtracted_quad_active = smear_subtracted_quad_active[300:800, 300:800] 
                #smear_subtracted_quad_active = smear_subtracted_quad_active[:, ::2] 
                
                
                active_region_all.append(smear_subtracted_quad_active)
                bias_subtracted_quad_storage = perform_bias_subtraction(quad_storage_active,
                                                                      quad_storage_tsoc)
                smear_subtracted_quad_storage = perform_smear_subtraction(bias_subtracted_quad_storage, int(int_time))
                #smear_subtracted_quad_storage =  smear_subtracted_quad_storage[:, ::2]
                #smear_subtracted_quad_storage =  smear_subtracted_quad_storage[300:800, 300:800]
                                                                       
                #bias_subtracted_quad_storage = bias_subtracted_quad_storage[10:500, 100:300] 
                storage_region_all.append(smear_subtracted_quad_storage)
                # Ok, let's compare the histoggrams of active ans storage reigon
                hist_dir = 'compare_hist_active_region_after_smear'
                save_hist_plot = os.path.join(save_dir, hist_dir, quads[i])
                if not os.path.exists(save_hist_plot):
                       os.makedirs(save_hist_plot)
                figure_name =  save_hist_plot + '/'+ string1+ '_darkt_current_rate.png'
                title = 'Histogram of Dark Current Rate\n ' + quads[i]+', ' + string1
                xlabel = 'Dark Current Rate (DN/msecs)'
                create_hist(smear_subtracted_quad_active/float(int_time), smear_subtracted_quad_storage/int(int_time),
                           title, figure_name, xlabel)
                create_image(smear_subtracted_quad_storage/float(int_time), title,'Dark Current')
                
                if i == 0:
                    quad_A_dark_current_image.append(np.median(filter_outlier_median(smear_subtracted_quad_active/float(int_time))))                    
                    quad_A_dark_current_storage.append(np.median(filter_outlier_median(smear_subtracted_quad_storage/float(int_time))))
                    quad_A_dark_current_image_unct.append(np.std(filter_outlier_median(smear_subtracted_quad_active/float(int_time))))
                    quad_A_dark_current_storage_unct.append(np.std(filter_outlier_median(smear_subtracted_quad_storage/float(int_time))))
                elif i==1:
                    quad_B_dark_current_image.append(np.median(filter_outlier_median(smear_subtracted_quad_active/float(int_time))))
                    quad_B_dark_current_storage.append(np.median(filter_outlier_median(smear_subtracted_quad_storage/float(int_time))))
                    quad_B_dark_current_image_unct.append(np.std(filter_outlier_median(smear_subtracted_quad_active/float(int_time))))
                    quad_B_dark_current_storage_unct.append(np.std(filter_outlier_median(smear_subtracted_quad_storage/float(int_time))))
    
    
    
    #        quad_A_image = np.concatenate((np.array(storage_region_all)[0, :, :], np.array(active_region_all)[0,: , :]), axis=0)
    #        quad_B_image = np.concatenate((np.array(storage_region_all)[1, :, :], np.array(active_region_all)[1,: , :]), axis=0)
    #        full_image = np.concatenate((quad_A_image, quad_B_image), axis=1)
    #        #quad_A_image = np.array(storage_region_all)[0, :, :]
    #        print(np.min(full_image))
    #        print(np.max(full_image))
    #        print(np.mean(full_image))
    #        threshold = 150
    #        threshold1= -50
    #        data = np.ma.masked_greater(full_image, threshold)
    #        data = np.ma.masked_less(data, threshold1)
    #        title = 'Dark Current Image\n' + string1
    #        create_image(data, title,'b')
    #        rows, cols = full_image.shape
    #        plt.hist(np.reshape(full_image, (rows*cols,1)))
    #        plt.show()
    #        quad_B_image = np.concatenate((np.array(storage_region_all)[0, :, :], np.array(active_region_all)[0,: , :]), axis=0)
    #        storage_region =np.concatenate((np.array(storage_region_all)[0, :, :], np.array(storage_region_all)[1, :, :] ), axis=1)
    #        full_image = np.concatenate((storage_region, image_region), axis=0)
        dframe1 = pd.DataFrame(
                      {'Int. Time' : all_int,
                       'Temp(k)' : all_temp,
                       'Dark_current_image_A' : quad_A_dark_current_image,
                       'Dark_current_image_unct_A' : quad_A_dark_current_image_unct,                       
                       'Dark_current_storage_A' : quad_A_dark_current_storage,
                       'Dark_current_storage_unct_A' : quad_A_dark_current_storage_unct,                       
                       'Dark_current_image_B' : quad_B_dark_current_image,
                       'Dark_current_image_unct_B' : quad_B_dark_current_image_unct, 
                       'Dark_current_storage_B' : quad_B_dark_current_storage,
                       'Dark_current_storage_unct_B' : quad_B_dark_current_storage_unct,
                       })
    
        csv_file_name = save_dir+'/'+'dark_current_Even.csv'
        dframe1.to_csv(csv_file_name)


# end of main function

if __name__ == "__main__":
    main()
    

