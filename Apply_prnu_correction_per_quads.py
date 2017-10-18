# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""
# This script uses the PRNU map derived from integration time sweep data to some
# of the images acquired during intesnity sweep.

from random import randint
import os
import math
from matplotlib.pyplot import cm 
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from create_PRNU_map_each_quads import  perform_bias_subtraction_ave, \
                                         create_striping_metric, \
                                         create_striping_metric_plot,\
                                         perform_smear_subtraction

def create_image_same_scale(image1, corrected_image, title, figure_name):
    #plt.figure()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    plt.subplots_adjust(wspace=0.05, hspace=0.9)
    fig.suptitle(title, fontsize=16, y=.75)
    print('mean_original', np.mean(image1))
    print('mean_corrected', np.mean(corrected_image))
    im = axes[0].imshow(image1, origin='lower')
    axes[0].set_title('Original Image')
    axes[0].autoscale(False)
    clim = im.properties()['clim']
    axes[1].imshow(corrected_image, origin='lower', clim=clim)
    axes[1].set_title('After Applying PRNU')
    axes[0].set_axis_off()
    axes[1].set_axis_off()
    plt.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.4, pad=0.01)
    plt.savefig(figure_name)

def calculate_fft(active_quad_A, corrected_image):
    """Compares FFT before and after PRNU correction
    """
    single_line = active_quad_A[300,:]
    fft_val = np.fft.fft(single_line)    
    fft_freq = np.fft.fftfreq(len(single_line))
    i = fft_freq>0
    
    single_line_corrected = corrected_image[300,:]
    fft_val_corrected = np.fft.fft(single_line_corrected)    
    fft_freq_corrected = np.fft.fftfreq(len(single_line_corrected))
    j = fft_freq_corrected>0
    fig, ax = plt.subplots(nrows=2,ncols=1,sharey=True)
    dc_val = fft_val[i]
    dc_comp = dc_val[0]
    print(np.abs(fft_val[i])[0])
    dc_val1 = fft_val_corrected[j]
    dc_comp1 = dc_val1[0]
    print(np.abs(dc_comp1))
    #cc
    ax[0].plot(fft_freq[i], np.abs(fft_val[i]/np.abs(dc_comp)),'r')
    ax[0].legend(['Original Image'])
    ax[0].set_title('FFT along Spectral Direction, Int time =125047 microsecs, Quad A')
    #ax[0].set_ylabel('Amplitude (DC Normalized)')
    ax[1].plot(fft_freq_corrected[j], np.abs(fft_val_corrected[j])/np.abs(dc_comp1),'g', label='PRNU Corrected Image')
    ax[1].legend(['After applying PRNU'])
    #ax[1].set_ylim(0,1500)
    fig.text(0.05, 0.5, 'Amplitude (DC Normalized)', va='center', rotation='vertical', fontsize=14)
    #ax[1].set_ylabel('Amplitude (DC Normalized)',fontsize=14,fontweight="bold")
    ax[1].set_xlabel('Spatial Frequency (HZ)')
    ax[0].set_xlim(0, 0.5)
    ax[1].set_xlim(0, 0.5)
    ax[0].grid(True, linestyle=':')
    ax[1].grid(True, linestyle=':')
    plt.show()
    
 
def calculate_fft_each_quad(active_quad):
    column_average = np.mean(active_quad, axis=1)
    fft_val = np.fft.fft(column_average)
    fft_freq = np.fft.fftfreq(len(column_average))
    i = fft_freq>0
    return fft_freq[i],np.abs(fft_val[i])/np.abs(fft_val[i])[0]       


def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    if 'Integration_Sweep' in file_path1:
        saturated_collects = ['FT6_SHORT_INT_0.dat.sav','FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']

    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
                     

    nominal_int_files = [items for items in all_int_files if items not in saturated_collects]
    orig_directory = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads_validation/Filtered_data'
    
    PRNU_directory = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads'
    if not os.path.exists(save_dir):
                os.makedirs(save_dir)
    
    for i in range(0, 3): # for the 4 quads
        i=2
        all_freqs = [ ]
        all_amplitudes = [ ]
       
        for data_files in nominal_int_files:
            print(data_files)
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)

            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                string1 = 'VA_'
                string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '

            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data

            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q

            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
            
            
            active_quad_A = perform_bias_subtraction_ave(quad_A[4:1028, 10:1034], tsoc_A[:, 2:])
            active_quad_A, smear_signal = perform_smear_subtraction(active_quad_A[9:1000, :], int_time)
            quad_dir = quads[i]

            PRNU_maps = quads[i] +'/'+'Final_PRNU\Filtered_PRNU'
            PRNU_file = [each for each in os.listdir(os.path.join(PRNU_directory, PRNU_maps)) \
                        if each.endswith('_Final_PRNU.csv')]
            PRNU_mask = os.path.join(PRNU_directory, PRNU_maps, PRNU_file[0])
           # PRNU = genfromtxt(PRNU_mask, delimiter=',')
            #print(PRNU.shape)
           # PRNU_corrected_image = np.true_divide(active_quad_A, PRNU)
            #calculate_fft(active_quad_A, PRNU_corrected_image)
             # spatial FFT
            freq_spatial, amplitude_spatial = calculate_fft_each_quad(active_quad_A)
           
            
            #plt.figure()
            plt.plot(freq_spatial,amplitude_spatial)
            plt.grid(True, linestyle=':')
            plt.title('Spectral Dimension Fourier Mode: '+ quads[i] )
            plt.xlabel('Frequency (HZ)')
            plt.ylabel('DC Normalized Amplitude')
            plt.xlim(0,0.6)
            plt.ylim(0,0.6)
            
            
            all_freqs.append(freq_spatial)
            all_amplitudes.append(amplitude_spatial)
            

           #*******************************************************************************************
            
            orig_image_dir = 'saved_quads'
            plot_dir = os.path.join(save_dir, quad_dir, orig_image_dir)
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
                

           #Let's plot the image of the active quad
            figure_name = plot_dir+'/'+ string1+ str(int_time)+'.png'
            image_size = ', size = '+ str(active_quad_A.shape)
            title = 'Active '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'+ image_size
            #create_image_same_scale(active_quad_A, PRNU_corrected_image, title, figure_name)
            
            #****************STRIPING METRIC SECTION**************************

            # let us plot the striping metric for each integration time
            # begin with histograms
            striping_directory_corrected = 'striping_metric_hist'
            striping_met_directory = os.path.join(save_dir, quad_dir, striping_directory_corrected)
            striping_met_directory_orig = os.path.join(orig_directory, quad_dir, striping_directory_corrected)
            if not os.path.exists(striping_met_directory):
                os.makedirs(striping_met_directory)

            if not os.path.exists(striping_met_directory_orig):
                os.makedirs(striping_met_directory_orig)
            figure_name = striping_met_directory +'/'+ string1+ str(int_time)+'.png'
            figure_name_orig = striping_met_directory_orig +'/'+ string1 +str(int_time)+'.png'

            title = 'Striping Metric In TEMPO FPS Level Data After PRNU, '+ quads[i]+', ' +string2 + str(int_time)+ r" $\mu$" +'secs'
            title_orig = 'Striping Metric In TEMPO FPS Level Data Before PRNU, '+ quads[i]+', ' +string2 + str(int_time)+ r" $\mu$" +'secs'
            nx_quad, ny_quad = active_quad_A.shape
            sample_rows = randint(0, nx_quad)
            sample_cols = randint(0, ny_quad)
           # striping_metric_all_rows, striping_metric_all_cols = create_striping_metric(sample_rows, sample_cols, PRNU_corrected_image, title, figure_name)
            #create_striping_metric(sample_rows, sample_cols, active_quad_A, title_orig, figure_name_orig)
#
            # now the actual plot
            striping_plot_directory = 'striping_metric_plots'
            striping_met_plot_directory = os.path.join(save_dir, quad_dir, striping_plot_directory)
            if not os.path.exists(striping_met_plot_directory):
                os.makedirs(striping_met_plot_directory)
            figure_name = striping_met_plot_directory +'/'+ string1+ str(int_time)+'.png'
            title = 'Striping metric plot after PRNU, '+quads[i] +', '+ string2+ str(int_time) + r" $\mu$" +'secs'
            #create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name)
        plt.show()  
        cc
       

if __name__ == "__main__":
    main()
    