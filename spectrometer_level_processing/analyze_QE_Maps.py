# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:49:31 2018

@author: nmishra
"""

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
from scipy.ndimage import median_filter



def filter_outlier_median(quads):
    """ Apart from the fixed mask, there are times when the outlier needs to be
    run in order to get good statistics. This will be used in conjunction to
    the fixed mask"""
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
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    return outlier_filtered_data

def make_quads_from_fits_file(file_path, data_path):
    data_path_name_split = data_path.split('_')
    int_time = data_path_name_split[-4]
    temp = data_path_name_split[-2]
    data_file = os.path.join(file_path, data_path)
    hdulist = fits.open(data_file)
    full_frame = hdulist[0].data
    return full_frame, int_time, temp



def make_full_frame_image(file_path, CCD_orientation):
    file_path_1 = file_path + '/' + CCD_orientation[0]
    file_path_2 = file_path+ '/' + CCD_orientation[1]
    all_four_quads = [ ]
    
    wavelength = '740nm'
    data_file_name_1 = 'FT6BS-017_01222016_130535_QE_' + wavelength + '_int_500_temp_253K_4.fits'
    data_file_name_2 = 'FT6BS-014_01212016_124519_QE_' + wavelength + '_int_500_temp_253K_4.fits'
    
    data_path_1 = os.path.join(file_path_1, data_file_name_1)
    full_frame, int_time, temp = make_quads_from_fits_file(file_path_1, data_path_1)
    total_coadds = full_frame.shape[0]
    quad_A = full_frame[0:total_coadds:2, :, :]
    avg_quad_A = np.mean(quad_A, axis=0)
    all_four_quads.append(avg_quad_A)    
    quad_B = full_frame[1:total_coadds:2, :, :]    
    avg_quad_B = np.mean(quad_B, axis=0)
    #create_image(avg_quad_B)
    all_four_quads.append(avg_quad_B)
    full_frame = avg_quad_A = avg_quad_B = None
    
    data_path_2 = os.path.join(file_path_2, data_file_name_2)
    full_frame, int_time, temp = make_quads_from_fits_file(file_path_2, data_path_2)
    total_coadds = full_frame.shape[0]
    quad_C = full_frame[0:total_coadds:2, :, :]    
    avg_quad_C = np.mean(quad_C, axis=0)
    all_four_quads.append(avg_quad_C)
    #create_image(avg_quad_C)
    quad_D = full_frame[1:total_coadds:2, :, :]    
    avg_quad_D = np.mean(quad_D, axis=0)
    all_four_quads.append(avg_quad_D)
    return all_four_quads, int_time, wavelength
    

def create_image(image_data, int_time, wave_len):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """
    plt.figure()
    fig_ax = plt.gca()
    #print(np.max(image_data))
    image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          origin='lower', interpolation='none')
      
    plt.title('Example of QE Image' + ', Int. time =' +str(int_time) +' ms, '  +' WL = ' + wave_len, fontsize=12)
    plt.xlabel('Spatial Pixel Indices')
    plt.ylabel('Spectral Pixel Indices')
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show()
    #plt.pause(0.1)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    #plt.close('all')


def perform_bias_subtraction(raw_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads = []
    num_quads, nx_quad, ny_quad = raw_quads.shape   
    for quads in range(0, num_quads):
        active_quad = raw_quads[quads, 2:1030, 10:1034]
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]
        spec_pix, spat_pix = active_quad.shape
        bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
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
        bias_subtracted_quad = np.reshape(bias_subtracted_quad, (spec_pix, spat_pix))
        bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
        bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
        all_quads.append(bias_subtracted_quad)
        
    return np.array(all_quads)

def calculate_dark_current(file_path):
    dark1 = np.loadtxt(file_path+'/' + 'Dark_Current1.csv', delimiter=",")
    dark2 = np.loadtxt(file_path+'/' + 'Dark_Current2.csv', delimiter=",")
    dark_current = (dark1+dark2) / 2
    return dark_current


def perform_smear_removal(active_quads, int_time):
    """Perform the SMEAR subtraction using Active Quad Method.
    The underlying assumption in smear subtraction is that the dark current
    in the storage region is really small and hence neglected from the analysis.
    typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor*tRO)
    tft = 8ms
    """
    all_quads = []
    #frame_transfer = 8
    tft = 8
    num_quads, nx_quad, ny_quad = active_quads.shape
#    print(active_quads.shape)
#    print(outlier_masks.shape)
#    cc
    for quads in range(0, num_quads):
        #print(quads)
        active_quad = active_quads[quads]        
        #print(ma.is_masked(active_quad))
        smear_factor = (tft / (int_time + tft))* active_quad.mean(axis=0)
        smear_subtracted_quad = active_quad - smear_factor[None, :]
        all_quads.append(smear_subtracted_quad)
        active_quad = None
    return all_quads



def remove_cross_talk(full_frame):
    """ Apply the cross talk correction. The correction factor was provided by
    is 0.0015 for quads within the same CCD and 0 outside the CCD
    """
    quad_a = full_frame[0, :, :] - 0.0015*full_frame[1, :, :]
    quad_b = full_frame[1, :, :] - 0.0015*full_frame[0, :, :]
    quad_c = full_frame[2, :, :] - 0.0015*full_frame[3, :, :]
    quad_d = full_frame[3, :, :] - 0.0015*full_frame[2, :, :]
    corrected_full_frame = [quad_a, quad_b, quad_c, quad_d]
    return corrected_full_frame

def read_prnu_files():
    """ The PRNU correction is based on each quad. These were derived using
    very bright and "homogenous" images acquired during the photon transfer
    integration time sweep
    """
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_map_Median'
    prnu_mask_a = np.genfromtxt(file_path +'/' + 'Quad A_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_b = np.genfromtxt(file_path +'/' + 'Quad B_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_c = np.genfromtxt(file_path +'/' + 'Quad C_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_d = np.genfromtxt(file_path +'/' + 'Quad D_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask = [prnu_mask_a, prnu_mask_b, prnu_mask_c, prnu_mask_d]
    return prnu_mask

def parse_prnu_file():
    """ Read the PRNU hdf file provided by BATC. It takes the spectrometer
        orientation including the overclocks.
        """
    hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
    file = h5py.File(hdf_name, 'r')
    prnu = file.get('prnu')
    prnu = np.array(prnu).transpose()
    quad_d = prnu[2:1030, 10:1034]
    quad_c = prnu[2:1030, 1078:2102]
    quad_a = prnu[1062:2090, 10:1034]
    quad_b = prnu[1062:2090, 1078:2102]
    prnu_map_lower = np.concatenate((quad_d, quad_c), axis=1)
    prnu_map_upper = np.concatenate((quad_a, quad_b), axis=1)
    prnu_map = np.concatenate((prnu_map_lower, prnu_map_upper), axis=0)
    return prnu_map



def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image
    """
    quad_a = full_frame[0, :, :]
    quad_b = full_frame[1, :, :]
    quad_c = full_frame[2, :, :]
    quad_d = full_frame[3, :, :]
    uv_ccd = np.concatenate((quad_d, np.fliplr(quad_c)), axis=1)
    visible_ccd = np.concatenate((np.flipud(quad_a), np.rot90(quad_b, 2)),
                                 axis=1)
    processed_image = np.concatenate((uv_ccd, visible_ccd), axis=0)
    return processed_image

def main():
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPS_QE'
    CCD_orientation = [ r'FT6BS-017_QE', r'FT6BS-014_QE']
    #dark_current = calculate_dark_current(file_path)
    #print(np.mean(filter_outlier_median(dark_current)))    
    #Note: FT6BS-017_QE : Visible  FT6BS-017_QE : UV
    all_four_quads, int_time, wavelength = make_full_frame_image(file_path, CCD_orientation)  
    #create_image(dark_current, int_time, wavelength)
    spectro_image = create_final_image(np.array(all_four_quads))
    #create_image(spectro_image, int_time, wavelength)
    bias_subtracted_quad = perform_bias_subtraction(np.array(all_four_quads))
    smear_removed_quad= perform_smear_removal(bias_subtracted_quad, int(int_time))
    cross_talk_removed_quad = remove_cross_talk(np.array(smear_removed_quad))
    
    
#    cross_talk_quad = np.array(cross_talk_removed_quad)[3, :, :]
#    cross_talk_quad [cross_talk_quad  >= 0.85*2**16] = 2**16
#    data = filter_outlier_median(cross_talk_quad)
#    plt.hist(data, bins=np.arange(min(data), max(data) + 20, 20), color='blue')
#    plt.grid(True, linestyle = ':')
#    plt.ylabel('Frequency')
#    plt.xlabel('Signal Counts (DN)')
#    print(100*np.std(data)/np.mean(data))    
#    plt.show()
#    cc
#    
    processed_image = create_final_image(np.array(cross_talk_removed_quad)) #- dark_current
    #processed_image = median_filter(processed_image, size=(5,5), mode='reflect')
      
    #processed_image[processed_image>25000] = np.mean(processed_image)
    create_image(processed_image, int_time, wavelength)
    index = [260, 800, 1500]
    plt.plot(processed_image[int(index[0]), :], 'b', label= 'Spectral Pixel Index = '+ str(index[0]))
    plt.plot(processed_image[int(index[1]), :], 'r', label= 'Spectral Pixel Index = '+ str(index[1]))
    plt.plot(processed_image[int(index[2]), :], 'lime', label= 'Spectral Pixel Index = '+ str(index[2]))
#    plt.plot(processed_image[int(index[3]), :], 'k', label= 'Spectral Pixel Index = '+ str(index[3]))
#    plt.plot(processed_image[int(index[4]), :], 'orange', label= 'Spectral Pixel Index = '+ str(index[4]))
   # plt.plot(processed_image[int(index[5]), :], 'magenta', label= 'Spectral Pixel Index = '+ str(index[5]))
    plt.title('Spatial Profile of a QE Image' + ', WL = ' + wavelength, fontsize=12)
    plt.grid(True, linestyle=':')
    plt.xlabel('Spatial Pixel Index (#)', fontsize=12)
    plt.ylabel('Signal Counts (DN)', fontsize=12)
    plt.legend(loc='best')
    plt.show()
    print(processed_image.shape)
    index = [100, 2000]
    spectral_range = np.arange(0, 2056)
    plt.plot(processed_image[:, int(index[0])],spectral_range, 'g', label= 'Spatial Pixel Index = '+ str(index[0]))
    plt.plot(processed_image[:, int(index[1])], spectral_range, 'm', label= 'Spatial Pixel Index = '+ str(index[1]))
#    plt.plot(processed_image[int(index[2]), :], 'lime', label= 'Spectral Pixel Index = '+ str(index[2]))
#    plt.plot(processed_image[int(index[3]), :], 'k', label= 'Spectral Pixel Index = '+ str(index[3]))
#    plt.plot(processed_image[int(index[4]), :], 'orange', label= 'Spectral Pixel Index = '+ str(index[4]))
   # plt.plot(processed_image[int(index[5]), :], 'magenta', label= 'Spectral Pixel Index = '+ str(index[5]))
    plt.title('Spectral Profile of a QE Image' + ', WL = ' + wavelength, fontsize=12)
    plt.grid(True, linestyle=':')
    plt.ylabel('Spectral Pixel Index (#)', fontsize=12)
    plt.xlabel('Signal Counts (DN)', fontsize=12)
    plt.legend(loc='best')
    plt.show()
    
    
    
    
    
if __name__ == "__main__":
    main()