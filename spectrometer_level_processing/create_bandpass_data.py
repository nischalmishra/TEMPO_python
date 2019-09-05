# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:15:00 2019

@author: nmishra

This function creates TEMPO spectral bandpass functions. The idea to create
a look up table for bandpass and derive an equivalent continous bandpass function
to define the bandpass

"""

import os
import h5py
import numpy as np
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt

#pylint: disable= E1101
#pylint: disable-msg=R0912
#pylint: disable-msg=R0914
#pylint: disable-msg=R0915


def read_pixel_to_wavelen_map(data_dir):
    """ read pixel to wavelength map
    """
    pixel_to_wavelen = np.genfromtxt(os.path.join(data_dir, 'Pixel_To_Wavelen_Map.csv'),
                                     delimiter=',')
    return pixel_to_wavelen

def calculate_dark_current(data_dir):
    """ Function to calculate mean dark current"""
    dark_data_dir = os.path.join(data_dir, 'Dark_Data')
    processed_image_all = []
    dark_data_files = [each for each in os.listdir(dark_data_dir)
                       if each.endswith('.h5')]
    for dark_files in dark_data_files:
        data_file = h5py.File(os.path.join(dark_data_dir, dark_files), 'r')
        processed_image = data_file.get('Processed_data')
        processed_image_all.append(processed_image)
    dark_current = np.mean(np.array(processed_image_all), axis=0)
    return dark_current


def gauss_function(x_val, a_val, x0_val, sigma_val):
    """ Fit Gaussian function
    """

    return  a_val *np.exp(-(x_val-x0_val)**2/(2*sigma_val**2))

def main():
    """
    Main function begins from here. It call all the required functions
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Spectral_Band_Pass_Data'
    pixel_to_wavelen = read_pixel_to_wavelen_map(data_dir)
    all_wavelengths_file = [each for each in os.listdir(data_dir)
                            if each.endswith('_nm')]
    print(all_wavelengths_file)


    max_row_all = [2025, 1964, 1913, 1863, 1736, 1559, 1382, 1206, 1063,
                   1005,  686, 509, 331, 66, 23]
    for files in np.arange(8, len(all_wavelengths_file)):
        print(all_wavelengths_file[files])
        bandpass_val = []
        bandpass_val_norm = []
        bandpass_wvl_all = []

        bp_data_path = os.path.join(data_dir, all_wavelengths_file[files],
                                    r'processed_h5')
        bp_save_dir = os.path.join(data_dir, all_wavelengths_file[files],
                                   r'bandpass_lookup_table')
        if not os.path.exists(bp_save_dir):
            os.makedirs(bp_save_dir)
        dark_current = calculate_dark_current(bp_data_path)

        #max_row = 2025 #for 297.8 nm
        #max_row = 1964 #for 310 nm
        #max_row = 1913 #for 320 nm
        #max_row = 1863 #for 330 nm
        #max_row = 1736 #for 355 nm
        #max_row = 1559 #for 390 nm
        #max_row = 1382 #for 425 nm
        #max_row = 1206 #for 460 nm
        #max_row = 1063 #for 488.2 nm
        #max_row = 1005 #for 541.8 nm
        #max_row = 686  #for 605 nm
        #max_row = 509  #for 640 nm
        #max_row = 331  #for 675 nm
        #max_row = 66   #for 727.5 nm
        #max_row = 23   #for 736 nm

        bandpass_files = sorted([each for each in os.listdir(bp_data_path)
                                 if each.endswith('.h5')])
        for bandpass_files in bandpass_files:
           # print(bandpass_files)
            bandpass_wvl = []
            data_file = h5py.File(os.path.join(bp_data_path, bandpass_files), 'r')
            processed_image = data_file.get('Processed_data')
            processed_image = processed_image - dark_current
#            row_mean = np.mean(processed_image, axis=1)
#            max_row_each = np.where(row_mean == max(np.mean(processed_image, axis=1)))[0]
#            max_row_each = max_row_each[0]
#            max_row = max_row_each
#            print(max_row)
#            cc
            max_row = max_row_all[files]
            #print(max_row)
            subset_image = processed_image[max_row-18 : max_row+18, 5:2041]
            normalization_factor = np.sum(subset_image)
            subset_image_normalized = subset_image/normalization_factor
            subset_pixel_to_wvl = pixel_to_wavelen[max_row-18 : max_row+18, 5:2041]
            # Now lets fit the gaussian function to find the peak wavelength
            # for each data
            for i in range(0, subset_pixel_to_wvl.shape[1]):
                #print(i)
                #n = len(subset_image_normalized[:, i])
                #mean_wvl = np.sum(subset_pixel_to_wvl[:, i]*subset_image_normalized[:, i])/n
                mean_wvl = np.mean(subset_pixel_to_wvl[:, i])
                std_wvl = np.std(subset_pixel_to_wvl[:, i])
                best_val_interpol, covarr = curve_fit(gauss_function, subset_pixel_to_wvl[:, i],
                                                      subset_image_normalized[:, i]/np.max(subset_image_normalized[:, i]),
                                                      p0=[1, mean_wvl, std_wvl])
#                plt.plot(subset_pixel_to_wvl[:, i],
#                          subset_image_normalized[:, i]/np.max(subset_image_normalized[:, i]),
#                          'ro')
#                plt.plot(subset_pixel_to_wvl[:, i],
#                          gauss_function(subset_pixel_to_wvl[:, i],
#                                         *best_val_interpol), 'b')
#                plt.grid(True, linestyle=':')
#                plt.xlabel('Wavelength(nm)', fontsize=12)
#                plt.ylabel('Normalized Response', fontsize=12)
#                plt.title('Example of Gaussian fit to the Laser Signal Input',
#                          fontsize=12)
#                plt.show()
#                cc
                bandpass_wvl.append(best_val_interpol[1])
                #bandpass_wvl = np.array(bandpass_wvl)
            bandpass_wvl_all.append(bandpass_wvl)
            bandpass_val.append(processed_image[max_row, 5:2041])
            bandpass_val_norm.append(processed_image[max_row, 5:2041]/normalization_factor)
        save_bandpass = os.path.join(bp_save_dir, 'bandpass_radiance_V2.csv')
        save_bandpass_wvl = os.path.join(bp_save_dir, 'bandpass_wvl_V2.csv')
        save_bandpass_norm = os.path.join(bp_save_dir, r'bandpass_radiance_normalized_V2.csv')
        np.savetxt(save_bandpass, np.array(bandpass_val), delimiter=",")
        np.savetxt(save_bandpass_wvl, np.array(bandpass_wvl_all), delimiter=",")
        np.savetxt(save_bandpass_norm, np.array(bandpass_val_norm), delimiter=",")
        print(all_wavelengths_file[files] + 'Done')
       #cc
if __name__ == "__main__":
    main()
