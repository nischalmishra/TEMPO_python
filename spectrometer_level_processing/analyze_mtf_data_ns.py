# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:12:07 2018

Code written to verify some of the TEMPO Bandpass requirements.


@author: nmishra
"""

import os
#from random import randint
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as tkr
#from scipy.interpolate import interp1d
#pylint: disable= E1101
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d
#from scipy.interpolate import spline
#from scipy import signal
#from scipy.optimize import curve_fit


#
#def perform_spline_interpol(x, y, x_smooth):
#    """Perform spline interpolation
#    """
#    y_smooth = spline(x, y, x_smooth)
#    return y_smooth


def gauss_function(wavelength, amplitude, center, width):
    """ Fit Gaussian function
    """
    return  amplitude * np.exp(-(wavelength-center)**2 / width)


def full_width_half_max(wavelen, data):
    """ Calculated FWHM
    """
    half_max = max(data) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - data[])
    inflection_point = np.sign(half_max - np.array(data[0:-1])) - \
                              np.sign(half_max - np.array(data[1:]))
    #plt.plot(wavelen[0:3839], inflection_point) #if you are interested
    #plt.show()
    #cc
    #find the left and right most indexes
    left_idx = np.where(inflection_point > 0)[0]
    #print(wavelen[left_idx])

    right_idx = np.where(inflection_point < 0)[-1]
    if wavelen[right_idx] - wavelen[left_idx] is None:
        fwhm_val = 0
    else:
        fwhm_val = wavelen[right_idx]-wavelen[left_idx]
        fwhm_val = fwhm_val[0]
    #print(FWHM)
    #print(wavelen[right_idx])
    #cc
    return fwhm_val, wavelen[left_idx][0], wavelen[right_idx][0]


def main():
    """ Main function begins from here
    """
    data_dir = r'C:\Users\nmishra\Desktop\390_nm_px\NS'
    file_name = 'mean_MTF_data_rev1.csv'
    data_file = os.path.join(data_dir, file_name)
    image_data = np.loadtxt(data_file, delimiter=",")
    
    
    
    image_data1 = image_data[0:len(image_data):9]
    image_data2 = image_data[1:len(image_data):9]
    image_data3 = image_data[2:len(image_data):9]
    image_data4 = image_data[3:len(image_data):9]
    image_data5 = image_data[4:len(image_data):9]
    image_data6 = image_data[5:len(image_data):9]
    image_data7 = image_data[6:len(image_data):9]
    image_data8 = image_data[7:len(image_data):9]
    image_data9 = image_data[8:len(image_data):9]
   
    
   
    
    #pixel_num = [2035, 2034, 2033, 2032, 2031, 2030, 2029, 2028, 2027]
    pixel_num = [19, 18, 17, 16, 15, 14, 13, 12, 11]

    #image_data = np.abs(image_data)
    #image_data = image_data[1036:1048, 1023]/np.max(image_data[1036:1048, 1023])

    # Ok, now let's smooth the data

#    pixel = np.arange(0, len(image_data))
#    resampled_pixel = np.arange(min(pixel), max(pixel), 0.1)
#    fs_ = 1/0.1
#    frequency = fs_ * np.arange(0, len(resampled_pixel))/len(resampled_pixel)*1/0.018
#    print(len(frequency))
#    fit_params = interp1d(pixel, image_data, kind='slinear')
#    image_data_smooth = fit_params(resampled_pixel)
#    fwhm, _, _ = full_width_half_max(resampled_pixel, image_data_smooth)
    #print
    #image_data_smooth = signal.detrend(image_data_smooth)

    #image_data_smooth= image_data
#    #print(fwhm)+

           
           
    plt.plot(image_data1/np.max(image_data1), 'b', linewidth=2, label='pixel # ' +str(pixel_num[0]))
    plt.plot(image_data2/np.max(image_data2), 'g', linewidth=2, label='pixel # ' +str(pixel_num[1]))    
    plt.plot(image_data3/np.max(image_data3), 'r', linewidth=2, label='pixel # ' +str(pixel_num[2]))
    plt.plot(image_data4/np.max(image_data4), 'orange', linewidth=2, label='pixel # ' +str(pixel_num[3]))
    plt.plot(image_data5/np.max(image_data5), 'm', linewidth=2, label='pixel # ' +str(pixel_num[4]))    
    plt.plot(image_data6/np.max(image_data6), 'purple', linewidth=2, label='pixel # ' +str(pixel_num[5]))    
    plt.plot(image_data7/np.max(image_data7), 'k', linewidth=2, label='pixel # ' +str(pixel_num[6]))    
    plt.plot(image_data8/np.max(image_data8), 'lime', linewidth=2, label='pixel # ' +str(pixel_num[7]))    
    plt.plot(image_data9/np.max(image_data9), 'brown', linewidth=2, label='pixel # ' +str(pixel_num[8]))    
    plt.legend(loc='best', fontsize=12)
    plt.title('pX : NS 390 nm, Single Pixel Spatial Response', fontsize=20)
    plt.ylabel('Single Pixel Normalized Count', fontsize=18)
    plt.xlabel('Image Number Counts (Temporal)', fontsize=18)
    plt.grid(True, linestyle=':')
    plt.xlim(0, len(image_data)/9+8)
    plt.show()


    plt.plot(resampled_pixel, image_data_smooth, 'b', label='Discrete Output')
    plt.title('COF : EW 540 nm, Spatial Response Function, Slit Width = ' + \
              str(round(fwhm, 2)) + ' pixels', fontsize=12)
    plt.ylabel('EW Spatial Response Function', fontsize=12)
    plt.xlabel('Number of spectral Pixels', fontsize=12)
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.show()

# Now for the MTF part.
    freq_n = 1/(2*0.018*fwhm) # Define Nyquist Frequency
    print('Nyquist_Freq=', freq_n)
    mtf_smooth = np.fft.fft(image_data_smooth)
    mtf_smooth = np.abs(mtf_smooth)

    #plt.closr(all)
    #mtf_smooth = mtf_smooth[0:len(mtf_smooth)/2]
    print(len(mtf_smooth))
    #cc
    #mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))
    dc_loc = np.where(mtf_smooth == max(mtf_smooth))[0]
    mtf_smooth = mtf_smooth[int(dc_loc):int(len(mtf_smooth))]
    frequency = frequency[int(dc_loc):int(len(frequency))]
    frequency = frequency - frequency[0]
    mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))
    resampled_freq = np.arange(np.min(frequency), np.max(frequency), 0.2)
    fit_params = interp1d(frequency, mtf_smooth, kind='slinear')
    nor_mtf_smooth = fit_params(resampled_freq)
    print(len(mtf_smooth))
    #mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))
    plt.plot(resampled_freq, nor_mtf_smooth, 'b.', label='EW MTF')
    plt.grid(True, linestyle=':')
    plt.axvline(x=freq_n, ymin=0, ymax=0.32, color='r', linestyle=':',
                linewidth=2, label='EW MTF Requirement')
    plt.plot(freq_n, 0.32, 'ro', markersize=8)
    plt.plot(freq_n, 0.51, 'bo', markersize=8, label='EW MTF @ Nyquist')
    plt.title('')
    plt.legend(loc='best', fontsize=14)
    plt.xlim(0, 18)
    plt.title('COF : EW 540 nm, EW MTF @ Nyquist = 0.51,' + \
              ' Allocation = 0.32', fontsize=12)
    plt.ylim(0, 1)
    plt.show()
    freq = np.fft.fftfreq(len(image_data_smooth))
    freq1 = freq[freq >= 0.0]
    loc = np.where(freq1 >= 0.0)
    mtf_smooth = mtf_smooth[loc]
    mtf_smooth = mtf_smooth/np.max(mtf_smooth)
    print(len(mtf_smooth))

    plt.plot(8*freq1/0.018, mtf_smooth, 'ro--', label='Discrete Output')
    plt.title('Normalized FFT of the Discrete Signal')
    plt.legend(loc='best')

    plt.grid(True, linestyle=':')
    plt.show()






if __name__ == "__main__":
    main()
