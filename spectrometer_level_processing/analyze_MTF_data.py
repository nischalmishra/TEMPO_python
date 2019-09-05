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
    data_dir = r'C:\Users\nmishra\Desktop\540_nm\EW'
    file_name = 'mean_MTF_data_rev1.csv'
    data_file = os.path.join(data_dir, file_name)
    image_data = np.loadtxt(data_file, delimiter=",")
    #image_data = np.abs(image_data)
#    image_data1 = image_data[0:len(image_data):9]
#    image_data2 = image_data[1:len(image_data):9]
    #image_data3 = image_data[2:len(image_data):9]
    #image_data6 = image_data[3:len(image_data):9]
    #image_data5 = image_data[4:len(image_data):9]
    #image_data6 = image_data[5:len(image_data):9]
    #image_data7 = image_data[6:len(image_data):9]
#    image_data8 = image_data[7:len(image_data):9]
#    image_data9 = image_data[8:len(image_data):9]

    image_data = np.abs(image_data)
    image_data = image_data[1036:1048, 1023]/np.max(image_data[1036:1048, 1023])
#    print(image_data[:])
    #cc
    #image_data = image_data[2:12]
    #image_data = image_data/np.max(image_data)
    pixel = np.arange(0, len(image_data))

    resampled_pixel = np.arange(min(pixel), max(pixel), 0.1)


    Fs = 1/0.1
#    frequency = np.arange(0, Fs/2+1, Fs/len(resampled_pixel))
#    #print(len(frequency))
    frequency = Fs * np.arange(0, len(resampled_pixel))/len(resampled_pixel)*1/0.018
    print(len(frequency))





    fit_params = interp1d(pixel, image_data, kind='slinear')
    image_data_smooth = fit_params(resampled_pixel)

    fwhm, _, _  = full_width_half_max(resampled_pixel, image_data_smooth)
    #print
    #image_data_smooth = signal.detrend(image_data_smooth)

    #image_data_smooth= image_data
#    #print(fwhm)+
#    #plt.plot(image_data1/np.max(image_data1), 'b', linewidth=2, label='pixel # 1020')
#    #plt.plot(image_data2/np.max(image_data2), 'g', linewidth=2, label='pixel # 1021')
#    #plt.plot(image_data3/np.max(image_data3), 'r', linewidth=2, label='pixel # 1022')
#    plt.plot(image_data4/np.max(image_data4), 'orange', linewidth=2, label='pixel # 1023')
#    plt.plot(image_data5/np.max(image_data5), 'm', linewidth=2, label='pixel # 1024')
#    plt.plot(image_data6/np.max(image_data6), 'purple', linewidth=2, label='pixel # 1025')
#    #plt.plot(image_data7/np.max(image_data7), 'k', linewidth=2, label='pixel # 1026')
#    #plt.plot(image_data8/np.max(image_data8), 'lime', linewidth=2, label='pixel # 1027')
#    #plt.plot(image_data9/np.max(image_data9), 'brown', linewidth=2, label='pixel # 1028')
#    #plt.legend(loc='best', fontsize=12)
#    plt.title('COF : NS 540 nm, Single Pixel Spatial Response', fontsize=20)
#    plt.ylabel('Single Pixel Normalized Count', fontsize=18)
#    plt.xlabel('Image Number Counts (Temporal)', fontsize=18)
#    plt.grid(True, linestyle=':')
#    plt.xlim(0, len(image_data)/9+8)
#    plt.show()


    #plt.plot(resampled_pixel-center, rect_pulse, 'ro-', label='Continuous Input')
    plt.plot(resampled_pixel, image_data_smooth, 'b', label='Discrete Output')
    #plt.axvline(x=center)
    plt.title('COF : EW 540 nm, Spatial Response Function, Slit Width = ' + \
              str(round(fwhm, 2)) + ' pixels', fontsize=12)
    plt.ylabel('EW Spatial Response Function', fontsize=12)
    plt.xlabel('Number of spectral Pixels', fontsize=12)
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.show()

    # creat a rectangular pulse


    #plt.show()

#
#
#    #image_data = fit_params(resampled_pixel)
#
##    print(len(image_data))
##    plt.plot(image_data)
##    plt.show()
##    input_data = [0]*9
##    input_data[4] = 1
##    input_data[5] = 1
##    input_data[3] = 1
#
#
##    plt.figure()
##    fig_ax = plt.gca()
##    mean_data[mean_data>=0.9*10*16383] = np.min(mean_data)
##    image = fig_ax.imshow(np.array(mean_data), cmap='nipy_spectral',
##                          origin='lower', interpolation='none')
##    plt.title ('Mean Image (61 images, MTF EW Test COF, Laser WL:390 nm)')
##    #plt.title(title)
##    divider = make_axes_locatable(fig_ax)
##    cax = divider.append_axes("right", size="5%", pad=0.05)
##    plt.colorbar(image, cax=cax)
##    plt.grid(False)
##    fig_ax.set_xlim(1015, 1035)
##    fig_ax.set_ylim(490, 500)
##    plt.show()
##
    #mean_data = mean_data[491:505, 1023] # EW
   # mtf_data = image_data
#    #MTF_data1 = mean_data[497, 980:1060]/np.max(mean_data[497, 1005:1040])


    #mean_data1 = image_data[495, 1010:1040]
    #mtf_data1 = mean_data1/np.max(mean_data1)
#    mean_data2 = image_data[497, 1010:1040]
#    mtf_data2 = mean_data2/np.max(mean_data2)
#    plt.plot(mtf_data, 'b')
#    #plt.plot(mtf_data1, 'r')
#    #plt.plot(mtf_data2, 'g')
##    #plt.plot(MTF_data1, 'r')
##    #plt.xlim(1005, 1040)
##    plt.grid(True, linestyle=':')
#    plt.show()
#
#
#    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
#    plt.show()
#    plt.close('all')
    #MTF_data = MTF_data/np.max(MTF_data)

    #cropped_image =  MTF_data/np.max(MTF_data[490:502, 1023])
    #cropped_image = MTF_data[:, 1023]
    #cropped_image = cropped_image[::-1]


    #plt.plot(cropped_image, 'ro-')

    #plt.xlim(490, 505)
    #plt.grid(True, linestyle=':')
    #plt.show()
    #n = int(len(cropped_image))
    #print(n)
    #k = np.arange(n)

    #Fs = 0.018000 # pixel pitch mm
    #T = n/Fs
    #frq = k/T
    #frq = frq[range(n//2)]
    #image_data = image_data[::-1]


    freq_N = 1/(2*0.018*fwhm)
    print('Nyquist_Freq=', freq_N)
#

    mtf_smooth = np.fft.fft(image_data_smooth)
    mtf_smooth = np.abs(mtf_smooth)



    #plt.closr(all)
    #mtf_smooth = mtf_smooth[0:len(mtf_smooth)/2]
    print(len(mtf_smooth))
    #cc
    #mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))
    DC_loc = np.where(mtf_smooth == max(mtf_smooth))[0]
    mtf_smooth = mtf_smooth[int(DC_loc):int(len(mtf_smooth))]
    frequency = frequency[int(DC_loc):int(len(frequency))]
    frequency = frequency - frequency[0]
    mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))

    resampled_freq = np.arange(np.min(frequency), np.max(frequency), 0.2)
    fit_params = interp1d(frequency, mtf_smooth, kind='slinear')
    nor_mtf_smooth = fit_params(resampled_freq)

    print(len(mtf_smooth))
    #mtf_smooth = np.abs(mtf_smooth)/np.max(np.abs(mtf_smooth))
    plt.plot(resampled_freq, nor_mtf_smooth, 'b.', label='EW MTF')
    plt.grid(True, linestyle=':')
    plt.axvline(x=freq_N, ymin=0, ymax=0.32, color='r', linestyle=':',
                linewidth=2, label='EW MTF Requirement')
    plt.plot(freq_N, 0.32, 'ro', markersize=8)
    plt.plot(freq_N, 0.51, 'bo', markersize=8, label='EW MTF @ Nyquist')
    plt.title('')
    plt.legend(loc='best', fontsize=14)
    plt.xlim(0, 18)
    plt.title('COF : EW 540 nm, EW MTF @ Nyquist = 0.51,' + \
              ' Allocation = 0.32', fontsize=12)
    plt.ylim(0, 1)
    plt.show()



    freq = np.fft.fftfreq(len(image_data_smooth))

    freq1 = freq[freq >= 0.0]
#   freq1 = freq1/max(freq)
#
#    #print(freq1)
#    #cc
    loc = np.where(freq1 >= 0.0)
    mtf_smooth = mtf_smooth[loc]

#
#        #MTF = MTF[range(n//2)]
#
#        #length = len(freq1)
#        #pixel =  np.arange(0, length)
   # resampled_pixel = np.arange(min(freq1), max(freq1), 0.01)
##    #print(resampled_pixel)
    #fit_params = interp1d(freq1, mtf_smooth, kind='linear')
    #mtf_smooth = fit_params(resampled_pixel)
    mtf_smooth = mtf_smooth/np.max(mtf_smooth)
    print(len(mtf_smooth))

    #resampled


    #resampled_freq = np.arange(np.min(freq2), np.max(freq2), 0.001)
    #fit_params = interp1d(freq2, nor_mtf, kind='slinear')
    #nor_mtf_smooth = fit_params(resampled_freq)
    #print(np.min(freq2), np.max(freq2))
    #cc
    #print(max(freq1)/(4*0.018))

    #mtf_val = mtf_smooth/mtf_smooth_i
#    resampled_pixel = 2*resampled_pixel/0.018

    plt.plot(8*freq1/0.018, mtf_smooth, 'ro--', label='Discrete Output')
   # plt.plot(freq2/0.018, mtf_smooth_i, 'bo--', label= 'Continuous Input')
    #plt.plot(resampled_freq,nor_mtf_smooth, 'mo-', label='Estimated MTF')
    plt.title('Normalized FFT of the Discrete Signal')
    plt.legend(loc='best')
    #plt.plot(freq1, np.true_divide(mtf_smooth,mtf_smooth_i), 'mo--')
    #plt.plot(mtf_val_input/np.max(mtf_val_input),'b*')
    plt.grid(True, linestyle=':')
    #plt.ylim(0,1)
    #plt.xlim(0, 0.5*max(freq1)/2*0.018)
    plt.show()
    #cc





    #plt.xlim(0, 60)

    plt.show()








if __name__ == "__main__":
    main()
