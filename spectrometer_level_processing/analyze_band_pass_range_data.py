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
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
#pylint: disable= E1101

def parse_wavelength_data(wavelen, wavelength_data):
    """ Parse the wavelenght file. The wavelength file was created manually
    from the spectrometer log file. It was a painful process
    """
    dframe = pd.read_csv(wavelength_data, sep=',')
    wavelengths = dframe[wavelen].dropna(how='all')
    return wavelengths.values


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
    return fwhm_val


def fit_gaussian_func(wavelength, radiance_data, wavelen_val):
    """ Function to fit the Gaussian fit to the wavelength Vs. Amplitude Data
    """
    wavelen_val = float(wavelen_val.strip('_')[0:3])
    radiance_data = radiance_data[:, 4]
    dframe = pd.DataFrame()
    dframe['wavelength'] = wavelength
    dframe['radiance_data'] = radiance_data
    dframe = dframe.groupby('wavelength', as_index=False).mean()
    dframe.reset_index(drop=True)
     # waveln example : 675_nm represents the header in the wavelength file
    wavelength = dframe['wavelength']
    radiance_data = dframe['radiance_data']
    radiance_data = radiance_data / np.max(radiance_data)
    width = 2
    center = wavelen_val
    amp = 1
    init_vals = [amp, center, width]
    best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data,
                                  p0=[init_vals])
    print(covarr)
    fitted_data = gauss_function(wavelength, *best_vals)
    fitted_data = fitted_data/np.max(fitted_data)
    #full_width_half_max(wavelength, radiance_data)
    plt.plot(wavelength, fitted_data, 'g', label=' Gaussian fit')
    plt.plot(wavelength, radiance_data, 'ro--', label='Wavelength averaged data')
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Image Counts (DN)')
    plt.title('Signal Vs Wavelength (736 nm, pixel # 500)')
    plt.legend(loc='best')
    #plt.ylim(0, 1)
    #plt.yscale('log')
    plt.show()


def fit_linear_interpolation(wavelength, radiance, wavelen_val):

    """ Function to fit the linear itnerpolation to the wavelength Vs. Amplitude Data
    """
    print(wavelen_val)
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength
    col_index = [20, 400, 1030, 1500, 2020]
    #plt.figure()
    color = ['red', 'blue', 'green', 'magenta', 'orange']
    count = 0
    performance = []
    for i in range(10, 2035):
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).min()
        dframe3.reset_index(drop=True)
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data/np.max(radiance_data)
        sampled_wavelength = np.arange(min(wavelength), max(wavelength), 0.0010)
        fit_params = interp1d(wavelength, radiance_data, kind='linear')
        fitted_val = fit_params(sampled_wavelength)
        x_val = np.where(fitted_val >= 0.5)
        initial_val = x_val[0][-1]
        initial_value = sampled_wavelength[initial_val]
        performance.append(initial_value)
        #print(initial_value)

    for i in col_index:
        #print(i)
        #i = 2044
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).min()
        dframe3.reset_index(drop=True)
        # waveln example : 675_nm represents the header in the wavelength file
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data/np.max(radiance_data)
        sampled_wavelength = np.arange(min(wavelength), max(wavelength), 0.0010)
        fit_params = interp1d(wavelength, radiance_data, kind='linear')
        fitted_val = fit_params(sampled_wavelength)
        plt.plot(wavelength, radiance_data, 'o-', color=color[count])
        #plt.plot(sampled_wavelength, fitted_val, color= COLOR[count])
        #plt.hold(True)
        count = count+1
    #plt.axhline(y=0.5, color='r', linestyle='--', label='TEMPO Requirement')
    mean_val = np.mean(np.array(performance))
    std_val = np.std(np.array(performance))
    per_val = 'Performance: '+ str(round(mean_val, 2))+ \
              ' nm'  +' +/- '+ str(round(std_val, 2)) + ' nm'
    plt.axvline(x=mean_val, color='lime', linestyle='-',
                label='TEMPO Performance', linewidth=2)
    plt.axvline(x=740, color='k', linestyle='-', label='TEMPO Requirement')
    plt.text(737.5, 0.63, 'Requirement: >740 nm', fontsize=10,
             fontweight='bold', style='italic',
             bbox={'facecolor':'yellow', 'alpha':0.6, 'pad':4})
    plt.text(737.5, 0.53, per_val, fontsize=10, fontweight='bold',
             style='italic', bbox={'facecolor':'yellow', 'alpha':0.6, 'pad':4})
    plt.legend(['Spatial Pixel 20', 'Spatial Pixel 400', 'Spatial Pixel 1030',
                'Spatial Pixel 1500', 'Spatial Pixel 2020',
                'Performance (50% pt.)', ' Requirement'], loc='best')
    plt.title('TEMPO Visiblr CCD, Short WL Spectral Range\n Half-Power Point = ' \
              + str(round(mean_val, 2))+ ' nm')
    plt.ylim(-0.2, 1.2)
    plt.xlim(736, 742)
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Normalized Image Counts (DN)', fontsize=12)
    plt.show()



def plot_bandpass_data(wavelength, radiance_data, wavelen_val):
    """ Plot few curves of spectral bandpass data
    """
    # waveln example : 675_nm represents the header in the wavelength file
    #print(wavelen_val)
    col_index = [10, 400, 1030, 1500, 2040]
    plt.figure()
    for i in col_index:
        plt.plot(wavelength, radiance_data[:, i]/max(radiance_data[:, i]),
                 '.', label='spatial index = '+ str(i))
    plt.title('TEMPO Spectral Bandpass Function, CW = '+ str(wavelen_val)+ ' nm\n'
              + '(Linear Scale)')
    plt.ylabel('Signal Counts (DN)')
    plt.xlabel('Wavelength (nm)')
    #plt.xticks(np.arange(np.min(round(wavelength, 1)), np.max(round(wavelength, 1)), 0.5))
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    #plt.hold(True)
    #plt.yscale('log')
    plt.show()
    #cc


def main():
    """ Main function begins from here
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectral_Range'
    data_wavelengths = sorted([each for each in os.listdir(data_dir)
                               if each.endswith('_nm_2')])
    for wavelen in data_wavelengths:
        # first let's read the spectral radiance data
        data_path = os.path.join(data_dir, wavelen, r'saved_quads\processed_image')
        radiance_file = os.path.join(data_path, r'bandpass_radiance.csv')
        #integrated_signal_file = os.path.join(data_path, r'integrated_radiance.csv')
        print(radiance_file)
        wavelen_val = wavelen[0:6]
        radiance_data = np.loadtxt(radiance_file, delimiter=",")
        # Now let's read the wavelength file. The wavelength file is a csv file
        #created from the image log file.
        wavelength_data = os.path.join(data_dir, r'wavelength_all.csv')
        wavelength = parse_wavelength_data(wavelen_val, wavelength_data)
        #integrated_signal = np.loadtxt(integrated_signal_file, delimiter = ",")
         ##the data may have some obvious outliers. Those need to be tossed out

         ##The wavelength data was acquired at finite spectral sampling. Hence
        ## a continuos curve is generated by fitting a Gaussian curve to the data


        fit_linear_interpolation(wavelength, radiance_data, wavelen)
       # plot_bandpass_data(wavelength, radiance_data, wavelen_val)
#        fit_gaussian_func(wavelength, radiance_data, wavelen)
        #plot_integrated_counts(integrated_signal,wavelength,  wavelen)
#    plot_FWHM_all_wavelength(data_dir)

if __name__ == "__main__":
    main()
