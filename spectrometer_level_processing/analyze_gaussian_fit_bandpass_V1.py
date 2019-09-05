# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:15:00 2019

@author: nmishra

This function creates TEMPO spectral bandpass functions. The idea to create
a look up table for bandpass and derive an equivalent continous bandpass function
to define the bandpass. The two functions being expolored on this version is the
Gaussian Flat top and modified Cauchy function. This is work under progress.

"""

import os
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import pandas as pd


#pylint: disable= E1101
#pylint: disable-msg=R0912
#pylint: disable-msg=R0914
#pylint: disable-msg=R0915




def flat_top_gaussian(wavelength, amp, center, width):
    """
    Fit the flat top Gaussian FUnction   """

    gauss_fxn = amp * np.exp(-(wavelength-center)**4/ 2*width**4)
    return gauss_fxn



def gaussian(wavelength, amp, center, width):
    """
    Fit the  Gaussian FUnction
    """

    gauss_fxn = amp * np.exp(-(wavelength-center)**2/ 2*width**2)
    return gauss_fxn


def cauchy_function(wavelength, center, sigma):
    """ FIt the cauchy function
    """
    return 1/(pi)*(7*sigma**2/(2*(wavelength-center)**4+ 2*sigma**2))


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
        FWHM = 0
    else:
        FWHM = wavelen[right_idx]-wavelen[left_idx]
        FWHM = FWHM[0]

    return FWHM


def smooth_data_interpol(wavelen, bp_data):
    """ smooth the data with linear interpolation if you are interested
    """
    sampled_wavelen = np.arange(min(wavelen), max(wavelen), 0.05)
    fit_params = interp1d(wavelen, bp_data, kind='linear')
    fitted_val = fit_params(sampled_wavelen)
    return sampled_wavelen, fitted_val


def fit_cauchy_function(bandpass_data, bandpass_wvl, save_dir, wavelen_val):
    """ Analyze Cauchy Function
    """
    cauchy_fit_params_all1 = []
     #cauchy_fit_params_all2 = []
     #dframe_cauchy = pd.DataFrame()
    for i in range(0, np.shape(bandpass_data)[1]):
        normalized_data = bandpass_data[:, i]/np.max(bandpass_data[:, i])
        wavelength = bandpass_wvl[:, i]
        width = np.std(wavelength)
        center = np.mean(wavelength)
        init_vals = [center, width]
        best_vals, covarr = curve_fit(cauchy_function, wavelength,
                                      normalized_data, p0=[init_vals])
        fitted_data_cauchy = cauchy_function(wavelength, *best_vals)

        cauchy_fit_params_all1.append(best_vals)
         #Cw_all2.append(best_vals1[1])
        fitted_data_cauchy = fitted_data_cauchy/np.max(fitted_data_cauchy)

        plt.figure(1)
        plt.plot(wavelength, fitted_data_cauchy, 'b.:', label=' Cauchy fit')
        plt.plot((wavelength), normalized_data, 'ro', label='Bandpass Data')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Image Counts (DN)')
        plt.title('Signal Vs Wavelength (Along 675 nm, pixel # 1500)')
        plt.legend(loc='best')

        plt.figure(2)
        plt.plot(wavelength, fitted_data_cauchy-normalized_data, 'k.',
                 label='fit- data')
        plt.ylim(-0.06, 0.06)
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Image Counts (DN)')
        plt.title('Difference Between Fit and Actual Data (Along 675 nm, pixel # 1500)')
        plt.legend(loc='best')
        plt.show()
        #cc

#     gauss_fit_params_all1 = np.array(gauss_fit_params_all1)
#     gauss_fit_params_all2 = np.array(gauss_fit_params_all2)
#     dframe_gauss = pd.DataFrame({'A1' : gauss_fit_params_all1[:, 0],
#                                 'CW1' : gauss_fit_params_all1[:, 1],
#                                 'Sigma1' : gauss_fit_params_all1[:, 2],
#                                  'A2' : gauss_fit_params_all2[:, 0],
#                                 'CW2' : gauss_fit_params_all2[:, 1],
#                                 'Sigma2' : gauss_fit_params_all2[:, 2]
#                                    })
#
#     dframe_gauss.to_csv(save_dir +'/' + 'Gaussian_fit_params_' + str(wavelen_val) + '_nm.csv'))

def  fit_double_gaussian_func(bandpass_data, bandpass_wvl, save_dir, wavelen_val):
    """ Analyze Gaussian Fits
    """
    gauss_fit_params_all1 = []
    gauss_fit_params_all2 = []
    dframe_gauss = pd.DataFrame()
    for i in range(0, np.shape(bandpass_data)[1]):
        normalized_data = bandpass_data[:, i]/np.max(bandpass_data[:, i])
        wavelength = bandpass_wvl[:, i]
        width = np.std(wavelength)
        center = np.mean(wavelength)
        amp = 1
        init_vals = [amp, center, width]
        best_vals, covarr = curve_fit(flat_top_gaussian, wavelength,
                                      normalized_data, p0=[init_vals])
        fitted_data_flat_top = flat_top_gaussian(wavelength, *best_vals)

        best_vals1, covarr1 = curve_fit(gaussian, wavelength,
                                        normalized_data, p0=[init_vals])

        fitted_data_gauss = gaussian(wavelength, *best_vals1)
         #print(best_vals, best_vals1)
        gauss_fit_params_all1.append(best_vals)
        gauss_fit_params_all2.append(best_vals1)
         #Cw_all2.append(best_vals1[1])
        fitted_data_flat_top = fitted_data_flat_top/np.max(fitted_data_flat_top)
        fitted_data_gauss = fitted_data_gauss/np.max(fitted_data_gauss)
       # final_fit = fitted_data_flat_top + fitted_data_gauss
        final_fit = fitted_data_flat_top
        final_fit = final_fit/np.max(final_fit)
        #%--------------------------------------------------------------%%%%%
        plt.figure(1)
        plt.plot(wavelength, final_fit, 'b.:', label=' Gaussian fit')
        plt.plot((wavelength), normalized_data, 'ro', label='Bandpass Data')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized Image Counts')
        plt.title('Signal Vs Wavelength (Along  425 nm, pixel # 1500)')
        plt.legend(loc='best')
    
        plt.figure(2)
        plt.plot(wavelength,final_fit-normalized_data,'k.', label= 'fit- data')
        plt.ylim(-0.06, 0.06)
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Residue')
        plt.title('Difference Between Fit and Actual Data (Along 425 nm, pixel # 1500)')
        plt.legend(loc='best')
        plt.show()
        cc

    gauss_fit_params_all1 = np.array(gauss_fit_params_all1)
    gauss_fit_params_all2 = np.array(gauss_fit_params_all2)
    dframe_gauss = pd.DataFrame({'A1' : gauss_fit_params_all1[:, 0],
                                 'CW1' : gauss_fit_params_all1[:, 1],
                                 'Sigma1' : gauss_fit_params_all1[:, 2],
                                 'A2' : gauss_fit_params_all2[:, 0],
                                 'CW2' : gauss_fit_params_all2[:, 1],
                                 'Sigma2' : gauss_fit_params_all2[:, 2]
                                })
    dframe_gauss.to_csv(save_dir +'/' + 'Gaussian_fit_params_' + str(wavelen_val) + '_nm.csv')


def main():
    """
    Main function begins from here. It call all the required functions
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Spectral_Band_Pass_Data'
    channel = 'UV'
    channel_data = os.path.join(data_dir, channel)
    all_wavelengths_file = sorted([each for each in os.listdir(channel_data)
                                   if each.endswith('_nm')])
    for files in range(6, len(all_wavelengths_file)):
        save_dir = os.path.join(channel_data, all_wavelengths_file[files])
        wavelen_val = all_wavelengths_file[files].strip('_nm')

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        bp_data_path = os.path.join(channel_data, all_wavelengths_file[files],
                                    r'bandpass_lookup_table')
        bandpass_data = np.genfromtxt(os.path.join(bp_data_path,
                                                   'bandpass_radiance_normalized_V1.csv'),
                                      delimiter=',')
        bandpass_wvl = np.genfromtxt(os.path.join(bp_data_path,
                                                  'bandpass_wvl_V1.csv'),
                                     delimiter=',')
        fit_double_gaussian_func(bandpass_data, bandpass_wvl, save_dir, wavelen_val)


if __name__ == "__main__":
    main()
