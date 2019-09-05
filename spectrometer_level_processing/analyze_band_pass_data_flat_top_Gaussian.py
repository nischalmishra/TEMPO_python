# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:12:07 2018

Code written to verify some of the TEMPO Bandpass requirements.

@author: nmishra
"""

import os
from random import randint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatteri



def parse_wavelength_data(wavelen, wavelength_data):
    """ Parse the wavelength file. The wavelength file was created manually
    from the spectrometer log file. It was a painful process
    """
    dframe = pd.read_csv(wavelength_data, sep=',')
    wavelengths = dframe[wavelen].dropna(how='all')
    return wavelengths


def fit_double_gaussian_func(wavelength, radiance_data, wavelen_val):
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
    mean_val1 = np.mean(wavelength)
    sigma_val1 = np.std(wavelength)
    mean_val2 = mean_val1
    sigma_val2 = sigma_val1
    best_vals, covarr = curve_fit(flat_top_gauss_function, wavelength, radiance_data,
                                  p0=[1, mean_val1, sigma_val1, 1, mean_val2, sigma_val2])
    #sampled_pixel_wvl= np.arange(min(wavelength), max(wavelength), 0.00001)
    fitted_data = flat_top_gauss_function(wavelength, *best_vals)
    fitted_data = fitted_data/np.max(fitted_data)

    sampled_pixel_wvl1 = np.linspace(min(wavelength), max(wavelength), len(fitted_data))
    plt.plot(sampled_pixel_wvl1, fitted_data, 'go-', label='Gaussian fit')
    plt.plot(wavelength, radiance_data, 'ro', label='Wavelength averaged data')
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Image Counts (DN)')
    plt.title('Signal Vs Wavelength (736 nm, pixel # 500)')
    plt.legend(loc='best')
    plt.show()



def rect_function(wavelength, center_wavelen):
    #center_wavelen = center_wavelen-0.005
    y = np.where((wavelength <= center_wavelen+0.8) & (wavelength >= center_wavelen-0.8), 1, 0.0)
    return y

def gauss_function(x, a, x0, sigma):
    """ Fit Gaussian function
    """

    return  a *np.exp(-(x-x0)**2/(2*sigma**2))

def flat_top_gauss_function(x, a, x0, sigma):
    """ Flat Top Gaussian function
    """
    return  a *np.exp(-(x-x0)**4/(2*sigma**4))


def super_gaussian_function(W, W0, A, K):
    """
    W0: center wavelength
    W : Wavelength
    A : Asymmetric Parameter. If A is not zero, slit function is asymetric
    k: if K is fixed at 2, slit function is standard Gaussian
    """
    K=4
    A = 0.03
    Sigma = np.std(W)
    print(Sigma/3.8)
    W0 = 541.7893654   
    return np.exp(-abs((W-W0)/(0.3+np.sign(W-W0)*A))**K)


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


def fit_super_gaussian_func(wavelength, radiance, wavelen_val):
    num_data = radiance.shape[1]
    #print(num_data)
    #cc
    dframe_gauss = pd.DataFrame()   
    estimates_super_gaussian = [ ]
    
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength

    FWHM_all = []
    for i in range(4, num_data-7):
        #print(i)
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i ]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).mean()
        dframe3.reset_index(drop=True)
       # print(dframe['wavelength'].head())
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data / np.max(radiance_data)
        mean_val = np.mean(wavelength)       
        
        #(W, W0, A, K)        
        best_vals, covarr = curve_fit(super_gaussian_function, wavelength, radiance_data)
                                      
        sampled_pixel_wvl = np.arange(min(wavelength), max(wavelength), 0.001)
        fitted_data_super_gauss = super_gaussian_function(sampled_pixel_wvl, *best_vals)
        
        plt.plot(wavelength, radiance_data, 'ro')
        plt.plot(sampled_pixel_wvl, fitted_data_super_gauss,'b', alpha=0.9)
        plt.show()
        

        estimates_super_gaussian.append(best_vals)
        #print(estimates_gauss)
 
 
    
    estimates_super_gaussian = np.array(estimates_super_gaussian)
    dframe_gauss = pd.DataFrame({'A1' : estimates_super_gaussian[:, 0],
                                 'W1' : estimates_super_gaussian[:, 1],
                                 'Sigma1' : estimates_super_gaussian[:, 2]                                                                                             
                                    })
   
    
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\All_FWHM_Gaussian_flat_top_UV_visible_same'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    csv_name = 'FWHM_Gaussian' + str(wavelen_val) + 'nm.csv'
    gauss_fit_params = 'Params_Gauss' + str(wavelen_val) + 'nm.csv'
    save_FWHM = save_dir + '/' + csv_name    
    np.savetxt(save_FWHM, FWHM_all, delimiter=",")
    dframe_gauss.to_csv(save_dir +'/' + gauss_fit_params)



def fit_gaussian_func(wavelength, radiance_data, wavelen_val):
    """ Function to fit the Gaussian fit to the wavelength Vs. Amplitude Data
    """
    wavelen_val = float(wavelen_val.strip('_')[0:3])


    radiance_data = radiance_data[:, 500]
    dframe = pd.DataFrame()
    dframe['wavelength'] = wavelength
    dframe['radiance_data'] = radiance_data
    dframe = dframe.groupby('wavelength', as_index=False).mean()
    dframe.reset_index(drop=True)


     # waveln example : 675_nm represents the header in the wavelength file
    wavelength = dframe['wavelength']
    radiance_data = dframe['radiance_data']
    radiance_data = radiance_data / np.max(radiance_data)
    #n = len(wavelength)
   # mean_val = sum(radiance_data * wavelength)/n
    #sigma_val =  np.sqrt(sum(radiance_data * (wavelength - mean_val)**2) / n)
    mean_val = np.mean(wavelength)
    sigma_val = np.std(wavelength)
    #print(mean_val, sigma_val)

    best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data,
                                  p0=[1, mean_val, sigma_val], W=0.2)
    #sampled_pixel_wvl= np.arange(min(wavelength), max(wavelength), 0.00001)
    fitted_data = gauss_function(wavelength, *best_vals)

    rect_data = rect_function(wavelength, wavelen_val)

    fitted_data = fitted_data/np.max(fitted_data)
    new_fit = np.convolve(rect_data, fitted_data)
    new_fit = new_fit/max(new_fit)
    sampled_pixel_wvl = np.linspace(min(wavelength), max(wavelength), len(new_fit))
    fitted_data = gauss_function(sampled_pixel_wvl, *best_vals)
    fitted_data = fitted_data/np.max(fitted_data)
    final_fit = new_fit+fitted_data
    final_fit = final_fit/np.max(final_fit)
    #plt.plot(sampled_pixel_wvl, final_fit, 'bo-', label='Flat Top Gaussian')
    #full_width_half_max(wavelength, radiance_data)

    plt.plot(sampled_pixel_wvl, final_fit, '*-', color='lime',
             label=' Flat Top Gaussian + Standard Gaussianfit')
    plt.plot(wavelength, radiance_data, 'ro', label='Peak Normalized Spectral Data')
   # plt.plot(wavelength, fitted_data, 'b*-', label ='Flat Top Gaussian Fit')
    #plt.plot(wavelength, rect_data, 'bo--', label= 'Rectangle Function')
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Normalized Image Counts (DN)', fontsize=12)
    plt.title('Signal Vs Wavelength (320 nm, pixel # 1000)', fontsize=12)
    plt.legend(loc='best')
    plt.show()
    dframe1 = pd.DataFrame()
    dframe1['wvl'] = wavelength
    dframe1['radiance_original'] = radiance_data
    dframe = pd.DataFrame()
    dframe['wvl'] = sampled_pixel_wvl
    dframe['radiance_fit'] = final_fit
    # Inner join to do find the same wavelength
    data = pd.merge(dframe, dframe1, how='inner', on=['wvl'])
    residue = data['radiance_fit'] - data['radiance_original']
    #mean_diff = 'Mean Diff = ', str(round(np.mean(fitted_data-radiance_data)))
    #std_dev = 'Stdev. = ', str(round(np.std(fitted_data-radiance_data)))
    #MSE = 'MSE = ', +str(round(MSE,2))
#    text1 = 'MSE = ' +str(round(MSE,2)) +\
#            '\nMean Diff = ' + str(round(np.mean(fitted_data-radiance_data), 2)) +\
#            '\n Stdev  = ' + str(round(np.std(fitted_data-radiance_data), 2))


    plt.plot(data['wvl'], residue, 'mo--', label='Flat Top Gaussian')
   # plt.plot(wavelength,fitted_data-radiance_data,'ko--', label='Standard Gaussian')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Resdiue [Fit-Data]', fontsize=12)
    plt.title('Residue of a signal fit Vs. Wavelength (675 nm, pixel #)', fontsize=12)
    plt.grid(True, linestyle=':')
    #plt.legend(loc='best')
    plt.show()

    #plt.ylim(0, 1)
    #plt.yscale('log')


def fit_gaussian_func_two_combo(wavelength, radiance_data, wavelen_val):
    """ Function to fit the Gaussian fit to the wavelength Vs. Amplitude Data
    """
    std_gauss = [ ]
    flat_top_gauss =  [ ]
    wavelen_val = float(wavelen_val.strip('_')[0:3])
    num_data = radiance_data.shape[1]
    radiance_data = radiance_data[:, 2040]
    dframe = pd.DataFrame()
    dframe['wavelength'] = wavelength
    dframe['radiance_data'] = radiance_data
    dframe = dframe.groupby('wavelength', as_index=False).mean()
    dframe.reset_index(drop=True)

     # waveln example : 675_nm represents the header in the wavelength file
    wavelength = dframe['wavelength']
    radiance_data = dframe['radiance_data']
    radiance_data = radiance_data / np.max(radiance_data)
    #n = len(wavelength)
   # mean_val = sum(radiance_data * wavelength)/n
    #sigma_val =  np.sqrt(sum(radiance_data * (wavelength - mean_val)**2) / n)
    mean_val = np.mean(wavelength)
    sigma_val = np.std(wavelength)
    #print(mean_val, sigma_val)

    best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data,
                                  p0=[1, mean_val, sigma_val])
    #sampled_pixel_wvl= np.arange(min(wavelength), max(wavelength), 0.00001)
    fitted_data_gauss = gauss_function(wavelength, *best_vals)
    print(best_vals)
  

    best_vals1, covarr1 = curve_fit(flat_top_gauss_function, wavelength, radiance_data,
                                    p0=[1, mean_val, sigma_val])

    print(best_vals1)
    
    fitted_data_gauss_flat_top = flat_top_gauss_function(wavelength, *best_vals1)
    #MSE = np.mean((radiance_data-gauss_function(wavelength, *best_vals)**2))

    #print ("Mean Squared Error: ", MSE)


    fitted_data = fitted_data_gauss + fitted_data_gauss_flat_top
    fitted_data = fitted_data/np.max(fitted_data)

    #sampled_pixel_wvl = np.linspace(min(wavelength), max(wavelength), len(new_fit))



    #plt.plot(sampled_pixel_wvl, final_fit, 'bo-', label='Flat Top Gaussian')
    #full_width_half_max(wavelength, radiance_data)


    plt.plot(wavelength, radiance_data, 'ro', label='Peak Normalized Spectral Data')
    plt.plot(wavelength, fitted_data, '*-', color='purple',
             label=' Flat Top Gaussian + Standard Gaussian fit')
   # plt.plot(wavelength, fitted_data, 'b*-', label ='Flat Top Gaussian Fit')
    #plt.plot(wavelength, rect_data, 'bo--', label= 'Rectangle Function')
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Normalized Image Counts (DN)', fontsize=12)
    plt.title('Signal Vs Wavelength (675 nm, pixel #' + str(num_data-7) + ')',
              fontsize=12)
    plt.legend(loc='best')
    plt.show()

    plt.figure()

    plt.plot(wavelength, fitted_data-radiance_data, 'ko--', label='Standard Gaussian')

    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Resdiue [Fit-Data]', fontsize=12)
    plt.title('Residue of a signal fit Vs. Wavelength (675 nm, pixel # 500)', fontsize=12)
    plt.grid(True, linestyle=':')
    #plt.legend(loc='best')
    plt.show()

    #plt.ylim(0, 1)
    #plt.yscale('log')


def plot_bandpass_data(wavelength, radiance_data, wavelen_val):
    """ Plot few curves of spectral bandpass data
    """
    # waveln example : 675_nm represents the header in the wavelength file
    #print(wavelen_val)
    wavelen_val = wavelen_val.strip('_nm')
    col_index = ([randint(5, 2000) for p in range(0, 3)])
    plt.figure()
    for i in col_index:
        plt.plot(wavelength, radiance_data[:, i], '.', label='spatial index = '+ str(i))
    plt.title('TEMPO Spectral Bandpass Function, CW = '+ str(wavelen_val)+ ' nm\n'
              + '(Linear Scale)')
    plt.ylabel('Signal Counts (DN)')
    plt.xlabel('Wavelength (nm)')
    #plt.xticks(np.arange(np.min(round(wavelength, 1)), np.max(round(wavelength, 1)), 0.5))
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    #plt.yscale('log')
    plt.show()



def calculate_FWHM_all(wavelength, radiance, wavelen_val):
    wavelen_val = float(wavelen_val.strip('_nm'))
    num_data = radiance.shape[1]
    #print(num_data)
    #cc
    dframe_gauss = pd.DataFrame()
    estimates_gauss = []
    estimates_flat_top_gauss = [ ]
    
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength

    FWHM_all = []
    for i in range(4, num_data-7):
        #print(i)
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i ]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).mean()
        dframe3.reset_index(drop=True)
       # print(dframe['wavelength'].head())
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data / np.max(radiance_data)
        mean_val = np.mean(wavelength)
        sigma_val = np.std(wavelength)
        best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data,
                                      p0=[1, mean_val, sigma_val])
        sampled_pixel_wvl = np.arange(min(wavelength), max(wavelength), 0.001)
        fitted_data_gauss = gauss_function(sampled_pixel_wvl, *best_vals)

        best_vals1, covarr1 = curve_fit(flat_top_gauss_function, wavelength, radiance_data,
                                        p0=[1, mean_val, sigma_val])

#        print(best_vals)
#        print(best_vals1)
#        cc
        estimates_gauss.append(best_vals)
        estimates_flat_top_gauss.append(best_vals1)
        #print(estimates_gauss)
        
        fitted_data_gauss_flat_top = flat_top_gauss_function(sampled_pixel_wvl, *best_vals1)
        fitted_data = fitted_data_gauss + fitted_data_gauss_flat_top
        fitted_data = fitted_data/np.max(fitted_data)
        FWHM_val = full_width_half_max(sampled_pixel_wvl, fitted_data)
        del dframe2
        del dframe3

        #FWHM_val = FWHM_val[0]
        FWHM_all.append(FWHM_val)
    FWHM_all = np.array(FWHM_all)
    estimates_gauss = np.array(estimates_gauss)
    estimates_flat_top_gauss = np.array(estimates_flat_top_gauss)
    dframe_gauss = pd.DataFrame({'A1' : estimates_gauss[:, 0],
                                 'W1' : estimates_gauss[:, 1],
                                 'Sigma1' : estimates_gauss[:, 2],
                                 'A2': estimates_flat_top_gauss[:, 0],
                                 'W2' : estimates_flat_top_gauss[:, 1],
                                 'Sigma2' : estimates_flat_top_gauss[:, 2]                                                            
                                    })
   
    
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\All_FWHM_Gaussian_flat_top'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    csv_name = 'FWHM_Gaussian' + str(wavelen_val) + 'nm.csv'
    gauss_fit_params = 'Params_Gauss' + str(wavelen_val) + 'nm.csv'
    save_FWHM = save_dir + '/' + csv_name    
    np.savetxt(save_FWHM, FWHM_all, delimiter=",")
    dframe_gauss.to_csv(save_dir +'/' + gauss_fit_params)
    cc
    




def calculate_FWHM_all_only_gaussian(wavelength, radiance, wavelen_val):
    wavelen_val = float(wavelen_val.strip('_nm'))
    num_data = radiance.shape[1]
    #print(num_data)
    #cc
    dframe_gauss = pd.DataFrame()
    estimates_gauss = []
    estimates_flat_top_gauss = [ ]
    
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength

    FWHM_all = []
    for i in range(4, num_data-7):
        #print(i)
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i ]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).mean()
        dframe3.reset_index(drop=True)
       # print(dframe['wavelength'].head())
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data / np.max(radiance_data)
        mean_val = np.mean(wavelength)
        sigma_val = np.std(wavelength)
        best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data,
                                      p0=[1, mean_val, sigma_val])
        sampled_pixel_wvl = np.arange(min(wavelength), max(wavelength), 0.001)
        fitted_data_gauss = gauss_function(sampled_pixel_wvl, *best_vals)

#        print(best_vals)
#        print(best_vals1)
#        cc
        estimates_gauss.append(best_vals)
        
        del dframe2
        del dframe3

        #FWHM_val = FWHM_val[0]
       
   
    estimates_gauss = np.array(estimates_gauss)
    estimates_flat_top_gauss = np.array(estimates_flat_top_gauss)
    dframe_gauss = pd.DataFrame({'A1' : estimates_gauss[:, 0],
                                 'W1' : estimates_gauss[:, 1],
                                 'Sigma1' : estimates_gauss[:, 2]                                                                                          
                                    })
   
    
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\All_FWHM_only_Gaussian'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    csv_name = 'FWHM_Gaussian_' + str(wavelen_val) + 'nm.csv'
    gauss_fit_params = 'Params_Gauss_' + str(wavelen_val) + 'nm.csv'
    save_FWHM = save_dir + '/' + csv_name    
    np.savetxt(save_FWHM, FWHM_all, delimiter=",")
    dframe_gauss.to_csv(save_dir +'/' + gauss_fit_params)
    cc



def plot_FWHM_all_wavelength(data_dir):
    """ Plot FWHM for all spectral bandpass
    """
    FWHM_loc = r'All_FWHM'
    FWHM_files = sorted([each for each in os.listdir(os.path.join(data_dir, FWHM_loc))
                         if each.endswith('nm.csv')])
    plt.figure(figsize=(12, 4))
    for data in FWHM_files:
        #print(data)

        bandpass_data = np.loadtxt(os.path.join(data_dir, FWHM_loc, data), delimiter=",")
        legends = data[5:10] + ' nm'
        #print(legends)

        plt.plot(bandpass_data, '.', markersize=8, label=legends)

    plt.axhline(y=0.6, color='r', linestyle='--', label='TEMPO Req.')
    plt.title('FWHM Vs. Spatial Index for Different Wavelengths')
    plt.xlabel('Spatial Index(#)')
    plt.ylabel('FWHM (nm)')
    plt.xlim(0, 2100)
    plt.ylim(0.55, 0.65)
    plt.grid(True, linestyle=':')
    plt.legend(loc='upper center', bbox_to_anchor=(0.55, 0.3), ncol=4)
    plt.show()


    # now plot these as a function of wavelengths too

    plt.figure(figsize=(8, 4))
    spatial_pixels_begin = [p for p in range(0, 700)]
    spatial_pixels_middle = [p for  p in range(700, 1400)]
    spatial_pixels_end = [p for p in range(1337, 2037)]
    print(len(spatial_pixels_begin))
    #cc
    all_wavelength = []
    all_pixels_begin = []
    all_pixels_middle = []
    all_pixels_end = []
    #data_all = []

    for data in FWHM_files:
        bandpass_data = np.loadtxt(os.path.join(data_dir, FWHM_loc, data), delimiter=",")
        wavelength = data[5:10]
        all_wavelength. append([float(wavelength)])
        all_pixels_begin.append(bandpass_data[spatial_pixels_begin])
        all_pixels_middle.append(bandpass_data[spatial_pixels_middle])
        all_pixels_end.append(bandpass_data[spatial_pixels_end])

   # let plot the data now
    wavelength = all_wavelength
    all_pixels_begin = np.array(all_pixels_begin)
    all_pixels_middle = np.array(all_pixels_middle)
    all_pixels_end = np.array(all_pixels_end)
    print(all_pixels_end.shape)
    #print(wavelength.shape)
    legends = ['Spatial index: 0 to 700', 'Spatial index: 701 to 1400',
               'Spatial index: 1401 onwards']


    for i in range(0, all_pixels_begin.shape[1]):
        #print(i)
        data_begin = all_pixels_middle[:, i]
        #data_middle = all_pixels_middle[:, i]
        #data_end= all_pixels_end[:, i]

        plt.plot(wavelength, data_begin, 'o-', linewidth=1, label='Spatial index: 0 to 700')
        #plt.plot(wavelength, data_middle,'mo--', linewidth=0.5, label='Spatial index: 701 to 1400')
        #plt.plot(wavelength, data_end,'g*--', linewidth=0.5, label='Spatial index: 1401 onwards')

    plt.axhline(y=0.6, color='r', linestyle='-', linewidth=2, label='TEMPO Requirement')
    plt.title('FWHM Vs. Wavelengths for different spatial pixels\n' +'(' + legends[1]+')')
    plt.xlabel('Wavelength(nm)')
    plt.ylabel('FWHM (nm)')
    plt.xlim(290, 800)
    plt.ylim(0.53, 0.65)
    plt.grid(True, linestyle=':')
    #plt.legend([legends[0]], loc='best')
    plt.show()



def main():
    """ Main function begins from here
    """
    data_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass'
    channel = 'Visible'
    data_channel = os.path.join(data_dir, channel)
    data_wavelengths = sorted([each for each in os.listdir(data_channel)
                               if each.endswith('727.5_nm')])
    for wavelen in data_wavelengths:
        # first let's read the spectral radiance data
        data_path = os.path.join(data_channel, wavelen, r'saved_quads\processed_image')
        radiance_file = os.path.join(data_path, r'bandpass_radiance_normalized.csv')
        #integrated_signal_file = os.path.join(data_path, 
                                               #r'integrated_radiance_second_approach.csv')
        print(radiance_file)
        #cc
        radiance_data = np.loadtxt(radiance_file, delimiter=",")
        # Now let's read the wavelength file. The wavelength file is a csv file
        #created from the image log file.
        wavelength_data = os.path.join(data_dir, r'wavelength_spectrum_analyzer.csv')
        wavelength = parse_wavelength_data(wavelen, wavelength_data)
        #integrated_signal = np.loadtxt(integrated_signal_file, delimiter=",")
         ##the data may have some obvious outliers. Those need to be tossed out
         ##The wavelength data was acquired at finite spectral sampling. Hence
        ## a continuos curve is generated by fitting a Gaussian curve to the data
        #fit_linear_interpolation(wavelength, radiance_data, wavelen)

       # print(wavelen)
        #fit_gaussian_func(wavelength, radiance_data, wavelen)
        #fit_gaussian_func_two_combo(wavelength, radiance_data, wavelen)
        calculate_FWHM_all_only_gaussian(wavelength, radiance_data, wavelen)
        #fit_super_gaussian_func(wavelength, radiance_data, wavelen)
        #calculate_FWHM_all(wavelength, radiance_data, wavelen) # flat top Gaussian
        #plot_integrated_counts(integrated_signal,wavelength,  wavelen)
#    plot_FWHM_all_wavelength(data_dir)

if __name__ == "__main__":
    main()
