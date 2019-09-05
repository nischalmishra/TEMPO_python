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
import matplotlib.ticker as tkr
import pandas as pd
#pd.options.display.float_format = '${:,.3f}'.format
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


def parse_wavelength_data(wavelen, wavelength_data):
    """ Parse the wavelenght file. The wavelength file was created manually
    from the spectrometer log file. It was a painful process
    """
    dframe = pd.read_csv(wavelength_data, sep=',')
    wavelengths = dframe[wavelen].dropna(how='all')
    return wavelengths


def gauss_function(x, a, x0, sigma):
    """ Fit Gaussian function
    """

    return  a *np.exp(-(x-x0)**2/(2*sigma**2))

def gauss_function_double(wavelength, amplitude, center, width):
    amplitude1 = amplitude
    center1= center
    width1= width+10
    
    return  amplitude * np.exp(-(wavelength-center)**2 / width) + \
            amplitude1 * np.exp(-(wavelength-center1)**4 / width1)
            


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

    #print(FWHM)
    return FWHM


def fit_double_gaussian_func(wavelength, radiance_data, wavelen_val):
    """ Fit double Gaussian Function
    """
    wavelen_val = float(wavelen_val.strip('_')[0:3])

    radiance_data = radiance_data[:, 4]
    dframe = pd.DataFrame()
    dframe['wavelength'] = wavelength
    dframe['radiance_data'] = radiance_data
    dframe = dframe.groupby('wavelength', as_index=False).mean()
    dframe.reset_index(drop=True)
    wavelength = dframe['wavelength']
    radiance_data = dframe['radiance_data']
    radiance_data = radiance_data / np.max(radiance_data)
    width = 2
    center = wavelen_val
    amp = 1
    init_vals = [amp, center, width]
    best_vals, covarr = curve_fit(gauss_function, wavelength, radiance_data, p0=[init_vals])

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
    
    
    

def fit_gaussian_func(wavelength, radiance, wavelen_val, file_dir):
    """ Function to fit the Gaussian fit to the wavelength Vs. Amplitude Data
    """
    save_dir = os.path.join(file_dir, 'Pixel_to_wavelen_map')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)      
    wavelen_val = float(wavelen_val.strip('_nm'))
    num_data = radiance.shape[1]    
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength
    dframe_gauss = pd.DataFrame()
    pixel_to_wavelen_map = []
    for i in range(4, num_data-7):
     # waveln example : 675_nm represents the header in the wavelength file
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
         fitted_data_gauss = fitted_data_gauss/np.max(fitted_data_gauss)     
         pixel_to_wavelen_map.append(best_vals)
         print(best_vals)
         
         plt.plot(sampled_pixel_wvl, fitted_data_gauss, 'g', label=' Gaussian fit')
         plt.plot(wavelength, radiance_data, 'ro--', label='Wavelength averaged data')
         plt.grid(True, linestyle=':')
         plt.xlabel('Wavelength (nm)')
         plt.ylabel('Image Counts (DN)')
         plt.title('Signal Vs Wavelength (736 nm, pixel # 500)')
         plt.legend(loc='best')
    
        #plt.ylim(0, 1)
#        #plt.yscale('log')
         plt.show()

    pixel_to_wavelen_map = np.array(pixel_to_wavelen_map)
    dframe_gauss = pd.DataFrame({'A' : pixel_to_wavelen_map[:, 0],
                                 'CW' : pixel_to_wavelen_map[:, 1],
                                 'Sigma' : pixel_to_wavelen_map[:, 2]                                                            
                                    })
    dframe_gauss.to_csv(save_dir +'/' + 'Gaussian_fit_params_' + str(wavelen_val) + '_nm.csv')
    
def fit_linear_interpolation(wavelength, radiance, wavelen_val):
    """ Function to fit the linear itnerpolation to the wavelength Vs. Amplitude Data
    """
    
    wavelen_val = float(wavelen_val.strip('_nm'))
    num_data = radiance.shape[1]
    dframe1 = pd.DataFrame()
    dframe1['wavelength'] = wavelength

    FWHM_all = []
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\All_FWHM'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    csv_name = 'FWHM_' + str(wavelen_val) + 'nm.csv'
    save_FWHM = save_dir +'/'+csv_name
    csv_name_bandpass = 'Bandpass_' + str(wavelen_val) + 'nm.csv'
    pixel_indices = [1000, 1004, 1404, 2004]
    for i in pixel_indices:
        sub_path = os.path.join(save_dir, 'spectral_bandpass_' + str(i))
        if not os.path.exists(sub_path):
            os.makedirs(sub_path)        
        #print(i)
        #i = 2044
        dframe2 = pd.DataFrame()
        dframe3 = pd.DataFrame()
        radiance_data = radiance[:, i]
        dframe2['rad'] = radiance_data
        dframe2.astype(float)
        dframe3['rad'] = dframe2['rad']
        dframe3['wavelength'] = dframe1['wavelength'].astype(float)
        dframe3 = dframe3.groupby('wavelength', as_index=False).mean()
        dframe3.reset_index(drop=True)
        

         # waveln example : 675_nm represents the header in the wavelength file
        wavelength = dframe3['wavelength']
        radiance_data = dframe3['rad']
        radiance_data = radiance_data/np.max(radiance_data)
        sampled_wavelength = np.arange(min(wavelength), max(wavelength), 0.00010)
        fit_params = interp1d(wavelength, radiance_data, kind='linear')
        fitted_val = fit_params(sampled_wavelength)
        FWHM_val = full_width_half_max(sampled_wavelength, fitted_val)
        #FWHM_val = FWHM_val[0]
        FWHM_all.append(FWHM_val)
        dframe3['rad'] = dframe3['rad']/np.max(dframe3['rad'])
        dframe3.to_csv(sub_path + '/' + csv_name_bandpass)

        plt.figure()
        plt.plot(wavelength, radiance_data,'ro--', label='Original Data')
        plt.plot(sampled_wavelength, fitted_val, 'b', label='Interpolated Data' )
        plt.grid(True, linestyle  = ':')
        plt.xlabel ('Wavelength (nm)')
        plt.ylabel('Image Counts (DN)')
        plt.title('Signal Vs Wavelength (736 nm, pixel # 500)')
        plt.legend(loc='best')
        plt.show()
#
#        #plt.yscale('log')
#        plt.show(block=False)
###    print(FWHM_all)
    
    FWHM_all = np.array(FWHM_all)
    #print(FWHM_all.shape)

    
    

    #np.savetxt(save_FWHM, FWHM_all, delimiter=",")




#    failed_pixels = len(FWHM_all[FWHM_all > 0.6])
#    text1 = '*** Req. not met on ' + str(failed_pixels)+' out of '+ str(len(FWHM_all)) +' pixels' + ' ('+ str(round(100* failed_pixels/len(FWHM_all), 2)) +'%)'
#
#
#    plt.figure(figsize=(10, 4))
#    plt.plot(FWHM_all, 'b*', markersize=1.5, label='FWHM')
#    plt.xticks(np.arange(0, 2048, 200))
#    plt.ylim(0.55, 0.63)
#    plt.axhline(y=0.6, color='r', linestyle='--', label='TEMPO Requirement')
#    plt.xlabel('Spatial Pixel Index (#)')
#    plt.ylabel('FWHM (nm)')
#    plt.suptitle('TEMPO Spectral Bandpass FWHM Over the Entire FOV (' + str(wavelen_val)+ ' nm)')
#    plt.title(text1, backgroundcolor='yellow', size=10, color='red')
#    plt.grid(True, linestyle=':')
#    plt.legend(loc='best')
#    plt.show()


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


def plot_integrated_counts(signal, wavelength, wavelen_val):
    """ Plot integrated counts for each line
    """
    wavelen_val = wavelen_val.strip('_nm')
    plt.figure(figsize=(10, 4))
    ax = plt.gca()
    #ax.get_yaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().set_major_formatter(tkr.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.plot(wavelength, signal, 'b.', markersize=4.5)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Image Counts (DN)')
    plt.title('Integrated Signal Counts for 22 rows Vs. Wavelength \n (Encompasses all the illuminated rows for the data collection of '
              + str(wavelen_val) + ' nm)')
    plt.grid(True, linestyle=':')
    plt.ylim(5000000, 80000000)
    plt.show()

def plot_FWHM_all_wavelength(data_dir):
    """ Plot FWHM for all spectral bandpass
    """
    FWHM_loc = r'All_FWHM_Gaussian_flat_top'
    FWHM_files = sorted([each for each in os.listdir(os.path.join(data_dir, FWHM_loc))
                         if each.startswith('FWHM_Gaussian')])
    
    plt.figure(figsize=(12, 4))
    for data in FWHM_files:
        
        bandpass_data = np.loadtxt(os.path.join(data_dir, FWHM_loc, data), delimiter=",")
        legends = data[13:18] + ' nm'
        

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
        radiance_file = os.path.join(data_path, r'bandpass_radiance_normalized_second_iteration.csv')
        #integrated_signal_file = os.path.join(data_path, r'integrated_radiance_second_approach.csv')
        print(radiance_file)
        #cc
        radiance_data = np.loadtxt(radiance_file, delimiter=",")
        # Now let's read the wavelength file. The wavelength file is a csv file
        #created from the image log file.
        wavelength_data = os.path.join(data_dir, r'wavelength_spectrum_analyzer.csv')
        wavelength = parse_wavelength_data(wavelen, wavelength_data)
        #integrated_signal = np.loadtxt(integrated_signal_file, delimiter = ",")
         ##the data may have some obvious outliers. Those need to be tossed out

         ##The wavelength data was acquired at finite spectral sampling. Hence
        ## a continuos curve is generated by fitting a Gaussian curve to the data


        fit_linear_interpolation(wavelength, radiance_data, wavelen)
        #plot_bandpass_data(wavelength, radiance_data, wavelen)
        #fit_gaussian_func(wavelength, radiance_data, wavelen, data_dir)
        #fit_double_gaussian_func(wavelength, radiance_data, wavelen)
        #plot_integrated_counts(integrated_signal,wavelength,  wavelen)
    #plot_FWHM_all_wavelength(data_dir)
#
if __name__ == "__main__":
    main()
