# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 13:45:49 2018

@author: nmishra
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
np.set_printoptions(threshold=np.inf)
#pylint: disable = E1101


def read_radiance_data():
    """ Read in the high res. radiance data from SAO
    """
    data_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass'
    dframe = pd.read_csv(data_dir +'/' + 'High_res_radiance_data.csv', delimiter=',')

    plt.plot(dframe['Wavelength'], dframe['Response'], 'k.', markersize=1,
             label='SAO 2010 Solar Irradiance Spectrum')
    plt.grid(True, linestyle=':')
    plt.title('Solar Irradiance Vs. Wavelength', fontsize=12)
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Solar Irradiance (' + r'$Wm^{-2}nm^{-1}$' +')', fontsize=12)
    plt.legend(loc='best')
    plt.savefig(data_dir +'/'+'SAO_Solar_Irradiance.png', dpi=100)
    plt.close('all')
    return dframe


def gauss_function(x_val, a_val, x0_val, sigma_val):
    """ Fit Gaussian function
    """

    return  a_val *np.exp(-(x_val-x0_val)**2/(2*sigma_val**2))

##### WRITE THIS FUNCTION FOR GAUSSIAN OLNY!!!!!!!!!!!!!!!!!!!!!!
def perform_spectral_interpolation_only_gaussian(gaussian_data):
    """ perform spectral intepolation of Gaussian fits"""
    dframe = pd.DataFrame()

    wavelength1 = gaussian_data[:, -1]

    sampled_wavelength1 = np.arange(min(wavelength1), max(wavelength1), 2)
    a1_val = gaussian_data[:, 0]
    sigma1 = gaussian_data[:, 1]

    # A1 first
    fit_params_a1 = interp1d(wavelength1, a1_val, kind='linear')
    fitted_val_a1 = fit_params_a1(sampled_wavelength1)
    # Now A2


    # Now Sigma1
    fit_params_sigma1 = interp1d(wavelength1, sigma1, kind='linear')
    fitted_val_sigma1 = fit_params_sigma1(sampled_wavelength1)

    # Now Sigma2


#    plt.plot(wavelength1, Sigma1, 'bo')
#    plt.plot(sampled_wavelength1, fitted_val_Sigma1, 'ro--', markersize=3)
#    plt.grid(True, linestyle=':')
#    plt.show()
    dframe = pd.DataFrame({'W1' : sampled_wavelength1,
                           'A1' : fitted_val_a1,
                           'Sigma1' : fitted_val_sigma1,
                          })

    return dframe.round(3)



def perform_spectral_interpolation(gaussian_data):
    """ perform spectral intepolation of Gaussian fits"""

    dframe = pd.DataFrame()
    wavelength1 = gaussian_data[:, -1]

    sampled_wavelength1 = np.arange(min(wavelength1), max(wavelength1), 2)
    wavelength2 = gaussian_data[:, -1]
    sampled_wavelength2 = np.arange(min(wavelength2), max(wavelength2), 2)
    a1_val = gaussian_data[:, 0]
    a2_val = gaussian_data[:, 1]
    sigma1 = gaussian_data[:, 2]
    sigma2 = gaussian_data[:, 3]

    # A1 first
    fit_params_a1 = interp1d(wavelength1, a1_val, kind='linear')
    fitted_val_a1 = fit_params_a1(sampled_wavelength1)
    # Now A2
    fit_params_a2 = interp1d(wavelength2, a2_val, kind='linear')
    fitted_val_a2 = fit_params_a2(sampled_wavelength2)

    # Now Sigma1
    fit_params_sigma1 = interp1d(wavelength1, sigma1, kind='linear')
    fitted_val_sigma1 = fit_params_sigma1(sampled_wavelength1)

    # Now Sigma2
    fit_params_sigma2 = interp1d(wavelength2, sigma2, kind='slinear')
    fitted_val_sigma2 = fit_params_sigma2(sampled_wavelength2)


#    plt.plot(wavelength1, Sigma1, 'bo')
#    plt.plot(sampled_wavelength1, fitted_val_Sigma1, 'ro--', markersize=3)
#    plt.grid(True, linestyle=':')
#    plt.show()
    dframe = pd.DataFrame({'W1' : sampled_wavelength1,
                           'W2' : sampled_wavelength2,
                           'A1' : fitted_val_a1,
                           'A2' : fitted_val_a2,
                           'Sigma1' : fitted_val_sigma1,
                           'Sigma2' : fitted_val_sigma2,
                          })

    return dframe.round(3)


def flat_top_gaussian(a1_val, a2_val, sigma1, sigma2, w1_val, w2_val, w_val):
    """
    Fit the flat top Gaussian FUnction
    """
    gauss1 = a1_val * np.exp(-(w_val - w1_val)**4/(2 * sigma1**2))
    gauss2 = a2_val * np.exp(-(w_val - w2_val)**4/(2 * sigma2**4))
    sum_gauss = gauss1 + gauss2
    return sum_gauss


def gauss_function_only(a1_val, sigma1, w1_val, w_val):
    """ Fit Gaussian function
    """

    return  a1_val *np.exp(-(w_val-w1_val)**2/(2*sigma1**2))

def create_spectral_bandpass_only_gaussian(dframe, radiance, file_path):
    """ Creates the spectral bandpass from Gaussian fit parameters
    """
    print(radiance)
    save_dir = os.path.join(file_path, 'spectral_bandpass_1400')
    print(save_dir)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    dframe1 = pd.DataFrame()
    for i in range(0, len(dframe['W1'])):
    #plt.plot(radiance['Wavelength'], radiance['Response']/np.max(radiance['Response']),
     #         'k--', markersize=2, label='SAO 2100 Solar Irradiance Spectrum')

   # for i in range(0, 5):
        a1_val = dframe['A1'][i]
        sigma1 = dframe['Sigma1'][i]
        w1_val = dframe['W1'][i]

        lower_range = w1_val - 1.92
        upper_range = w1_val + 1.92

        wavelens = np.arange(lower_range, upper_range, 0.01)
        #wavelens = ran
        bandpass = [gauss_function_only(a1_val, sigma1, w1_val, wavelens)
                    for wavelens in np.arange(lower_range, upper_range, 0.01)]

        dframe1['Wavelength'] = wavelens
        dframe1['Response'] = bandpass/np.max(bandpass)
        #dframe1 = dframe1.round(3)
        dframe1.round(4).to_csv(save_dir + '/' + 'bandpass_' + str(round(w1_val, 2))+'_nm.csv')
        plt.plot(wavelens, bandpass/np.max(bandpass), 'r.--')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized Spectral Response')
        plt.title('TEMPO Spectral Bandpass (WL = ' + str(round(w1_val, 2)) + ' nm)')
        plt.ylim(0, 1.1)
        plt.xlim(lower_range, upper_range)
        #plt.show()
   # plt.show()
        # Now let us save the spectral bandpass data and spectral bandpass plot
        plt.savefig(save_dir + '/' + 'bandpass_' + str(round(w1_val, 2))+'_nm.png', dpi=100)
        plt.close('all')



def create_spectral_bandpass(dframe, radiance, file_path):
    """ Creates the spectral bandpass from Gaussian fit parameters
    """

    save_dir = os.path.join(file_path, 'spectral_bandpass_1400')
    print(save_dir)
    print(radiance)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    dframe1 = pd.DataFrame()
    for i in range(0, len(dframe['W1'])):
    #plt.plot(radiance['Wavelength'], radiance['Response']/np.max(radiance['Response']),
              #'k--', markersize=2, label='SAO 2100 Solar Irradiance Spectrum')

   # for i in range(0, 5):
        a1_val = dframe['A1'][i]
        a2_val = dframe['A2'][i]
        sigma1 = dframe['Sigma1'][i]
        sigma2 = dframe['Sigma2'][i]
        w1_val = dframe['W1'][i]
        w2_val = dframe['W2'][i]


        lower_range = w1_val - 1.92
        upper_range = w1_val + 1.92

        wavelens = np.arange(lower_range, upper_range, 0.01)
        #wavelens = ran
        bandpass = [flat_top_gaussian(a1_val, a2_val, sigma1, sigma2, w1_val,
                                      w2_val, wavelens)
                    for wavelens in np.arange(lower_range, upper_range, 0.01)]

        dframe1['Wavelength'] = wavelens
        dframe1['Response'] = bandpass/np.max(bandpass)
        #dframe1 = dframe1.round(3)
        dframe1.round(4).to_csv(save_dir + '/' + 'bandpass_' + str(round(w1_val, 2))+'_nm.csv')
        plt.plot(wavelens, bandpass/np.max(bandpass), 'r.--')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized Spectral Response')
        plt.title('TEMPO Spectral Bandpass (WL = ' + str(round(w1_val, 2)) + ' nm)')
        plt.ylim(0, 1.1)
        plt.xlim(lower_range, upper_range)
        #plt.show()
   # plt.show()
        # Now let us save the spectral bandpass data and spectral bandpass plot
        plt.savefig(save_dir + '/' + 'bandpass_' + str(round(w1_val, 2))+'_nm.png', dpi=100)
        plt.close('all')



def  perform_point_interpolation(sub_sample_wvl, sub_sample_rad, center_wv):
    """ This function performs the interpolation of spectral bandpass
    """
    # let us define spectral resolution

    print(center_wv)
    dframe = pd.DataFrame()

    sampled_wvl = np.arange(min(sub_sample_wvl), max(sub_sample_wvl), 2)
    fit_params = interp1d(sub_sample_wvl, sub_sample_rad, kind='slinear')
    fitted_val = fit_params(sampled_wvl)
    dframe['wavelength'] = sampled_wvl
    dframe['rad'] = fitted_val
    return dframe

def create_spectral_bandpass_interpol(interpol_wavelen, interpol_rad, center_wvl,
                                      save_dir):
    """ Creates a spectral bandpass
    """

    save_dir = os.path.join(save_dir, r'look_up_table')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)


    center_wvl1 = np.arange(min(center_wvl), max(center_wvl), 2)




    for j in np.arange(0, interpol_wavelen.shape[1]):
        #print(j)
        dframe = pd.DataFrame()
        wavelen = interpol_wavelen[:, j]

        radiance = interpol_rad[:, j]
        sampled_wvl = np.arange(min(wavelen), max(wavelen), 0.01)
        fit_params = interp1d(wavelen, radiance, kind='slinear')
        fitted_val = fit_params(sampled_wvl)
        #peak_val = np.where(fitted_val==max(fitted_val))[0]
        #print(peak_val)
        #peak_shift = sampled_wvl[peak_val] - CW1[j]


#        if peak_shift >0:
#            sampled_wvl = sampled_wvl - peak_shift
#        elif peak_shift <0:
#            sampled_wvl = sampled_wvl + peak_shift
#        else:
#            sampled_wvl = sampled_wvl
#
#        print(sampled_wvl[peak_val] - CW1[j])

        dframe['Wavelength'] = sampled_wvl
        dframe['Radiance'] = fitted_val
        dframe.round(4).to_csv(save_dir + '/' + 'bandpass_' + \
                               str(round(center_wvl1[j], 2))+'_nm.csv')
        plt.plot(sampled_wvl, fitted_val/np.max(fitted_val), 'g.--')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized Spectral Response')
        plt.title('TEMPO Spectral Bandpass (WL = ' + str(round(center_wvl1[j], 2)) + ' nm)')
        plt.ylim(0, 1.1)
        plt.xlim(min(wavelen), max(wavelen))
        #plt.show()

        # Now let us save the spectral bandpass data and spectral bandpass plot
        plt.savefig(save_dir + '/' + 'bandpass_' + str(round(center_wvl1[j], 2))+'_nm.png',
                    dpi=100)
        plt.close('all')





def calculate_in_band_irradiance(file_path_gaussian, file_path_interpol, radiance_file):
    """ This function performs the spectral convolution of radiance and the spectral bandpass.
    The spectral bandpass functions are saved in the directory in the form of csv file
    """
    save_dir = os.path.join(file_path_gaussian, 'Radiance_data')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    wavelength_solar = radiance_file['Wavelength']
    radiance = radiance_file['Response']#/max(radiance_file['Response'])
    band_pass_dir_gauss = r'spectral_bandpass_1400'
    band_pass_dir_interpol = r'look_up_table'

    band_pass_files_gaussian = ([each for each in os.listdir(os.path.join(file_path_gaussian,
                                                                          band_pass_dir_gauss))
                                 if each.endswith('_nm.csv')])
    band_pass_files_interpol = ([each for each in os.listdir(os.path.join(file_path_interpol,
                                                                          band_pass_dir_interpol))
                                 if each.endswith('_nm.csv')])


    cw_all = []
    dframe_radiance = pd.DataFrame()

    irradiance_all_gauss = []
    irradiance_all_interpol = []
    for bandpass_gauss, bandpass_interpol in zip(band_pass_files_gaussian,
                                                 band_pass_files_interpol):
        #print(bandpass_gauss, bandpass_interpol)


        # start with spectral integration with gaussian RSRs
        wavelength_bp_gauss = pd.read_csv(os.path.join(file_path_gaussian,
                                                       band_pass_dir_gauss,
                                                       bandpass_gauss))['Wavelength'].round(2)
        rsr_val_gauss = pd.read_csv(os.path.join(file_path_gaussian,
                                                 band_pass_dir_gauss,
                                                 bandpass_gauss))['Response']
        # interested only in solar wavelength within the bandpass
        wavelength_bp_interpol = pd.read_csv(os.path.join(file_path_interpol,
                                                          band_pass_dir_interpol,
                                                          bandpass_interpol))['Wavelength'].round(2)
        rsr_val_interpol = pd.read_csv(os.path.join(file_path_interpol,
                                                    band_pass_dir_interpol,
                                                    bandpass_interpol))['Radiance']
        mean_val_gauss = np.mean(wavelength_bp_gauss)
        sigma_val_gauss = np.std(wavelength_bp_gauss)
        best_vals_gauss, covarr = curve_fit(gauss_function, wavelength_bp_gauss, rsr_val_gauss,
                                            p0=[1, mean_val_gauss, sigma_val_gauss])
        print(covarr)
        mean_val_interpol = np.mean(wavelength_bp_interpol)
        sigma_val_interpol = np.std(wavelength_bp_interpol)
        best_vals_interpol, covarr = curve_fit(gauss_function,
                                               wavelength_bp_interpol, rsr_val_interpol,
                                               p0=[1, mean_val_interpol, sigma_val_interpol])
        offset = best_vals_interpol[1] - best_vals_gauss[1]


        wavelength_bp_gauss = wavelength_bp_gauss + offset
        cw_all.append(best_vals_interpol[1])




        solar_wvl_index = np.where((wavelength_solar >= min(wavelength_bp_gauss)) &
                                   (wavelength_solar <= max(wavelength_bp_gauss)))[0]
        solar_radiance = radiance[solar_wvl_index]
        solar_wvl = round(wavelength_solar[solar_wvl_index], 2)
        fit_params_rsr_gauss = interp1d(wavelength_bp_gauss, rsr_val_gauss,
                                        kind='slinear')
        fitted_rsr_gauss = fit_params_rsr_gauss(solar_wvl)
        numerator_gauss = np.sum(solar_radiance * fitted_rsr_gauss *solar_wvl)
        denominator_gauss = np.sum(fitted_rsr_gauss * solar_wvl)
        irradiance_gauss = numerator_gauss/denominator_gauss
        irradiance_all_gauss.append(irradiance_gauss)

        # Now the look up table RSR

        # interested only in solar wavelength within the bandpass
        solar_wvl_index = np.where((wavelength_solar >= min(wavelength_bp_interpol)) &
                                   (wavelength_solar <= max(wavelength_bp_interpol)))[0]
        solar_radiance = radiance[solar_wvl_index]
        solar_wvl = round(wavelength_solar[solar_wvl_index], 2)
        fit_params_rsr_interpol = interp1d(wavelength_bp_interpol, rsr_val_interpol, kind='slinear')
        fitted_rsr_interpol = fit_params_rsr_interpol(solar_wvl)
        numerator_interpol = np.sum(solar_radiance * fitted_rsr_interpol *solar_wvl)
        denominator_interpol = np.sum(fitted_rsr_interpol *solar_wvl)
        irradiance_interpol = numerator_interpol/denominator_interpol
        irradiance_all_interpol.append(irradiance_interpol)

        #plt.figure(figsize=(12, 5))
#        plt.plot(solar_wvl, solar_radiance/np.max(solar_radiance), 'k.--',
                  #markersize=3, linewidth=1.5)
#        #plt.plot(solar_wvl, irradiance_gauss, 'bo--')
#        #plt.plot(CW_interpol, irradiance_interpol,'m*--')
#        plt.plot(wavelength_bp_gauss, rsr_val_gauss, 'ro--', markersize=3)
#        plt.plot(wavelength_bp_interpol, rsr_val_interpol, 'go--', markersize=3)
#        plt.grid(True, linestyle=':')
#        plt.xlabel('Wavelength(nm)', fontsize=12)
#        plt.ylabel('Solar Irradiance (' + r'$Wm^{-2}nm^{-1}$' +')', fontsize=12)
#        plt.legend(['Solar Irradiance', 'Gaussian Bandpass','Look Up Table Bandpass'])
#        plt.title('Solar Irradiance Vs. Wavelength\n (SAO 2010 Solar Irradiance Spectrum)',
#                  fontsize=12)
#        #plt.pause(1.5)
#        plt.show()
#        plt.close('all')


    dframe_radiance['CW'] = np.array(cw_all)
    dframe_radiance['Gaussian'] = np.array(irradiance_all_gauss)
    dframe_radiance['Look_up_table'] = np.array(irradiance_all_interpol)
    irradiance_all_gauss = np.array(irradiance_all_gauss)
    irradiance_all_interpol = np.array(irradiance_all_interpol)

    dframe_radiance['difference'] = 100*(irradiance_all_gauss- irradiance_all_interpol)\
                                        /irradiance_all_gauss
    dframe_radiance.to_csv(save_dir+'/'+'radiance_simulation_gaussian.csv')


    plt.figure(figsize=(12, 5))
    #plt.plot(wavelength_solar, radiance, 'k', markersize=3, linewidth=1)
    #plt.show()
    plt.plot(np.array(cw_all), np.array(irradiance_all_gauss), 'r.--',
             markersize=10, linewidth=1)
    plt.plot(np.array(cw_all), np.array(irradiance_all_interpol), 'g.--',
             markersize=10, linewidth=1)

#    plt.plot(np.array(CW_all), 100*(irradiance_all_gauss- irradiance_all_interpol)\
#                                     /irradiance_all_gauss,'b*' )
    #plt.hist(np.array(irradiance_all_gauss)-np.array(irradiance_all_interpol),
              #bins=20, color='blue' )
    plt.grid(True, linestyle=':')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Solar Irradiance (' + r'$Wm^{-2}nm^{-1}$' +')', fontsize=12)
    #plt.legend(['Solar Irradiance', 'Gaussian Bandpass','Look Up Table Bandpass'])
    plt.text(600, -3, 'Spatial Pixel Index : 1400', fontsize=12,
             bbox=dict(facecolor='yellow', edgecolor='red'))
    #plt.legend(['Gaussian Bandpass','Look Up Table Bandpass'])
    plt.legend(['(Gaussian - Look Up Table)/Look Up Table'])

#    plt.suptitle('Spatial Pixel Index : 500)', fontsize=12,
#                 bbox=dict(facecolor='yellow', edgecolor='red'))
    plt.title('Solar Irradiance Vs. Wavelength\n (SAO 2010 Solar Irradiance Spectrum)',
              fontsize=12)
    #plt.ylim(-5, 5)
    plt.xlim(295, 740)
    plt.show()

    #return irradiance_all

#        print(CW)
#        print(irradiance)
#





def main():
    """ This is the main function. It calls all the other intermediate function
    """
#    pixel_to_wavelen_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\
#                              Spectral_Band_pass\Pixel_to_wavelen_map'

    file_path = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\
                  All_FWHM_only_Gaussian'
    radiance_file = read_radiance_data()
    file_path_2 = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\
                    All_FWHM\spectral_bandpass_1400'

    #start with Gaussian Bandpass
#    data_names = [each for each in os.listdir(file_path)
#                         if each.startswith("Params_Gauss")]
#
#
#    sample_data = []
#    for data_files in data_names[9:]:
#        #print(data_files)
#
#        wavelen_suffix = data_files.split('_')[-1]
#
#        pixel_to_wvl_map_data = sorted([each for each in os.listdir(pixel_to_wavelen_dir)
#                                        if each.endswith(wavelen_suffix)])
#
#        gaussian_files = os.path.join(file_path, data_files)
#
#        dframe = pd.read_csv(gaussian_files)
#        #dframe = dframe[['A1', 'A2', 'Sigma1', 'Sigma2']]
#        dframe = dframe[['A1', 'Sigma1']] # for Gaussian only
#        pixel_to_wav_map = os.path.join(pixel_to_wavelen_dir, pixel_to_wvl_map_data[0])
#        dframe1 = pd.read_csv(pixel_to_wav_map)
#        dframe['CW'] = dframe1['CW']
#        dframe = dframe.iloc[1400]
#        sample_data.append(dframe.values)
     # for flat top Gaussian
#    #gaussian_values = perform_spectral_interpolation(np.array(sample_data))

#    gaussian_values = perform_spectral_interpolation_only_gaussian(np.array(sample_data))
#
##
##    # Let us now create a spectral bandpass
#    #create_spectral_bandpass(gaussian_values, radiance_file, file_path) # flat top Gaussian
#    create_spectral_bandpass_only_gaussian(gaussian_values, radiance_file, file_path)
#
#
##    #Make sure that the center wavelength of Gaussians are the same
##    sample_val = []
##    data_names_interpol = sorted([each for each in os.listdir(file_path_2)
##                                 if each.endswith('csv')])
##    interpol_wavelen = []
##    interpol_rad = [ ]
##
##    for i in range(0, 64):
##        sub_sample_wvl = []
##        sub_sample_rad = []
##
##        for files in data_names_interpol[9:]:
##
##            interpol_rsr = os.path.join(file_path_2, files)
##            dframe = pd.read_csv(interpol_rsr, usecols=["wavelength", "rad"])
##
##            wavelength = dframe['wavelength'][i]
##            rad = dframe['rad'][i]
##            sub_sample_wvl.append(wavelength)
##            sub_sample_rad.append(rad)
##        dframe = perform_point_interpolation(sub_sample_wvl, sub_sample_rad,
                                               #np.array(sample_data)[:,-1])
##        interpol_rad.append(dframe['rad'].values)
##        interpol_wavelen.append(dframe['wavelength'].values)
##    create_spectral_bandpass_interpol(np.array(interpol_wavelen),
                                        #np.array(interpol_rad),
                                        #np.array(sample_data)[:,-1], file_path_2)
#    cc
##
#
##
###
##    # let us now perfrom spectral convolution with high res. radiance data
    calculate_in_band_irradiance(file_path, file_path_2, radiance_file)


if __name__ == "__main__":
    main()
