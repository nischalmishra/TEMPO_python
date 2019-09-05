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
#import scipy.signal
import pandas as pd
from scipy.signal import convolve as conv2
from skimage import color, data, restoration

#pylint: disable= E1101
#pylint: disable-msg=R0912
#pylint: disable-msg=R0914
#pylint: disable-msg=R0915


def gaussian(wavelength, amp, center, width):
    """
    Fit the  Gaussian FUnction
    """

    gauss_fxn = amp * np.exp(-(wavelength-center)**2/ 2*width**2)
    return gauss_fxn



def peak_allign_spectral_bp_LSA_data(bandpass_data, bandpass_wvl, LSA_data):
    """
    Function to allign peak of LSA data to peak of bandpass_data
    """

    LSA_wvl_all = []
    LSA_wvl = LSA_data[:, 0]
    LSA_data = LSA_data[:, 1]
    init_vals = [1, np.mean(LSA_wvl), np.std(LSA_wvl)]
    best_vals2, covarr2 = curve_fit(gaussian, LSA_wvl,
                                    LSA_data, p0=[init_vals])
    for i in range(0, np.shape(bandpass_wvl)[1]):

        normalized_data = bandpass_data[:, i]/np.max(bandpass_data[:, i])
        wavelength = bandpass_wvl[:, i]
        #plt.plot(wavelength, normalized_data, 'ro--', label='Bandpass Data')
        #plt.plot(LSA_wvl, LSA_data, color='blue', label='LSA Data before peak allignment')
        width = np.std(wavelength)
        center = np.mean(wavelength)
        amp = 1
        init_vals = [amp, center, width]
        best_vals1, covarr1 = curve_fit(gaussian, wavelength,
                                        normalized_data, p0=[init_vals])

        if best_vals1[1] > best_vals2[2]:
            LSA_wvl_new = LSA_wvl + (best_vals1[1]- best_vals2[1])
        else:
            LSA_wvl_new = LSA_wvl - (best_vals1[1]- best_vals2[1])
        LSA_wvl_all.append(LSA_wvl_new)
        #LSA_wvl=[]
        #LSA_data=[]

#    plt.plot(LSA_wvl, LSA_data, color='green', label='LSA Data after peak allignment')
#    plt.grid(True, linestyle=':')
#    plt.legend(loc='best')
#    plt.xlabel('Wavelength (nm)')
#    plt.ylabel('Normalized Response')
#    plt.title('Spectral Bandpass and LSA Data (Along 425 nm, pixel # 500)')
#    plt.xlim(best_vals1[1]-0.9, best_vals1[1]+0.9)
#    plt.show()
#    cc

    return np.array(LSA_wvl_all).T
def perform_richardson_lucy_deconv(bandpass_data, bandpass_wvl, LSA_wvl,
                            LSA_data, save_dir, wavelen_val):
    """ Function to perform the spectral deconvolution of Laser Linewidth from
    the spectral bandpass data
    """
    #print(LSA_wvl.shape)
    #print(LSA_data.shape)
   # print(bandpass_data.shape)
    for i in range(1000, 1001):#np.shape(bandpass_data)[1]):
        normalized_data = bandpass_data[:, i]/np.max(bandpass_data[:, i])
        wavelength = bandpass_wvl[:, i]       
        bp_data = np.array([wavelength, normalized_data]).T
        bp_data = bp_data[bp_data[:, 0].argsort()]
        wavelength = bp_data[:, 0]
        normalized_data = bp_data[:, 1]
       
    
        #LSA_data = LSA_data[max_loc-480 :max_loc+480]
        if len(LSA_data)>len(normalized_data):  
            factor = len(LSA_data)/len(normalized_data)
            wvl_LSA = LSA_wvl[:, i]           
            LSA_data = LSA_data[::int(factor)+1]           
            wvl_LSA = wvl_LSA[::int(factor)+1]
            #LSA_data = LSA_data-0.010
            LSA_data = LSA_data/np.max(LSA_data)           
            CW = np.where(LSA_data==max(LSA_data))[0]
            val = int(25)
            wvl_LSA = wvl_LSA[int(CW)-val : int(CW)+val]
            LSA_data = LSA_data[int(CW)-val : int(CW)+val]-0.01
            LSA_data = LSA_data/np.max(LSA_data)
            
           
            
          
           
        
        
#        plt.figure(); plt.plot(wvl_LSA, LSA_data,'r')
#        plt.plot(wavelength, normalized_data,'b.')
#        plt.grid(True, linestyle=':')
#        plt.show()
        
        
        normalized_data1 = conv2(normalized_data , LSA_data, mode='same')
        normalized_data1 = normalized_data1/np.max(normalized_data1)        
        deconvolved_RL = restoration.richardson_lucy(np.array([normalized_data]), np.array([LSA_data]), iterations=2)
        deconvolved_RL = deconvolved_RL.T
        deconvolved_RL =  deconvolved_RL.real
        print(deconvolved_RL/np.max(deconvolved_RL))
        
#        cc
        plt.figure()
        
        
        np.savetxt(r"C:\Users\nmishra\Desktop\spectral_deconv\foo.csv", deconvolved_RL, delimiter=",")
        plt.plot(wavelength, normalized_data,'b')
        plt.plot(wavelength, deconvolved_RL, 'r')
        plt.plot(wvl_LSA, LSA_data,'go--')        
        plt.show()
        cc

       
       
       
      
   
      
def perform_spectral_deconv(bandpass_data, bandpass_wvl, LSA_wvl,
                            LSA_data, save_dir, wavelen_val):
    """ Function to perform the spectral deconvolution of Laser Linewidth from
    the spectral bandpass data
    """
    #print(LSA_wvl.shape)
    #print(LSA_data.shape)
   # print(bandpass_data.shape)
    for i in range(1000, 1001):#np.shape(bandpass_data)[1]):
        normalized_data = bandpass_data[:, i]/np.max(bandpass_data[:, i])
        wavelength = bandpass_wvl[:, i]       
        bp_data = np.array([wavelength, normalized_data]).T
        bp_data = bp_data[bp_data[:, 0].argsort()]
        wavelength = bp_data[:, 0]
        normalized_data = bp_data[:, 1]
       
    
        #LSA_data = LSA_data[max_loc-480 :max_loc+480]
        if len(LSA_data)>len(normalized_data):  
            factor = len(LSA_data)/len(normalized_data)
            wvl_LSA = LSA_wvl[:, i]           
            LSA_data = LSA_data[::int(factor)+1]           
            wvl_LSA = wvl_LSA[::int(factor)+1]
            LSA_data = LSA_data-0.010
            LSA_data = LSA_data/np.max(LSA_data)
            CW = np.where(LSA_data==max(LSA_data))[0]
            val = int(25)
            wvl_LSA = wvl_LSA[int(CW)-val : int(CW)+val]
            LSA_data = LSA_data[int(CW)-val : int(CW)+val]-0.01
            LSA_data = LSA_data/np.max(LSA_data)
            #fit_params = interp1d(wvl_LSA,LSA_data)
            #resamp_data = fit_params(wavelength)
            #LSA_data = resamp_data/np.max(resamp_data)
           
        
        
#        plt.figure(); plt.plot(wavelength,LSA_data,'r')
#        plt.plot(wavelength, normalized_data,'b.')
#        plt.show()
#        cc
##       
#          
       
        yfreq = np.fft.fft(normalized_data)        
        max_amp = yfreq.max()
        winfreq = np.fft.fft(LSA_data)
#        plt.plot((winfreq),'b.')
#        plt.plot(yfreq,'r.')
#        plt.show()
#        cc
        #The problem is that the LSA has near-zero amplitude at high
        #frequencies, so you're blowing up the high-frequency content of the noisy
        #signal when you divide in the frequency domain.
        
#        Basically, you just:
#            1) convert to the frequency domain
#           2) replace any amplitudes below some threshold with that threshold in the
#            signal you're dividing by (the window, in your case)
#            3) pad the lengths to match
#            4) divide (the actual deconvolution)
#            5) convert back to the time domain
        eps = 1.0e-7
        winfreq[winfreq < eps] = eps
        padded =eps *np.ones_like(yfreq)
        padded[:winfreq.size] = winfreq
        newfreq = yfreq / padded
        newfreq *= max_amp / newfreq.max()
        deconv1 = np.fft.ifft(newfreq)        
        deconv2= deconv1.real
        deconv2 = deconv2+0.11
        deconv2 = deconv2/np.max(deconv2)
        
        #deconv2[deconv2 <= 0.0] = 0.0  
        #deconv2 = deconv2/np.max(deconv2)
        #error = normalized_data-deconv2
        fwhm = full_width_half_max(wavelength, normalized_data) 
        fwhm1 = full_width_half_max(wavelength, deconv2)
        #fwhm=0
        #fwhm1=0       
        print(fwhm1, fwhm)
        
        
#        plt.figure()
#        plt.plot(wavelength, normalized_data, 'r.:', label='Original (' + 'FWHM = ' +str(round(fwhm, 3 )) + ' nm)')
#        plt.plot(wavelength, deconv1, 'b--*', label='Deconvolved (' + 'FWHM = ' +str(round(fwhm1, 3 )) + ' nm)')
#        plt.plot(wvl_LSA , LSA_data, 'k.--', label='LSA Spectra')
#        plt.grid(True, linestyle=':')
#        plt.xlabel('Wavelength (nm)')
#        plt.ylabel('Normalized Counts')
#        plt.title('TEMPO Spectral Bandpass (Along '+ str(wavelen_val) +' nm, pixel# ' +str(i)+')')
#        plt.xlim(float(wavelen_val)-1.5, float(wavelen_val)+1.5)        
#        plt.legend(loc='best')
#        plt.show()
        
        
       
        plt.figure()
        plt.plot(wavelength, normalized_data, 'r.:', label='Original (' + 'FWHM = ' +str(round(fwhm, 3 )) + ' nm)')
        plt.plot(wavelength, (deconv2)/np.max(deconv2), 'b.--', label='Deconvolved (' + 'FWHM = ' +str(round(fwhm1, 3 )) + ' nm)')
        plt.plot(wvl_LSA, LSA_data, 'k', label='LSA Spectra')
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Normalized Counts')
        plt.xlim(float(wavelen_val)-1.5, float(wavelen_val)+1.5)  
        plt.title('TEMPO Spectral Bandpass (Along '+ str(wavelen_val) +' nm, pixel# ' +str(i)+')')
        plt.legend(loc='best')
        plt.show()
        cc
        
        
        # SANITY CHECK
        if len(normalized_data)>len(LSA_data):
            diff_len = len(normalized_data)- len(LSA_data)
            zero_pads = np.zeros(diff_len)
            LSA_data = np.concatenate((LSA_data, zero_pads))
        data = np.fft.fft(normalized_data)
        noise = np.fft.fft((LSA_data))
        muxed_signal = np.multiply(data, noise)
        convolved_signal = np.fft.fftshift(np.real(np.fft.ifft(muxed_signal)))        
        
        
        yfreq1 = np.fft.fft(convolved_signal)       
        #max_amp1 = yfreq1.max()
        winfreq1 = np.fft.fft(LSA_data)
        #winfreq1[winfreq1 < 0.001] = 0.001
        #padded1 = 0.001*np.ones_like(yfreq1)
        #padded1[:winfreq1.size] = winfreq1
        newfreq1 = yfreq1/winfreq1
        #newfreq1 *= max_amp1 / newfreq1.max()
        deconv1 = np.fft.fftshift(np.real(np.fft.ifft(newfreq1)))        
        deconv2= deconv1.real
        deconv2 = deconv2/np.max(deconv2)
        #deconv2[deconv2 <= 0.0] = 0.0  
       
        
        
        plt.figure(); #plt.plot(wavelength, convolved_signal/np.max(convolved_signal),'b.')
        plt.plot(wavelength, deconv2/np.max(deconv2),'r.--')
        plt.plot(wavelength, normalized_data,'b.--')
        plt.show()
        cc
   
    return





def flat_top_gaussian(wavelength, amp, center, width):
    """
    Fit the flat top Gaussian FUnction   """

    gauss_fxn = amp * np.exp(-(wavelength-center)**4/ 2*width**4)
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
        plt.plot(wavelength, final_fit-normalized_data, 'k.', label='fit- data')
        plt.ylim(-0.06, 0.06)
        plt.grid(True, linestyle=':')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Residue')
        plt.title('Difference Between Fit and Actual Data (Along 425 nm, pixel # 1500)')
        plt.legend(loc='best')
        plt.show()

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
    plt.ioff()
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Spectral_Band_Pass_Data'
    channel = 'VIS'
    channel_data = os.path.join(data_dir, channel)
    all_wavelengths_file = sorted([each for each in os.listdir(channel_data)
                                   if each.endswith('_nm')])
    #del all_wavelengths_file[4]  
    
    
    for files in range(2, len(all_wavelengths_file)):
        save_dir = os.path.join(channel_data, all_wavelengths_file[files])
        wavelen_val = all_wavelengths_file[files].strip('_nm')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        bp_data_path = os.path.join(channel_data, all_wavelengths_file[files],
                                    r'bandpass_lookup_table')
        bandpass_data = np.genfromtxt(os.path.join(bp_data_path,
                                                   'bandpass_radiance_normalized_V1.csv'),
                                      delimiter=',')
        bandpass_wvl = np.genfromtxt(os.path.join(bp_data_path, 'bandpass_wvl_V1.csv'),
                                     delimiter=',')
        LSA_data_name = [each for each in os.listdir(save_dir)
                         if each.startswith('LSA')]
        LSA_data = np.genfromtxt(os.path.join(save_dir, LSA_data_name[0]),
                                 delimiter=',', skip_header=1)
        LSA_wvl = peak_allign_spectral_bp_LSA_data(bandpass_data, bandpass_wvl, LSA_data)
        perform_spectral_deconv(bandpass_data, bandpass_wvl, LSA_wvl,
                                LSA_data[:, 1], save_dir, wavelen_val)
        #fit_double_gaussian_func(bandpass_data, bandpass_wvl, save_dir, wavelen_val)

if __name__ == "__main__":
    main()
