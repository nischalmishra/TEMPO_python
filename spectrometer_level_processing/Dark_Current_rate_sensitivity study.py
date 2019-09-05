# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 07:58:14 2019

@author: nmishra
"""
import os
import numpy as np
import matplotlib.pyplot as plt
#from random import randint
import pandas as pd
from scipy.interpolate import interp1d



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

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def interpol_RMM(dframe2, wavelen):
   # print(wavelen)
    RMM_wvl = dframe2['Wavelength'].values
   # print(RMM_wvl)
    RMM_signal = dframe2['Total Signal'].values
    fit_params = interp1d(RMM_wvl, RMM_signal, kind='slinear')
    tempo_signal = fit_params(wavelen)
    #plt.plot (wavelen, tempo_signal,'ro')
    #plt.show()
    return tempo_signal/16.5
    
    

def main():
    """
    Main Function
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\processed_corrected_linearity\Dark_data'
    #save_dir = os.path.join(data_dir, 'Dark_current_sensitivity_study')

    dframe = pd.read_excel(r'C:\Users\nmishra\Desktop\pixel_to_wavelen.xlsx',
                           sheetname='Sheet1', index_col=None, header=None)
    dframe2 =  pd.read_excel(r'C:\Users\nmishra\Desktop\pixel_to_wavelen.xlsx', sheetname='RMM')    
    wavelen = dframe.values
    RMM_signal = interpol_RMM(dframe2, wavelen)
    
    dark_data_file = [each for each in os.listdir(data_dir) if each.endswith('.csv')]   
    result_array = np.empty([0, 90])
    unct_all = []
    mean_all = []
    dframe2 = pd.DataFrame()
    #numbers = [randint(0, 99) for p in range(0, 99)]
    #print(numbers)
    for i in range(len(dark_data_file), 0, -1):
        dark_data = np.genfromtxt(os.path.join(data_dir, dark_data_file[i-1]), delimiter=',')    
        #print(dark_data.shape)
        dark_data1 = dark_data[:, 1000]
        #dark_data1 = np.mean(dark_data, axis=1)
        #dark_data_unct = 100*(np.std(dark_data, axis=1)/np.sqrt(63))/dark_data1
        #dframe2 = pd.DataFrame(dark_data_unct)      
        #dframe2 = dframe2.rolling(window=20).median()
        #unct_est = dframe2.values
        #unct_est[unct_est>1.75]=np.median(unct_est)
       
       
        #dark_data_unct = running_mean(dark_data_unct, 5)
        #plt.plot
        
#        dc_rates = dark_data[:,1000]
#        plt.plot([abs(t - s )for s, t in zip(dc_rates, dc_rates[1:])], 'r', label='Spatial Pixel # 1000')
#        #plt.plot(dark_data[1500,:]/93, 'b', label='Wvl : 630.4 nm')
#                 
#        plt.title('TEMPO Dark Current Rates Vs. Wavelength', fontsize=14)
#        plt.ylabel('Dark Current Rates (DN/msec)', fontsize=14)
#        plt.xlabel('Wavelength (nm)', fontsize=14)
#        plt.legend(loc='best')
#        plt.grid(linestyle=':')
#        plt.show()
#        cc
        if i == 100:
            result_array = dark_data1

        else:
           #print(result_array.shape)
            result_array = np.dstack((result_array, dark_data1))          
            mean_array = np.mean(result_array, axis=2)
            std_array = np.std(result_array, axis=2)/np.sqrt(63)
            #unct = 100* (std_array/mean_array)
            unct_all.append(std_array)
            mean_all.append(mean_array/93)

    unct_all = np.array(unct_all)
    unct_all = np.squeeze(unct_all)
    mean_all = np.array(mean_all)
    mean_all = np.squeeze(mean_all)
    
    plt.figure(1)
    plt.plot(unct_all[:, 50], 'm.--', label='301.45 nm' )
    plt.plot(unct_all[:, 100], 'k.--', label='311.35 nm')
    plt.hold(True)
    plt.plot(unct_all[:, 400], 'b.--', label='370.70 nm')
    plt.plot(unct_all[:, 1000], 'g.--', label='489.33 nm')
    plt.plot(unct_all[:, 1500], 'r.--', label='630.4 nm')
    plt.title('Dark Current Std. Deviation Estimates Vs. Number of Images (Pixel #1000)', fontsize=12)
    plt.ylabel('1 Sigma Standard Deviation (DN)',  fontsize=12)
    plt.xlabel('Number of Dark Images from TVac Test')
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.show()
    
    
    plt.figure(2)
    plt.plot(100*unct_all[:, 50]/(RMM_signal[50]), 'm.--', label='301.45 nm' )
    plt.plot(100*unct_all[:, 100]/(RMM_signal[100]), 'k.--', label='311.35 nm')
    plt.hold(True)
    plt.plot(100*unct_all[:, 400]/(RMM_signal[400]), 'b.--', label='370.70 nm')
    plt.plot(100*unct_all[:, 1000]/(RMM_signal[1000]), 'g.--', label='489.33 nm')
    plt.plot(100*unct_all[:, 1500]/(RMM_signal[1500]), 'r.--', label='630.4 nm')
    plt.title('Dark Current Uncertainty Estimates Vs. Number of Images', fontsize=12)
    plt.ylabel('% Uncertainty Estimate (100* Std/Mean)',  fontsize=12)
    plt.xlabel('Number of Dark Images from TVac Test')
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.show()
    
    plt.figure(3)
    plt.plot(mean_all[:, 50], 'm.--', label='301.45 nm' )
    plt.plot(mean_all[:, 100], 'k.--', label='311.35 nm')
    plt.hold(True)
    plt.plot(mean_all[:, 400], 'b.--', label='370.70 nm')
    plt.plot(mean_all[:, 1000], 'g.--', label='489.33 nm')
    plt.plot(mean_all[:, 1500], 'r.--', label='630.4 nm')
    plt.title('Dark Current Rates Vs. Number of Images (Pixel # 1000)', fontsize=12)
    plt.ylabel('Dark Current Rates (DN/ms)',  fontsize=12)
    plt.xlabel('Number of Dark Images from TVac Test')
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.show()
    
    
    #read the pixel to wavelen map sample
    plt.figure(4)
    #print(len(RMM_signal))
    #print(len(unct_all[2,:]))
    #temporal_unct = unct_all/RMM_signal.T
    #print(temporal_unct.shape)    
    temporal_unct = 100*unct_all[40, :]/RMM_signal.T    
    plt.plot(np.squeeze(wavelen.T), np.squeeze(temporal_unct), 'r.--', label= 'Pixel # 1000, 20 images, 63 coadds')
    plt.title('% Radiometric Uncertainty Due To Dark Current Correction  Vs. Wavelength', fontsize=12)
    plt.ylabel('% Radiometric Uncertainty', fontsize=12)
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.legend(loc='best')
    plt.grid(linestyle=':')

    plt.show()
    #cc




if __name__ == "__main__":
    main()
