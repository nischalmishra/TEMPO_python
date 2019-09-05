# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:12:41 2018

@author: nmishra
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.interpolate import interp1d

def main():
    # Lets read the bandpass data
    bandpass_data = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\Pixel_to_wavelen_map'
    spatial_range = np.arange(4, 2041) # only 2037 spatial pixels are illuminated   
    # Lets create a big dataframe
    
    dframe_297_8nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_297.8nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_310nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_310.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_320nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_320.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_330nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_330.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_355nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_355.0nm.csv' ), usecols =['CW'], delimiter=',') 
    dframe_390nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_390.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_425nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_425.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_460nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_460.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_488_2nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_488.2nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_541_8nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_541.8nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_605nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_605.0nm.csv' ), usecols =['CW'], delimiter=',') 
    dframe_640nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_640.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_675nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_675.0nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_727_5nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_727.5nm.csv' ), usecols =['CW'], delimiter=',')
    dframe_736nm = pd.read_csv(os.path.join(bandpass_data,'Gaussian_fit_params_736.0nm.csv' ), usecols =['CW'], delimiter=',')
    
   
    dframe_UV = pd.concat([dframe_297_8nm, dframe_310nm,dframe_320nm, dframe_330nm, 
                           dframe_355nm, dframe_390nm, dframe_425nm, dframe_460nm,
                           dframe_488_2nm], axis=1)
    dframe_VIS = pd.concat([dframe_541_8nm, dframe_605nm, dframe_640nm,
                                dframe_675nm, dframe_727_5nm, dframe_736nm],
                                axis=1)
    spectral_indices_UV = [29, 91, 142, 192, 319, 496, 673, 849, 992]
    resampled_pixel_UV= np.arange(np.min(spectral_indices_UV), np.max(spectral_indices_UV), 1)
    spectral_indices_VIS = [1050, 1369, 1546, 1724, 1989, 2032]
    resampled_pixel_VIS= np.arange(np.min(spectral_indices_VIS), np.max(spectral_indices_VIS), 1)   
    wavelen_resamp_all_UV = []
    wavelen_resamp_all_VIS = []
    spatial_len = len(dframe_297_8nm['CW'])
    
    
    for i in np.arange(0, spatial_len):
        wavelen_data = dframe_UV['CW'].loc[i].values
        fit_params = interp1d(spectral_indices_UV, wavelen_data, kind='slinear')  
        wavelen_data_resamp = fit_params(resampled_pixel_UV)
        wavelen_resamp_all_UV.append(wavelen_data_resamp)
   
    file_name_wavelen = bandpass_data+'/'+'pixel_to_wavelen_map_UV.csv'
    file_name_spectral_indices = bandpass_data+'/'+'spectral_indices_UV.csv' 
    np.savetxt(file_name_wavelen, wavelen_resamp_all_UV, delimiter=",", fmt='%1.3f')
    np.savetxt(file_name_spectral_indices , resampled_pixel_UV, delimiter=",")
     
    for i in np.arange(0, spatial_len):
        wavelen_data = dframe_VIS['CW'].loc[i].values
        fit_params = interp1d(spectral_indices_VIS, wavelen_data, kind='slinear')  
        wavelen_data_resamp = fit_params(resampled_pixel_VIS)
        wavelen_resamp_all_VIS.append(wavelen_data_resamp)
        
    file_name_wavelen = bandpass_data+'/'+'pixel_to_wavelen_map_VIS.csv'
    file_name_spectral_indices = bandpass_data+'/'+'spectral_indices_VIS.csv'        
    np.savetxt(file_name_wavelen, wavelen_resamp_all_VIS, delimiter=",", fmt='%1.3f')
    np.savetxt(file_name_spectral_indices , resampled_pixel_VIS, delimiter=",")
         
   
    
    
    
    
        
                     
       
        
        
    
    
    
    
        
        
                          
    
    
if __name__ == "__main__":
    main()