# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:30:40 2018

@author: nmishra
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from nptdms import TdmsFile


TDMS_DATA_PATH = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\LSA Data 27-28Jul17\07272017'
PLOT_SAVE_DIR = os.path.join(TDMS_DATA_PATH, 'saved_data')
if not os.path.exists(PLOT_SAVE_DIR):
    os.makedirs(PLOT_SAVE_DIR)
FILE_NAME =[each for each in os.listdir(TDMS_DATA_PATH)
             if each.endswith('07272017_174636_WL675.tdms')]
print(len(FILE_NAME))
#cc

k=0

for files in range(0, len(FILE_NAME)):   
    full_file_path = os.path.join(TDMS_DATA_PATH, FILE_NAME[files])   
    #print(os.path.getsize(full_file_path))
    if os.path.getsize(full_file_path) > 10:
        amp = []
        wavelen = []
       
        
        lsa_file_name = FILE_NAME[files].split('.tdms')[0]
        wavelength_val = lsa_file_name.split('WL')[1]
        wavelength_val = wavelength_val.replace('_', '.')      
        print(full_file_path)
        tdms_file = TdmsFile(full_file_path)
        tdms_groups = tdms_file.groups()       
        plt.figure(figsize=(7, 5))
        dframe = pd.DataFrame()        
        for i in range(0, len(tdms_groups)-1):
            MessageData_channel_1 = tdms_file.object(str(tdms_groups[i]), 'Wavelength')
            MessageData_channel_2 = tdms_file.object(str(tdms_groups[i]), 'Amplitude')
            wavelength = MessageData_channel_1.data                        
            amplitude = MessageData_channel_2.data
            amplitude = amplitude#/np.max(amplitude)
            max_amp = np.where(amplitude == max(amplitude))[0]
            max_amp = max_amp[0]
            max_wvl = wavelength[max_amp]
            PEAK_WVL = 675
            if max_wvl> PEAK_WVL:
                wavelength = wavelength - (max_wvl-PEAK_WVL)
            else:
                wavelength = wavelength - (max_wvl-PEAK_WVL)
                
            wavelen.append(wavelength)            
            dframe['wavelength'] = wavelength
            #dframe['Amp'] = amplitude
            #dframe.to_csv(r'C:\Users\nmishra\Desktop\LSA_check.csv', sep=’,’)
            #cc
            
            plt.plot(wavelength, amplitude/np.max(amplitude), label=str(tdms_groups[i]))
            #plt.show()
            #cc            
            #plt.hold(True)
            plt.xlim(float(wavelength_val)-0.25, float(wavelength_val)+0.25)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Peak Normalized Response')
            #plt.yscale('log')
            plt.grid(True, linestyle=':')
            amp.append(amplitude/max(amplitude))
            
        #plt.xlim(max_wvl-0.18, max_wvl+0.18)
        #plt.ylim(-0.05, 0.75)
        #plt.legend(['*Different time stamps'],scatterpoints=0,loc='best', fancybox=True,
                   #framealpha=1, shadow=True, borderpad=1)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr = '** Multiple timestamps'
        #plt.ylim(np.mean(wavelength)-0.01, np.mean(wavelength+0.01))
        #plt.text(max_wvl-0.17, 0.02, textstr, bbox=props)
        plt.title('Normalized Response Vs. Wavelength \n Laser Spectrum Analyzer' +
                  ' (WL = ' +str(wavelength_val) + ' nm)')
        #fig_name = PLOT_SAVE_DIR  +'/' + 'LSA_'+wavelength_val+'_nm'+str(k)+'.png'
        #plt.savefig(fig_name, bbox_inches="tight")
        #plt.show() 
        #cc
        #plt.close('all')
        k= k+1

        nx, ny = np.array(wavelen).shape
        wavelen = np.reshape(np.array(wavelen), [nx*ny, 1])
            
#        nx, ny = np.array(amp).shape
        amp = np.reshape(np.array(amp), [nx*ny, 1]) 
        
        dframe1 = pd.DataFrame()
        dframe2 = pd.DataFrame()
        dframe = pd.DataFrame()
        
        dframe['wavelen'] = wavelen.flatten()        
        dframe['Amp'] = amp.flatten()
        
        dframe['wavelen'] = np.round(dframe['wavelen'], 2)
        dframe['Amp'] = np.round(dframe['Amp'], 2)    
        dframe1 = dframe.groupby(['wavelen']).mean()
        dframe1 = dframe1/dframe1.max()    
        ax = dframe1.plot(style='.',label='Mean', color='red') 
        print(len(dframe1))
        plt.legend(['Mean'],loc='best')        
        plt.show()
        #cc
        csv_name = PLOT_SAVE_DIR  +'/' + 'LSA_'+str(PEAK_WVL)+'_nm.csv'
        dframe1.to_csv(csv_name, sep=',', header=True, index=True)
