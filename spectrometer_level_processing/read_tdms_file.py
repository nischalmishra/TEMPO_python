# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:30:40 2018

@author: nmishra
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from nptdms import TdmsFile
import os
from matplotlib.pyplot import cm

tdms_data_path = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\LSA Data 27-28Jul17\07272017'
file_name = r'07272017_132129_WL450.tdms'

full_file_path = os.path.join(tdms_data_path, file_name)


tdms_file = TdmsFile(full_file_path)
tdms_groups = tdms_file.groups()
amp = []
wavelen = []
wavelength_val = full_file_path.split('.')[0].split('_')[-1][2:]
color=iter(cm.rainbow(np.linspace(0,1,len(tdms_groups))))
dframe = pd.DataFrame()
dframe1 = pd.DataFrame()
dframe2 = pd.DataFrame()

for i in range(0, len(tdms_groups)):
    MessageData_channel_1 = tdms_file.object(str(tdms_groups[i]),'Wavelength')
    MessageData_channel_2 = tdms_file.object(str(tdms_groups[i]),'Amplitude')    
    wavelength = MessageData_channel_1.data    
    wavelen.append(wavelength)
    #print(wavelen)
    amplitude = MessageData_channel_2.data
    amplitude = amplitude#/np.max(amplitude)    
    max_amp = np.where(amplitude == max(amplitude))[0]    
    max_amp= max_amp[0]
    max_wvl = wavelength[max_amp] 
#    plt.figure(num=i, figsize=(17,8))
#    plt.plot(wavelength, amplitude,'r.--', label = 'Peak Wavelength = ' +str(np.round(max_wvl, 3))+ ' nm')
#    #plt.xlim(673, 677)
#    plt.xlabel('Wavelength (nm)')
#    plt.ylabel('Normalized Amplitude')
#    plt.yscale('log')
#    plt.title('Amplitude Vs. Wavelength data from Laser Spectrum Analyzer\n' + 'Wavelength = ' +str(wavelength_val) + ' nm'+ ', Time = ' + str(tdms_groups[i]) )    
#    plt.grid(True, linestyle =':')
#    plt.legend(loc='lower left')
#    plt.pause(2)
#    plt.show(block=False)
#    plt.close('all')
    amp.append(amplitude)
#plt.show()
nx, ny = np.array(wavelen).shape
print(nx, ny)
cc     
#wavelen = np.reshape(np.array(wavelen), [nx*ny, 1])
#
#nx = None
#ny= None
#
#nx, ny = np.array(amp).shape
#amp = np.reshape(np.array(amp), [1, nx*ny])
#
#
#dframe['wavelen'] = wavelen.flatten()
#
#dframe['Amp'] = amp.flatten()
##print(dframe.shape)
#dframe['wavelen'] = np.round(dframe['wavelen'], 3)
#dframe['Amp'] = np.round(dframe['Amp'], 3)
#
#
#
#dframe2 = 100*dframe.groupby(['wavelen']).std()/dframe.groupby(['wavelen']).mean()
##dframe2 = dframe2/dframe2.max()
#dframe1 = dframe.groupby(['wavelen']).mean()
#dframe1 = dframe1/dframe1.max()
#
#ax = dframe2.plot(label='Sum', color='red', marker='o')
##dframe1.plot(label='Mean', color='blue', ax=ax, marker='*')
##plt.legend(['Sum','Mean'],loc='best')
#
#plt.show()
#cc
#dframe.reset_index()
#
#
#plt.plot(dframe['wavelen'], dframe['Amp'],'b.')
#plt.show()
#cc

plt.plot(wavelen, amp,'.', markersize=0.5)
plt.grid(True, linestyle =':')
#plt.xlim(674, 676)
plt.title('Amplitude Vs. Wavelength plot from Spectrum Analyzer (' + str(wavelength_val) + ' nm)')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Amplitude')
plt.show()
cc
#
##plt.show(block=False)
#wavelen = np.array(wavelen)
#amp = np.array(amp)
#rows, cols = wavelen.shape
#wavelen = np.reshape(wavelen, (rows* cols, 1))
#amp = np.reshape(amp, (rows* cols, 1))
#final_data = np.concatenate((wavelen, amp), axis=1)
#save_bandpass = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.07.30\saved_quads\675_nm\saved_quads\processed_image' + '/'+'spec_analyzer_normalized.csv'
#np.savetxt(save_bandpass, final_data, delimiter=",")
#df = pd.read_csv(r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.07.30\saved_quads\675_nm\saved_quads\processed_image\spec_analyzer_normalized.csv')
#df.columns = ['wvl', '675 nm']
#df['wvl'] = np.round(df['wvl'], 3)
#df['675 nm'] = np.around(df['675 nm'], 6)
##print(df)
#df1 = df
##df.groupby('wvl').mean()
##df.groupby('wvl').mean().plot(color='red', marker='.', linewidth=3, label='675 nm')
#plt.plot(df['wvl'], df['675 nm'], 'b.', markersize=0.1, label='before smoothing')
#plt.xlim(674, 676)
#plt.xlabel('Wavelength')
#plt.ylabel('Amplitude')
#plt.title('Amplitude Vs. Wavelength plot from Spectrum Analyzer')
#plt.grid(True, linestyle=':')
#plt.legend(loc='best')
#plt.show(block=False)


# Now read the intensity file
data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.07.30\saved_quads\675_nm\saved_quads\processed_image'
file_name ='bandpass_radiance.csv'
bandpass_data = np.loadtxt(os.path.join(data_dir, file_name), delimiter = ",")
window  = 15
test_data = bandpass_data[:, 100]
weights = np.repeat(1.0, window)/window
sma = np.convolve(test_data, weights, 'valid')
plt.figure(2)
plt.plot(test_data,'b.', markersize=3, label='before smoothing')
#plt.plot(bandpass_data[:, 100], 'r.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 500], 'b.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 900], 'g.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 1300], 'm.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 1500], 'y.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 1600], 'k.', markersize=1, label='before smoothing')
#plt.plot(bandpass_data[:, 2000], 'orange',marker='.', markersize=1, label='before smoothing')
#plt.plot(sma, 'r.--', linewidth=3, label='675 nm (pixel =700)')
plt.legend(loc='best')
plt.grid(True, linestyle=':')
plt.ylabel('Counts (DN)')
plt.xlabel('# Number of observations (Temporal sequence)')
plt.title('Counts Vs. Number of observations')

plt.show(block=False)
plt.figure(3)
plt.plot(bandpass_data[20, :]/np.max(bandpass_data[20, :]), 'm')
plt.plot(bandpass_data[35, :]/np.max(bandpass_data[35, :]), 'lime')
plt.plot(bandpass_data[50, :]/np.max(bandpass_data[50, :]), 'r')
plt.plot(bandpass_data[44, :]/np.max(bandpass_data[44, :]), 'navy')
plt.plot(bandpass_data[70, :]/np.max(bandpass_data[70, :]), 'b')
plt.plot(bandpass_data[73,:] /np.max(bandpass_data[73, :]),  'cyan')
plt.plot(bandpass_data[80, :]/np.max(bandpass_data[80, :]), 'k')
plt.plot(bandpass_data[100, :]/np.max(bandpass_data[100, :]), 'purple')
plt.plot(bandpass_data[110, :]/np.max(bandpass_data[110, :]), 'orange')
plt.plot(bandpass_data[120, :]/np.max(bandpass_data[120, :]), 'g')
plt.grid(True, linestyle=':')
plt.ylabel('Counts (DN)')
plt.xlabel('Spatial Pixel Index (#)')
plt.title('Spatial profile @ different illuminations for spectral bandpass 675 nm\n (Spectral Index = 1714)')
plt.show()





