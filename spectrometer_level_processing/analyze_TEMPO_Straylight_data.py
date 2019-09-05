# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:49:55 2018

@author: nmishra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
DATA_DIR = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC\saved_quads\saved_plots_straylight'

INT_FILE = np.loadtxt(os.path.join(DATA_DIR, 'int_time.csv'), delimiter=',')

PIXEL_TO_WAVELEN_MAP = np.loadtxt(os.path.join(DATA_DIR, 'pixel_to_wavelen_map.csv'),
                                  delimiter=',')
WAVELEN_ALL = PIXEL_TO_WAVELEN_MAP[:, 500]
TEMPO_DATA_ALL = sorted([each for each in os.listdir(DATA_DIR)
                         if each.endswith('.sav.csv')])


REF_WAVELEN = [10, 50, 100, 300, 700, 1000, 1500, 1750, 2000]
COLOR = ['red', 'blue', 'green', 'black', 'magenta', 'orange', 'purple',
         'yellow', 'cyan']
COUNT = 0
for i in REF_WAVELEN:
    wavelen_val = WAVELEN_ALL[i]
    wavelen_title = str(round(wavelen_val, 2)) + ' nm'
    all_straylight = []
    dframe = pd.DataFrame()
    for data_path in TEMPO_DATA_ALL:
        straylight_data = np.loadtxt(os.path.join(DATA_DIR, data_path), delimiter=',')
        straylight_data = straylight_data[i]
        all_straylight.append(straylight_data)

    dframe['Int_time'] = INT_FILE
    dframe['Straylight'] = all_straylight
    dframe = dframe.groupby(['Int_time']).mean()

    plt.plot(dframe, 'o', color=COLOR[COUNT], label=wavelen_title)
    plt.grid(True, linestyle=':')
    plt.xlabel('Integration Time (ms)', fontsize=12)
    plt.ylabel('StrayLight Estimates (DN)', fontsize=12)
    plt.title('StrayLight Estimates Vs Integration Time', fontsize=12)
    plt.legend(loc='best')
    save_dir = os.path.join(DATA_DIR, 'few_saved_data')
    file_name = save_dir+'/'+ 'straylight_'+ wavelen_title+'.png'
    plt.savefig(file_name, dpi=100, bbox_inches='tight')
    COUNT = COUNT+1
    save_dir = os.path.join(DATA_DIR, 'few_saved_data')
    file_name = save_dir+'/'+ 'straylight_'+ wavelen_title+'.csv'
    dframe.to_csv(file_name, sep = ',', float_format ='%1.2f', header=True)

plt.show()



