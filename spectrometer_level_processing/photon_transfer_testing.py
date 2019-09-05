# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 08:24:53 2018

@author: nmishra
"""

import os
import numpy as np
import pandas as pd

image_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\LSA Data 27-28Jul17\07282017'
save_dir = r'C:\Users\nmishra\Desktop\Photon_transfer'


data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.tdms')])

print(*data_path_all, sep= "\n")
cc

actual_files = pd.read_csv (os.path.join(save_dir,'Photon_Transfer.csv'))
log_files_list = actual_files['Scenes'].values

Difference = set(data_path_all).difference(log_files_list)
print(*Difference, sep= "\n")

