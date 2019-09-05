# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:10:42 2019

@author: nmishra
"""

import os
import h5py
import numpy as np

IMAGE_DIR = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\processed_h5\Mean_Processed_Data'
full_frame = np.genfromtxt(os.path.join(IMAGE_DIR,'Stitched_Radiance_VIS_V2.csv'), delimiter=',')
HF_FILE_NAME = os.path.join(IMAGE_DIR, 'Stitched_Radiance_VIS_V2.h5')
with h5py.File(HF_FILE_NAME, 'w') as hf:
    file_header = 'Processed_Data'
    hf.create_dataset(file_header, data=full_frame)
    

