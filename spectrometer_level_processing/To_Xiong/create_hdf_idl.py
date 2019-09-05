# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:54:10 2019

code to save the idl file as hdf file

@author: nmishra
"""


import os
#import numpy as np
import h5py
from scipy.io.idl import readsav

# pylint: disable=E1101

def read_idl_file(data):
    """Read IDL variables. Dave's IDL code to read spectrometer data runs
     much faster than my Python script. Also if there are .sav files,
     I can use them in  future rather than reading the whole file again
     """
    read_variable = readsav(data)
    try:
        full_frame_image = read_variable.quads
    except KeyError:
        full_frame_image = read_variable.q
    return full_frame_image

IMAGE_DIR = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads'
OUTPUT_DIR = os.path.join(IMAGE_DIR, 'saved_hdf_input_files')
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
DATA_PATH_ALL = sorted([each for each in os.listdir(IMAGE_DIR)
                        if each.endswith('.sav')])
QUADS = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
for data_path in DATA_PATH_ALL:
    HF_FILE_NAME = os.path.join(OUTPUT_DIR, data_path.strip('saved_quads').strip('.sav')+'.h5')
    data_file = os.path.join(IMAGE_DIR, data_path)
    full_frame = read_idl_file(data_file)
    num_quads = full_frame.shape[0]
    with h5py.File(HF_FILE_NAME, 'w') as hf:
        for i in range(0, num_quads):
            file_header = QUADS[i]
            hf.create_dataset(file_header, data=full_frame[i, :, :])
            #os.remove(HF_FILE_NAME)
            #cc
            
