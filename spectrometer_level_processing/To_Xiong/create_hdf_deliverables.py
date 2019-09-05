# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:55:43 2019

@author: nmishra
"""


import os
import datetime
import numpy as np
import h5py

DELIVERABLE_DIR = r'C:\Users\nmishra\Desktop\SAO_deliverable'
DELIVERABLE_FILE_NAME = [each for each in os.listdir(DELIVERABLE_DIR)
                         if each.endswith('.csv')]

#print(*DELIVERABLE_FILE_NAME, sep='\n')

META_DATA = {'Bad Pixel Mask': 'Unitless [1 is bad pixel, O is good pixel]',
             'Cross Talk Factor': 'Scaling Factor (Unitless)',
             'FPA Ref Temp' : 'Temp (Deg C)',
             'FPA gain scaling':'Temp (Deg C) and Scaling Factor (Unitless)',
             'FPA DC Scaling' : 'Scaling Factor (Unitless)',
             'FPA Gain Scaling' : 'Scaling Factor (Unitless)',
             'FPE Ref Temp' : 'Temp (Deg C)',
             'Electronics Gain': 'e-/DN',
             'Non-linearity correcttion': 'DN',
             'Wavelength': 'nm',
             'Radiance Cal. Coeffs':'ph.sec/sec.cm2.nm.sr.DN',
             'PRNU': 'Unitless (Scaling Factor)',
             'Trailing Overclocks': 'DN'
            }

TODAYS_DATE = datetime.datetime.now().strftime('%m_%d_%Y')
HF_NAME = os.path.join(DELIVERABLE_DIR, 'TEMPO_LaRC_Image_Processing_' + TODAYS_DATE + '_V2.h5')
with h5py.File(HF_NAME, 'w') as hf:
    for data_files in DELIVERABLE_FILE_NAME:
        file_header = data_files.strip('.csv')
        deliverable = os.path.join(DELIVERABLE_DIR, data_files)
        data_read = np.genfromtxt(deliverable, delimiter=',')
       # print(data_read[1:3])
        hf.create_dataset(file_header, data=data_read)
        hf.attrs.update(META_DATA)
print('DONE')
