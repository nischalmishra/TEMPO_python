# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:49:55 2018

@author: nmishra
"""

import os
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
DATA_DIR = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC\saved_quads\saved_plots\Linearity\Mean_Stats\Variance_inlcuded'
FILE_NAME_ACTIVE_EVEN = os.path.join(DATA_DIR, 'quad_B_active_even.csv')
FILE_NAME_ACTIVE_ODD = os.path.join(DATA_DIR, 'quad_B_active_odd.csv')
INT_FILE = os.path.join(DATA_DIR, 'int_time.csv')

QUAD_NAME1 = 'Quad A (Visible Channel, Sp. Index 1:1024)'
QUAD_NAME2 = 'Quad B (Visible Channel, Sp. Index 1025:2048)'
QUAD_NAME3 = 'Quad C (UV Channel, Sp. Index 1025:2048)'
QUAD_NAME4 = 'Quad D (UV Channel, Sp. Index 1:1024)'


active_data_odd = pd.read_csv(FILE_NAME_ACTIVE_ODD, delimiter=",")
integration_time = pd.read_csv(INT_FILE, delimiter=",")
integration_time.columns = ['Int. time']
active_data_odd.columns = ['DN']
length = active_data_odd.shape[0]

#active_data_odd = active_data_odd.values
#print(tsoc_data_odd)



active_data_even = pd.read_csv(FILE_NAME_ACTIVE_EVEN, delimiter=",")
active_data_even.columns = ['DN']
integration_even = pd.read_csv(INT_FILE, delimiter=",")

#active_data_even = active_data_even.values
#print(tsoc_data_even)
#
active_DNs_odd = []
active_DNs_even = []
dframe1 = pd.DataFrame()
dframe2 = pd.DataFrame()
dframe1['Int_time'] = integration_time['Int. time'].values
dframe1['Signal_Odd'] = active_data_odd['DN'].values
dframe1['Signal_Even'] = active_data_even['DN'].values
dframe1 = dframe1.groupby(['Int_time']).mean()
print(dframe1.head)
active_DNs_odd = np.array(active_DNs_odd)
active_DNs_even = np.array(active_DNs_even)

output_file_name_odd = os.path.join(DATA_DIR, 'avg_stats/quad_B_all.csv')
#output_file_name_even = os.path.join(data_dir,'avg_stats/quad_C_active_even.csv')


dframe1.to_csv(output_file_name_odd, float_format='%1.3f', sep=",")
