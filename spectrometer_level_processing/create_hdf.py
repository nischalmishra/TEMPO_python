# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:34:45 2017

@author: nmishra
    This function creates Processed TEMPO SPECTROMETER IMAGE. All these
    routines are stored in a script called PROCESS_TEMPO_Spectrometer_Data.py and
    imported to this code via Python import. Eg. @ line 27.

    Processing Steps
    --------------------------------------------
    1. Read RAW Image files (IDL saved variable)
        Note : Dave's IDL script reads CCSDS packets faster than my Python code.
        Hence rae image files are read from the .sav files are IDL script
    2. Offset Removal (Trailing Overclocks)
    3. Non-linearity correction via Look Up Table (Options available)
    4. Smear Removal via Active Area Method (Options available)
    5. Cross-Talk Removal (Options avaialable)
    6. PRNU Correction (Options available)
    7. Create Image from Quads (Options available)
    8. Dark Current Removal (Options available)

"""
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt

#import pandas as pd

#from random import randint

#from scipy.interpolate import interp1d
#pylint: disable= E1101



def main():
    """
    This is the main function that does all the processing. It does all the analysis portion

    """
    diffuser_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\processed_corrected_linearity'
    light_data = '2018_06_27_00_38_52_570602.img.sav.csv'
    vis_lamp_data = os.path.join(diffuser_data_dir, light_data)
    vis_dark_data = os.path.join(diffuser_data_dir, 'Dark_data','2018_06_27_00_18_29_409450.img.sav.csv')
    vis_lamp_data_file = np.genfromtxt(vis_lamp_data, delimiter=',')
    vis_dark_data_file = np.genfromtxt(vis_dark_data, delimiter=',')
    print(vis_dark_data_file.shape)
    data_all  = np.round(vis_dark_data_file, 3)
    #plt.plot(vis_dark_data_file[100,:],'r')
    #plt.plot(vis_lamp_data_file[100,:],'b')
    #plt.show()
   
    hf_name = os.path.join(diffuser_data_dir,'To_Dave','TEMPO_VIS_Lamp.h5')
    #Write into h5file    
    hf = h5py.File(hf_name,'w')
    hf.create_dataset('Processed_VIS_Lamp_Bank', data=np.round(vis_lamp_data_file, 3))
    hf.create_dataset('Processed_VIS_Dark', data=np.round(vis_dark_data_file, 3))
#    hf.create_dataset('radiance_UVLamp', data=radiance_dc_corrected_UV)
#    hf.create_dataset('radiance_VISLamp', data=radiance_dc_corrected_UV)






if __name__ == "__main__":
    main()
