# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 09:42:05 2017

@author: nmishra
"""
import os
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import pickle

def find_saturated_pixels(quads):
    """ This function identifies the saturated pixels in the quads. 
    Ideally saturation level is 2^N-1. But for additional threshold has been    
    added for analysis purpose

    """
    #print(np.max(quads))
    nx_quad = quads.shape[0]
    ny_quad = quads.shape[1]
    data = np.reshape(quads, (nx_quad*ny_quad, 1))
    
    
    saturated_pixels = np.where(data == 2**14-1)[0]
    #print(len(saturated_pixels))
    print(saturated_pixels[0:10])
    input("Press enter to exit")
    
    
    mask = np.zeros((nx_quad * ny_quad, 1)).astype(int)
    mask[saturated_pixels] = 1
    mask = np.reshape(mask, (nx_quad, ny_quad))
    sat_mask = deepcopy(mask)
    
   
    return  sat_mask


def main():    
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    if 'Integration_Sweep' in file_path1:
    #if 'Bar_target' in file_path1:
        saturated_collects = [ 'FT6_LONG_INT_130018.dat.sav','FT6_SHORT_INT_0.dat.sav',
                              'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']

        all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
   
  
        nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files] # for cross talk image
    
  # # nominal_int_files = 'SPI_20160912100114585_ACLMP-0X2E-5_DCLMP-0X2E-5_ASMPL-0X8-7_DSMPL-0X8-7_FT6_SHORT_INT_15009.dat.sav'                     
        save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Saturation_pixels'
    
    for i in range(0, 4): # for the 4 quads        
              
        for data_files in nominal_int_files:
            print(data_files)
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)
            all_full_frame = IDL_variable.q
            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                string1 = 'VA_'
                string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '           
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            color = ['red','blue', 'green','magenta']
            quad_full_frame = all_full_frame[:, i, :, :]
            avg_quad = np.mean(quad_full_frame[:, :, :], axis=0) 
            active_quad = avg_quad[4:1028, 10:1034]
            saturated_pixels = find_saturated_pixels(active_quad)
            variable_name = save_dir+'/'+ string1 + str(int_time)
            with open(variable_name, 'wb') as data:
                pickle.dump(saturated_pixels, data)
                
            
if __name__ == "__main__":
    main()           