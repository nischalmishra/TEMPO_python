# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:34:45 2017

@author: nmishra
"""


import os
import numpy as np

import matplotlib.pyplot as plt
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            perform_bias_subtraction,\
                                            perform_linearity_correction,\
                                            perform_smear_offset_correction,\
                                            apply_PRNU_correction
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_raw_image(full_frame):
    quad_A = full_frame[0,:,:]
    quad_B = full_frame[1,:,:]
    quad_C = full_frame[2,:,:]
    quad_D = full_frame[3,:,:]
    lower_quads = np.concatenate((quad_D, np.fliplr(quad_C)), axis=1)
    upper_quads = np.concatenate((np.flipud(quad_A), np.rot90(quad_B,2)),axis=1)
    full_frame = np.concatenate((lower_quads, upper_quads), axis=0)
    plt.imshow(full_frame, origin='lower', cmap='bwr')
    plt.xlim(0, 2112)
    plt.ylim(0, 2092)
    plt.show()
    cc
    



def main():
    """
    Read in the saved IDL variables and makes quad images and do all the data processing
    """
   
    
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.06.28\saved_quads\test'
    data_path_all = [each for each in os.listdir(file_path)
                     if each.endswith('.sav')]
    for data_path in data_path_all:
        data_file = os.path.join(file_path, data_path)
        full_frame = read_idl_file(data_file)
       
        # Ok, lets break the full_frame into quads        
        num_quads, spectral_dims, spatial_dims = full_frame.shape
        print(num_quads,spectral_dims , spatial_dims)
       
        
        quads = ['Quad D', 'Quad C', 'Quad A', 'Quad B']
#        quad_D = full_frame[0:int(spectral_dims/2), 0:int(spatial_dims/2)]
#        quad_C = full_frame[0:int(spectral_dims/2), int(spatial_dims/2):]
#        quad_A = full_frame[int(spectral_dims/2):, 0:int(spatial_dims/2)]
#        quad_B = full_frame[int(spectral_dims/2):, int(spatial_dims/2):]
        
        # Ok, let's seperate out the active quads, trailing overclocks and smear overclocks
        
        # Begin with active quads
            
        plot_raw_image(full_frame)
        
       
        # Now the trailing overclocks
        tsoc_D = full_frame[2:1030, 1034:1056]
        tsoc_C = full_frame[2:1030, 1056:1078]
        tsoc_A = full_frame[1062:2090,1034:1056]
        tsoc_B = full_frame[1062:2090, 1056:1078]

        num_coadds= 10 # you need to check the telemetry files
        integ_time = 139.9609

        # Now the smear overclocks        
         
        smear_oc_D = full_frame[1030:1046,10:1034]
        smear_oc_C = full_frame[1028:1044, 1078:2102]
        smear_oc_A = full_frame[1046:1062, 10:1034]
        smear_oc_B = full_frame[1046:1062, 1078:2102]
        
        # now for each quad, let subtract out the offset
        
        bias_subtracted_D = perform_bias_subtraction(active_D, tsoc_D)
        bias_subtracted_C = perform_bias_subtraction(active_C, tsoc_C)
        bias_subtracted_A = perform_bias_subtraction(active_A, tsoc_A)
        bias_subtracted_B = perform_bias_subtraction(active_B, tsoc_B)
        
  
        # now for each quad perfrom linearization
        linearized_D = perform_linearity_correction(bias_subtracted_D, quads[0], num_coadds)
        linearized_C = perform_linearity_correction(bias_subtracted_C, quads[1], num_coadds)        
        linearized_A = perform_linearity_correction(bias_subtracted_A, quads[2], num_coadds)
        linearized_B = perform_linearity_correction(bias_subtracted_B, quads[3], num_coadds)
        
        # now for each quad remove SMEAR
        smear_corr_D = perform_smear_offset_correction(linearized_D, integ_time)
        smear_corr_C = perform_smear_offset_correction(linearized_C, integ_time)
        smear_corr_A = perform_smear_offset_correction(linearized_A, integ_time)
        smear_corr_B = perform_smear_offset_correction(linearized_B, integ_time)
     
        # ok, now let's combine quads to one full frame again
        
        lower_quads = np.concatenate((smear_corr_D, smear_corr_C), axis=1)
        upper_quads = np.concatenate((smear_corr_A, smear_corr_B), axis=1)
        full_frame_image = np.concatenate((lower_quads, upper_quads), axis=0)
        
        # Apply PRNU Correction
        full_frame_image = apply_PRNU_correction(full_frame_image)
     
        ax = plt.gca()
        image = ax.imshow(full_frame_image, cmap='nipy_spectral', origin='lower')
        #plt.title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)            
        plt.colorbar(image, cax= cax)            
        plt.grid(False) 
        plt.show()
        
        
        ax = plt.gca()
        image = ax.imshow(linearized_D - smear_corr_D, cmap='nipy_spectral', origin='lower')
        #plt.title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)            
        plt.colorbar(image, cax= cax)            
        plt.grid(False) 
        plt.show()
        cc
        
        

        
        
        
        # now subtract out the smear using smear overclocks
        
        
        
        
        # To Do: Reallign the outlier mask, PRNU map
            
        

if __name__ == "__main__":
    main()