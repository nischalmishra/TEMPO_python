# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:49:31 2018

@author: nmishra
"""

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            perform_bias_subtraction,\
                                            apply_linearity_correction,\
                                            perform_smear_removal,\
                                            remove_cross_talk,\
                                            create_final_image,\
                                            read_outlier_mask,\
                                            parse_telemetry_file,\
                                            read_prnu_files
                                            #parse_prnu_file,\


def make_quads_from_fits_file(file_path, data_path):
    data_path_name_split = data_path.split('_')
    int_time = data_path_name_split[-4]
    temp = data_path_name_split[-2]
    data_file = os.path.join(file_path, data_path)
    hdulist = fits.open(data_file)
    full_frame = hdulist[0].data
    return full_frame, int_time, temp


def make_full_frame_image(file_path, CCD_orientation):
    file_path_1 = file_path+'/'+CCD_orientation[0]
    file_path_2 = file_path+'/'+CCD_orientation[1]
    all_four_quads = [ ]
    
    wavelength = '740nm'
    data_file_name_1 = 'FT6BS-017_01222016_130535_QE_'+wavelength+'_int_500_temp_253K_4.fits'
    data_file_name_2= 'FT6BS-014_01212016_124519_QE_'+wavelength+'_int_500_temp_253K_4.fits'
    
    data_path_1 = os.path.join(file_path_1, data_file_name_1)
    full_frame, int_time, temp = make_quads_from_fits_file(file_path_1, data_path_1)
    total_coadds = full_frame.shape[0]
    quad_A = full_frame[0:total_coadds:2, :, :]
    avg_quad_A = np.mean(quad_A, axis=0)
    all_four_quads.append(avg_quad_A)    
    quad_B = full_frame[1:total_coadds:2, :, :]    
    avg_quad_B = np.mean(quad_B, axis=0)
    #create_image(avg_quad_B)
    all_four_quads.append(avg_quad_B)
    full_frame = avg_quad_A = avg_quad_B = None
    
    data_path_2 = os.path.join(file_path_2, data_file_name_2)
    full_frame, int_time, temp = make_quads_from_fits_file(file_path_2, data_path_2)
    total_coadds = full_frame.shape[0]
    quad_C = full_frame[0:total_coadds:2, :, :]    
    avg_quad_C = np.mean(quad_C, axis=0)
    all_four_quads.append(avg_quad_C)
    #create_image(avg_quad_C)
    quad_D = full_frame[1:total_coadds:2, :, :]    
    avg_quad_D = np.mean(quad_D, axis=0)
    all_four_quads.append(avg_quad_D)
    return all_four_quads
    

def create_image(image_data):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """
    plt.figure()
    fig_ax = plt.gca()
    print(np.max(image_data))
    image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          origin='lower', interpolation='none')
   
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show()
    #plt.pause(0.1)
    plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    #plt.close('all')

def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image
    """
    quad_a = full_frame[0, :, :]
    quad_b = full_frame[1, :, :]
    quad_c = full_frame[2, :, :]
    quad_d = full_frame[3, :, :]
    uv_ccd = np.concatenate((quad_d, np.fliplr(quad_c)),
                            axis=1)
    visible_ccd = np.concatenate((np.flipud(quad_a), np.rot90(quad_b, 2)),
                                 axis=1)
    processed_image = np.concatenate((uv_ccd, visible_ccd), axis=0)
    return processed_image

def main():
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\FPS_QE'
    CCD_orientation = [ r'FT6BS-017_QE', r'FT6BS-014_QE']
    #Note: FT6BS-017_QE : Visible  FT6BS-017_QE : UV
    all_four_quads = make_full_frame_image(file_path, CCD_orientation)   
    #spectro_image = create_final_image(np.array(all_four_quads))
    bias_subtracted_quad = perform_bias_subtraction(all_four_quads)
    create_final_image(bias_subtracted_quad)
    
    
    
    
    

if __name__ == "__main__":
    main()
