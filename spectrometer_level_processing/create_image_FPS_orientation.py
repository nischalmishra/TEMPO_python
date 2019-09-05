# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 07:36:25 2019

Function to create TEMPO Image in Spectrometer Orientation

@author: nmishra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                          parse_telemetry_file,\
                                          perform_coaddition_correction

def make_full_frame_from_raw_fpe(full_frame):
    """ This function makes the full frame image including the overclocks from
    the raw binary data.

    ******Inputs***********
    1. Raw binary stream
    2. nx_quad : num of pixels in x-direction
    3. ny_quad : num of pixels in y-direction

    *******Outputs********
    full_frame[x,y]
    full frame glued together such that the spatial direction, or field angle
    direction, or position along the slit is the column #, or x direction
    AND the dispersion direction or spectral direction is the row #, or the
    y direction. This is like Fig. 5-2 of BATC SER2450865, TEMPO Processing
    Algorithm for Image Reconstruction

    Note : IDL displays columns and rows of matrix. Python diplays rows and
    colums. This function has been made longer purposely  to ensure that I
    understand how quads A, B, C and D are organized in a full frame image.

    TO DO : make the code shorter and optimize the numpy array to make it run
    faster. Not an immediate priority though.

    """
    
    quad_a = full_frame[0]
    quad_b = np.fliplr(full_frame[1])
    quad_d = np.flipud(full_frame[3])
    quad_c = np.rot90(full_frame[2], 2)
    lower = np.concatenate((quad_d, quad_c), axis=1)
    upper = np.concatenate((quad_a, quad_b), axis=1)
    full_frame = np.concatenate((upper, lower), axis=0)
    return full_frame
  
  
    
def create_image(image_data):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """
    #print(figure_name)
    #plt.figure()
    fig_ax = plt.gca()
    image_data = np.array(image_data)
    #image_data[image_data1800] = np.min(image_data)
    #image_data = np.abs(image_data)
    image = fig_ax.imshow(image_data, cmap='nipy_spectral',
                          origin='lower', interpolation='none')
    #image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          #origin='lower', interpolation='none')
    
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.show()
    #plt.pause(0.1)
    #plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    plt.close('all')

    
def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.27'
    image_dir = os.path.join(data_dir, r'saved_quads')
    image_save_dir = os.path.join(image_dir, 'FPS_orientation_raw')
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')   
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)

    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.sav')])
    print('Total data = ', len(data_path_all))

    for data_path in data_path_all:
        data_file = os.path.join(image_dir, data_path)
        print(data_path)
        
#        cc
        telemetry_file_name = data_path.split('.')[0]
       #################### Telemetry Statistics###############################
        # Now lets parse the telemetry file to get the required information###
        #Note; The function can return the CCD temp and FPE temp. But these are
        #not needed on the image processing algorithm.. atleast for the time being

        telemetry_file = os.path.join(telemetry_dir, telemetry_file_name+'.h5')
        if os.path.exists(telemetry_file):
            coadds, int_time, fpe_temp, fpa_temp = parse_telemetry_file(telemetry_dir,
                                                                        telemetry_file_name)
           # print('FPE Temp. = ', round(fpe_temp, 2), 'FPA Temp = ',
                 # round(fpa_temp, 2))
            #print('Integ. Time =' , int_time)
            #print(coadds)
            
        else:
            print('Telemetry file missing')
            coadds = 10
            int_time = 6000
            fpe_temp = 43
            fpa_temp = -21.5
        
       ##############Image Processing begins here##############################
       
        #print(coadds)
    
        full_frame = read_idl_file(data_file)
        full_frame = perform_coaddition_correction(full_frame, coadds)
        full_frame =  make_full_frame_from_raw_fpe(full_frame)
        #create_image(full_frame)
        processed_file_name = image_save_dir +'/'+ data_path+'.csv'        
        np.savetxt(processed_file_name, full_frame,fmt='%1.3f',  delimiter=",")
       
if __name__ == "__main__":
    main()        
        
        
        
        
        
        
        
        
        
        
        