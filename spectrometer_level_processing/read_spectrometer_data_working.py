# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:34:45 2017

@author: nmishra
    This function creates a  Processed TEMPO SPECTROMETER IMAGE.
    Processing Steps
    1. Read RAW Image files (IDL saved variable)
    2. Offset Removal
    3. Non-linearity correction
    4. Cross-Talk Removal
    5. PRNU Correction
    6. Create Image from Quads

"""


import os
import pickle
import numpy as np

import matplotlib.pyplot as plt
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            perform_bias_subtraction,\
                                            apply_linearity_correction,\
                                            perform_smear_removal,\
                                            remove_cross_talk,\
                                            create_final_image,\
                                            parse_telemetry_file,\
                                            parse_prnu_file,\
                                            apply_prnu_correction

from mpl_toolkits.axes_grid1 import make_axes_locatable

def create_image(image_data):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """
    plt.figure()
    fig_ax = plt.gca()
    image = fig_ax.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title('Difference between Ball and Locally Processed Image' )
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    plt.show()
    plt.close('all')

def create_image_both(my_image1, my_image2, my_image3):
    
    fig, axes = plt.subplots(nrows=1, ncols=1)
    im = axes.imshow(my_image1, origin='lower')
    clim=im.properties()['clim']
    fig.colorbar(im, shrink=0.8)
    plt.title('Original Image')
    plt.show()
    
    fig, axes = plt.subplots(nrows=1, ncols=1)
    axes.imshow(my_image2, origin='lower',clim=clim)
    fig.colorbar(im, shrink=0.8)
    plt.title('Processed Image + Ball PRNU')
    plt.show()
    
    fig, axes = plt.subplots(nrows=1, ncols=1)
    axes.imshow(my_image3, origin='lower',clim=clim)
    fig.colorbar(im, shrink=0.8)
    plt.title('Processed Image + Local PRNU')
    plt.show()

def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2017.06.30'
    image_dir = os.path.join(data_dir, 'saved_quads')
    save_dir = os.path.join(image_dir, 'processed_image')
    telemetry_dir = os.path.join(data_dir, 'processed/h5')
    linearity = 1 # option to turn on/off linearity
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    data_path_all = [each for each in os.listdir(image_dir)
                     if each.endswith('.sav')]
    for data_path in data_path_all:
        data_file = os.path.join(image_dir, data_path)
        telemetry_file_name = data_path.split('.')[0]

       #################### Telemetry Statistics###############################
        # Now lets parse the telemetry file to get the required information###
        #Note; The function can return the CCD temp and FPE temp. But these are
        #not needed on the image processing algorithm.. atleast for the time being

        if os.path.exists(telemetry_dir):
            coadds, int_time = parse_telemetry_file(telemetry_dir,
                                                    telemetry_file_name)
        else:
            coadds = 10
            int_time = 118

       ##############Image Processing begins here##############################

        full_frame = read_idl_file(data_file)
        uv_ccd = np.concatenate((full_frame[3, :, :],
                                 np.fliplr(full_frame[2, :, :])),
                                axis=1)
        visible_ccd = np.concatenate((np.flipud(full_frame[0, :, :]),
                                      np.rot90(full_frame[1, :, :], 2)),
                                     axis=1)
        image = np.concatenate((np.concatenate((full_frame[3, 2:1028, 10:1034],
                                np.fliplr(full_frame[2, 2:1028, 10:1034])),
                                axis=1), 
                                np.concatenate((np.flipud(full_frame[0, 2:1028, 10:1034]),
                                np.rot90(full_frame[1, 2:1028, 10:1034], 2)),
                                     axis=1)), axis=0)    
        #num_quads, spectral_dims, spatial_dims = full_frame.shape
        quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']

        ##########################OFFSET REMOVAL###############################
        # Input : Full_Frame TEMPO IMAGE
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
        #to store storage summation information

        bias_removed_quads = perform_bias_subtraction(full_frame)
        print(bias_removed_quads.shape)

        #create_image(create_final_image(np.array(full_frame)))
        #---------------------------------------------------------------------

       #########################NON LINEARITY CORRECTION#######################
        # Input : Bias Subtracted Active regions
        # Output : Linearized Active Region Quads
        if linearity:
            linearized_a = apply_linearity_correction(bias_removed_quads[0, :, :],
                                                      quads[0], coadds)
            linearized_b = apply_linearity_correction(bias_removed_quads[1, :, :],
                                                      quads[1], coadds)
            linearized_c = apply_linearity_correction(bias_removed_quads[2, :, :],
                                                      quads[2], coadds)
            linearized_d = apply_linearity_correction(bias_removed_quads[3, :, :],
                                                      quads[3], coadds)
            linearized_quads = np.array([linearized_a, linearized_b,
                                         linearized_c, linearized_d])
        else:
            linearized_quads = bias_removed_quads
        #----------------------------------------------------------------------
        #linearized_quads = bias_removed_quads

        ##########################SMEAR REMOVAL################################
        # Input : linearized quads ( all quads together)
        # Output : SMEAR offset corrected Quad
        smear_removed_quads = perform_smear_removal(linearized_quads, int_time)
        #----------------------------------------------------------------------

        ##########################CROSS-TALK REMOVAL###########################
        # Input : smear removed quads (all quads together)
        # Output : cross talk removed quads
        cross_talk_removed_quads = remove_cross_talk(smear_removed_quads)
        #----------------------------------------------------------------------
        processed_image = create_final_image(np.array(cross_talk_removed_quads))
        prnu_map = parse_prnu_file()
        final_image = processed_image/prnu_map
        
        # Now for my prnu
        
        my_image = apply_prnu_correction(processed_image)
        create_image_both(image, np.array(coadds*final_image), my_image*coadds)        
        diff = my_image-final_image
        create_image(diff)

        ######################Save the processed image variable################

        # Note : The final images is saved as python variable and read when needed
       # https://stackoverflow.com/questions/6568007/how-do-i-save-and-restore-multiple-variables-in-python

        file_name = data_path.split('.')[0]
        variable_name = save_dir +'/'+ file_name
        with open(variable_name, 'wb') as data:
            pickle.dump(processed_image, data)
        print('Hurrah!')
        create_image(np.array(processed_image))

if __name__ == "__main__":
    main()
