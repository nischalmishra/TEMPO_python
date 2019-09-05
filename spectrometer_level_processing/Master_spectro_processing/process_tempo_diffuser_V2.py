# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:34:45 2017

@author: nmishra
    This function creates Processed TEMPO SPECTROMETER IMAGE. All the functinalities
    are callled from a seperate script called PROCESS_TEMPO_Spectrometer_Data.py and
    imported to this code via Python import. Eg. @ line 39. Input image is expected
    in the from of hdf file and quads are arranged in the order of readout.
    If you are having to deal with the TEMPO  data packets in CCSDS format (.img)
    extension, run extract_spectrometer_data first and it create a raw hdf file
    in the deisred format.
    quadA = image[0, :, :]
    quadB = image[1, :, :]
    quadC = image[2, :, :]
    quadD = image[3, :, :]
    Although, most of the image processing is mosstly based off octants, each image
    processing steps the inputs are in the form of quadrants. Light Dark data are
    treated the same way as the normal 'Light' Data

    Processing Steps
    --------------------------------------------
    1. Read RAW Image files (hdf5 files)
    2. Coaddition Correction
    3. Ping-Pong Identification
    3. Offset Removal (Trailing Overclocks)
    4. Non-linearity correction via Look Up Table )
    5. FPE Temp Correction
    6. Cross-Talk Removal
    7. Smear Removal via Active Area Method
    8. Gain Application

"""
import os
import time

import numpy as np
#pylint: disable= E1101
#pylint: disable-msg=R0912
#pylint: disable-msg=R0914
#pylint: disable-msg=R0915
import h5py
from process_tempo_level0_v2 import read_h5_file,\
                                    parse_telemetry_file,\
                                    perform_coaddition_correction,\
                                    identify_ping_pong_fps,\
                                    identify_ping_pong_image,\
                                    perform_bias_subtraction,\
                                    apply_non_linearity_correction,\
                                    perform_temp_correction,\
                                    remove_cross_talk,\
                                    perform_smear_removal,\
                                    perform_gain_correction,\
                                    create_final_image

def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    start_time = time.time()
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_UV_Lamp'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.26'
    image_dir = os.path.join(data_dir, r'saved_quads\saved_hdf_input_files\Dark_Data')
    image_save_dir = os.path.join(data_dir, r'saved_quads\processed_h5\Dark_Data')
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')
    linearity = 1 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    temp_correction = 1
    apply_gain = 1
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)

    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.h5')])
   # print('Total data = ', len(data_path_all)) # sanity check

    for data_path in data_path_all:
        data_file = os.path.join(image_dir, data_path)
        #print(data_path) # sanity check
       #################### READ INPUT IMAGE ##############################
       # Input files are in the form of .sav files from IDL.
        #It can be in any form such as #.csv or hdf.
        #Just need to write a reader script in the process tempo_level0.py

        full_frame = read_h5_file(data_file)

       #################### COADDTION CORRECTION ##############################
        # Now lets parse the telemetry file to get the required information###
        #Note; The function can return the CCD temp and FPE temp which  are
        # for the image processing algorithm

        telemetry_file_name = data_path.split('.')[0]
        telemetry_file = os.path.join(telemetry_dir, telemetry_file_name+'.h5')
        #print(telemetry_file)
        
        if os.path.exists(telemetry_file):
            coadds, int_time, fpe_temp, fpa_temp = parse_telemetry_file(telemetry_dir,
                                                                        telemetry_file_name)
        else:
            print('Telemetry file missing')
            coadds = 63
            int_time = 93
            fpe_temp = 43
            fpa_temp = -21.5
            #print(fpa_temp)
        #cc
        full_frame = perform_coaddition_correction(full_frame, coadds)
        active_area = full_frame[:, 2:1030, 10:1034]
        trailing_oclk = full_frame[:, 2:1030, 1034:1056]

        ####################PING PONG INDENTIFICATION ##########################
        # Ping (Odd) and Pong (Even) are identified based on the magnitude of
        #trailing serial overclocks for the image and values from FPS testing

        ping_pong_phase_image = identify_ping_pong_image(trailing_oclk)
        ping_pong_phase_fps = identify_ping_pong_fps()

        ##########################OFFSET REMOVAL###############################
        # Input : Active Area and trailing overclocks
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024.

        bias_removed_quads = perform_bias_subtraction(active_area, trailing_oclk)
        #print(bias_removed_quads.shape)


       #########################NON LINEARITY CORRECTION#######################
        # Input : Bias Subtracted Active regions, ping pong phase. Look up table
        #is read from the function itself
        # Output : Linearized Active Region Quads

        if linearity:
            linearized_quads = apply_non_linearity_correction(bias_removed_quads,
                                                              ping_pong_phase_image,
                                                              ping_pong_phase_fps)
        else:
            linearized_quads = bias_removed_quads

        #########################FPE TEMP CORRECTION ##########################
        if temp_correction:
            temp_corrected_quads = perform_temp_correction(linearized_quads,
                                                           fpe_temp,
                                                           ping_pong_phase_image,
                                                           ping_pong_phase_fps)
        else:
            temp_corrected_quads = linearized_quads

        ##########################SMEAR REMOVAL################################
        # Input : linearized quads ( all quads together)
        # Output : SMEAR offset corrected Quad

        if cross_talk:
            cross_talk_removed_quads = remove_cross_talk(temp_corrected_quads)
        else:
            cross_talk_removed_quads = temp_corrected_quads

        if smear:
            smear_removed_quads = perform_smear_removal(np.array(cross_talk_removed_quads),
                                                        int_time)
        else:
            smear_removed_quads = cross_talk_removed_quads


        if apply_gain:
            gain_applied_quads = perform_gain_correction(smear_removed_quads,
                                                         ping_pong_phase_image,
                                                         ping_pong_phase_fps)
        else:
            gain_applied_quads = smear_removed_quads

        #----------------------------------------------------------------------
        ##########################TEMPO IMAGE FORMATION\#######################

        # Input : linearized quads (all quads together)
        # Output : FULL TEMPO IMAGE (VIS in the bottom and UV on the top)

        processed_image = create_final_image(np.array(gain_applied_quads))

        #####################CONVERSION TO COUNT RATES ########################
        # Input : PROCESSED TEMPO IMAGE and Integration time
        # Output : IMAGE RATES (DN/ms)

        processed_image = processed_image/int_time

        hf_file_name = os.path.join(image_save_dir, data_path.strip('.img.h5') +'.h5')
        with h5py.File(hf_file_name, 'w') as hf_id:
            hf_id.create_dataset('Processed_data', data=processed_image)
            hf_id.close()
        elapsed = time.time() - start_time
        print(elapsed)
        #cc
        #processed_file_name = image_save_dir +'/'+ data_path+'_V1.csv'
        #np.savetxt(processed_file_name, processed_image, fmt='%1.3f', delimiter=",")
        #cc

if __name__ == "__main__":
    main()
