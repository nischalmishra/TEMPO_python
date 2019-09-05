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
#pylint: disable= E1101
import h5py
from process_tempo_level0_V1 import read_idl_file,\
                                 read_h5_file,\
                                 parse_telemetry_file,\
                                 perform_coaddition_correction,\
                                 identify_ping_pong_fps,\
                                 identify_ping_pong_image,\
                                 perform_bias_subtraction,\
                                  apply_non_linearity_correction,\
                                 perform_temp_correction,\
                                 remove_cross_talk,\
                                 perform_smear_removal,\
                                 create_final_image
                                

def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.27'
    image_dir = os.path.join(data_dir, r'saved_quads/saved_hdf_input_files')
    image_save_dir = os.path.join(data_dir, r'saved_quads\processed_h5\Updated_Linearity')
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')
    linearity = 1 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    temp_correction = 1
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)

    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.h5')])
    print('Total data = ', len(data_path_all))

    for data_path in data_path_all[10:]:
        data_file = os.path.join(image_dir, data_path)
        print(data_path)
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
        if os.path.exists(telemetry_file):
            coadds, int_time, fpe_temp, fpa_temp = parse_telemetry_file(telemetry_dir,
                                                                        telemetry_file_name)
        else:
            print('Telemetry file missing')
            coadds = 100
            int_time = 6000
            fpe_temp = 43
            fpa_temp = -21.5
        full_frame = perform_coaddition_correction(full_frame, coadds)

        ####################PING PONG INDENTIFICATION ##########################
        # Ping (Odd) and Pong (Even) are identified based on the magnitude of
        #trailing serial overclocks

        ping_pong_phase_image = identify_ping_pong_image(full_frame)
        ping_pong_phase_fps = identify_ping_pong_fps()

        ##########################OFFSET REMOVAL###############################
        # Input : Full_Frame TEMPO IMAGE
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
        #to store storage summation information

        bias_removed_quads = perform_bias_subtraction(full_frame)        
        #print(bias_removed_quads.shape)
      

       #########################NON LINEARITY CORRECTION#######################
        # Input : Bias Subtracted Active regions
        # Output : Linearized Active Region Quads
        # Linearity correction run extremely slow.So decided to run it for each
        #quad at a time
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

        #----------------------------------------------------------------------
        ##########################TEMPO IMAGE FORMATION\#######################

        # Input : linearized quads ( all quads together)
        # Output : FULL TEMPI IMAGE ( VIS in the bottom and UV on the top)

        processed_image = create_final_image(np.array(smear_removed_quads))

        #####################CONVERSION TO COUNT RATES ########################
        # Input : PROCESSED TEMPO IMAGE and Integration time
        # Output : IMAGE RATES (DN/ms)

        processed_image = processed_image/int_time
       
        hf_file_name = os.path.join(image_save_dir, data_path.strip('.img.h5') +'_2.h5')
        with h5py.File(hf_file_name, 'w') as hf:
            hf.create_dataset('Processed_data', data=processed_image)
            hf.close()
        #print(int_time)
        #cc
        processed_file_name = image_save_dir +'/'+ data_path+'.csv'
        np.savetxt(processed_file_name, processed_image, fmt='%1.3f', delimiter=",")
        cc
if __name__ == "__main__":
    main()
