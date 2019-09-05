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
#pylint: disable= E1101
#import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
#pylint --disable=error1,error2

def create_image(image_data, int_time):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """
    #print(figure_name)
    #plt.figure()
    fig_ax = plt.gca()
    image_data = np.array(image_data)
    #image_data[image_data1800] = np.min(image_data)
    #image_data = np.abs(image_data)
    image = fig_ax.imshow(image_data, cmap='jet',
                          origin='lower', interpolation='none')
    #image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          #origin='lower', interpolation='none')
    plt.title('Spectral Registration Hg Ar Check, Int. time = ' + str(round(int_time,2)) + ' ms')
    plt.xlabel('---------> Spatial Direction', fontsize=12)
    plt.ylabel('---------> VIS to UV', fontsize=12)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.show()
    cc

    #plt.pause(0.1)
    #plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    #plt.close('all')


def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    #start_time = time.time()
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Spectral_reg_spatial_disc'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.26'
    image_dir = os.path.join(data_dir, r'saved_quads\saved_hdf_input_files')
    image_save_dir = os.path.join(data_dir, r'saved_quads\processed_h5\Jmp_discontinuity')
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')
    bias_subtraction = 1
    linearity = 1 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    temp_correction = 1
    apply_gain = 1
    int_time_all = []
    jmp_disc_all = []
    quad_d_all = []
    dframe = pd.read_excel(r'C:\Users\nmishra\Desktop\pixel_to_wavelen.xlsx',
                           sheetname='Sheet1', index_col=None, header=None)
    wavelen = dframe.values
    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.h5')])    

    for data_path in np.arange(0, len(data_path_all)):
        data_file = os.path.join(image_dir, data_path_all[data_path])
        #print(data_file)
       #################### READ INPUT IMAGE ##############################
       # Input files are in the form of .sav files from IDL.
        #It can be in any form such as #.csv or hdf.
        #Just need to write a reader script in the process tempo_level0.py

        full_frame = read_h5_file(data_file)      
        

       #################### COADDTION CORRECTION ##############################
        # Now lets parse the telemetry file to get the required information###
        #Note; The function can return the CCD temp and FPE temp which  are
        # for the image processing algorithm

        telemetry_file_name = data_path_all[data_path].split('.')[0]
        telemetry_file = os.path.join(telemetry_dir, telemetry_file_name+'.h5')
        if os.path.exists(telemetry_file):
            coadds, int_time, fpe_temp, fpa_temp = parse_telemetry_file(telemetry_dir,
                                                                        telemetry_file_name)
        else:
            print('Telemetry file missing')
            coadds = 63
            int_time = 93
            fpe_temp = 43
            fpa_temp = -21.5
            print(fpa_temp)

        #print(int_time)
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

        if bias_subtraction:
            bias_removed_quads = perform_bias_subtraction(active_area, trailing_oclk)
        else:
            bias_removed_quads = active_area
       
        #cc

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
        create_image(processed_image, int_time)      
        
        #####################CONVERSION TO COUNT RATES ########################
        # Input : PROCESSED TEMPO IMAGE and Integration time
        # Output : IMAGE RATES (DN/ms)
#        plt.figure()
        wvl = wavelen[::-1]
        plt.plot(wvl, processed_image[:, 500], 'b', label='Pixel# 500')
        #plt.axhline(y=16363, color='r', linestyle='dotted', label='14 bit saturation')
        #plt.ylim(0, 17000)
        plt.title('Processed Image Counts Vs Wavelength (int.time =' +
                  str(round(int_time, 2))+ ' ms)', fontsize=12)
        plt.ylabel('Signal Count (e-)', fontsize=12)
        plt.xlabel('Wavelength (nm)', fontsize=12)
        plt.legend(loc='best')
        plt.grid(linestyle=':')
        plt.show()
        cc

        left_quad = processed_image[:, 1000:1024]
        right_quad = processed_image[:, 1025:1049]
        ten_pixels_left = np.mean(left_quad, axis=1)
        ten_pixels_right = np.mean(right_quad, axis=1)
        pixel_loc = [1045, 1135, 1984, 2027, 2051]
        pixel_val = pixel_loc[4]        
        #for pixel_val in pixel_loc:
        #wavelen = wavelen[2010]
        quad_d_all.append(ten_pixels_right[pixel_val])
        jmp_discontinuity = 100*(ten_pixels_left[pixel_val] - ten_pixels_right[pixel_val])/\
                            ten_pixels_right[pixel_val]
        #jmp_discontinuity = (ten_pixels_left[pixel_val] - ten_pixels_right[pixel_val])
        int_time_all.append(int_time)
        jmp_disc_all.append(jmp_discontinuity)
        #file_name2 = image_save_dir +'/'+ 'jmp_disc_all_119.98ms.csv'
        #np.savetxt(file_name2, np.array(jmp_disc_all), fmt='%1.3f', delimiter=",")
       # cc


    #file_name = image_save_dir +'/'+ 'jmp_disc_480.8nm.csv'
    #np.savetxt(file_name, np.array(jmp_disc_all), fmt='%1.3f', delimiter=",")
    #file_name1 = image_save_dir +'/'+ 'int_time_all.csv'
    #np.savetxt(file_name1, np.array(int_time_all), fmt='%1.3f', delimiter=",")

    wvl = wavelen[::-1][pixel_val]
    print(pixel_val)
    plt.figure()
    plt.plot(processed_image[pixel_val, :], 'b', label='Pixel# 500 along ' + \
             str(round(wvl[0], 2)) + 'nm') 
    plt.title('Processed Image Counts Vs Spatial Indices (int.time =' +
          str(round(int_time, 2))+ ' ms)', fontsize=12)
    plt.ylabel('Signal Count (e-)', fontsize=12)
    plt.xlabel('Spatial Pixel Indices ', fontsize=12)
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.show()
    
    plt.figure()    
    plt.plot(int_time_all, jmp_disc_all, 'ro', label='Along ' + str(round(wvl[0], 2))+' nm')
    plt.title('Pct. Difference Between Quad C and Quad D Vs. In.time', fontsize=12)
    plt.ylabel('% Diff [(quadD-quadC)/quadC])', fontsize=12)
    plt.xlabel('Integration time (ms)', fontsize=12)
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.ylim(-5, 25)
    #plot_dir = r'C:\Users\nmishra\Desktop\Quad_Disc\V1'
   # fig_name = os.path.join(plot_dir,' disc_'+ str(round(wvl[0], 2))+ 'nm_int_time.png')
    #plt.savefig(fig_name, dpi=100)

    plt.figure()
    plt.plot(quad_d_all, jmp_disc_all, 'ro', label='Along ' + str(round(wvl[0], 2))+ ' nm')
    plt.axvline(x=6000, color='b', linestyle='dotted', label='6000 e-')
    plt.title('Pct. Difference Between Quad C and Quad D Vs. Signal (Quad C)', fontsize=12)
    plt.ylabel('% Diff [(quadD-quadC)/quadC])', fontsize=12)
    plt.xlabel('Electrons (e-)', fontsize=12)
    plt.legend(loc='best')
    plt.grid(linestyle=':')
    plt.ylim(-5, 25)
    plt.show()
    #cc
    #fig_name = os.path.join(plot_dir,' disc_'+ str(round(wvl[0], 2))+ 'nm_signal.png')
    #plt.savefig(fig_name, dpi=100)
    #cc

if __name__ == "__main__":
    main()
