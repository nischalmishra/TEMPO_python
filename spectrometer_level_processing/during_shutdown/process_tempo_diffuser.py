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



import matplotlib.pyplot as plt
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            parse_telemetry_file,\
                                            perform_coaddition_correction,\
                                            perform_bias_subtraction,\
                                            apply_linearity_correction,\
                                            perform_smear_removal,\
                                            remove_cross_talk,\
                                            create_final_image
                                            #read_prnu_files,\
                                            #apply_electronics_gain
                                            #parse_prnu_file,\
                                            #calculate_mean_offset_diff,\
                                            #read_outlier_mask,\


from mpl_toolkits.axes_grid1 import make_axes_locatable

def filter_outlier_median(quads):
    """ Apart from the fixed mask, there are times when the outlier needs to be
    run in order to get good statistics. This will be used in conjunction to
    the fixed mask"""
    if np.array(quads).ndim == 3:
        ndims, nx_quad, ny_quad = quads.shape
    elif np.array(quads).ndim == 2:
        ndims = 1
        nx_quad, ny_quad = quads.shape
    else:
        nx_quad = 1
        ndims = 1
        ny_quad = len(quads)
    hist_data = np.reshape(quads, (ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    return outlier_filtered_data


def create_image_active_region(full_frame):
    """ Creates the image of active region of image
    """
#    uv_ccd = np.concatenate((full_frame[3, :, :],
#                             np.fliplr(full_frame[2, :, :])),
#                            axis=1)
#    visible_ccd = np.concatenate((np.flipud(full_frame[0, :, :]),
#                                  np.rot90(full_frame[1, :, :], 2)),
#                                 axis=1)
    image = np.concatenate((np.concatenate((full_frame[3, 2:1030, 10:1034],
                                            np.fliplr(full_frame[2, 2:1030, 10:1034])),
                                           axis=1),
                            np.concatenate((np.flipud(full_frame[0, 2:1030, 10:1034]),
                                            np.rot90(full_frame[1, 2:1030, 10:1034], 2)),
                                           axis=1)), axis=0)
    return image

def create_image(image_data, title, figure_name):
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
    plt.title(title)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close('all')


def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.13.27'
    image_dir = os.path.join(data_dir, r'saved_quads')
    save_dir = os.path.join(image_dir, 'processed_image')
    image_save_dir = os.path.join(image_dir, 'saved_plots_modified')
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')
    linearity = 1 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
#    all_image = []
#    var_all = []
#    signal_all = []
    #dark_current = 0
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
#    data_path_all = sorted([each for each in os.listdir(image_dir)
#                            if each.startswith('2017_07_30_00_24_59_38064')
#                            and  each.endswith('.sav')])
#
    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.sav')])

#    #dark_data = data_path_all[1:len(data_path_all):2]


    print('Total data = ', len(data_path_all))
#    resultFile = open (os.path.join(save_dir,'file_name_all.csv'),'w')
#
#    for results in data_path_all[0:15]:
#         resultFile.write(results+"\n")
#    resultFile.close()
    count = 0
    #peak_loc_all = []
    for data_path in data_path_all[181:]:
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
            coadds = 100
            int_time = 6000
            fpe_temp = 43
            fpa_temp = -21.5

       ##############Image Processing begins here##############################

        #print(coadds)

        full_frame = read_idl_file(data_file)
        full_frame = perform_coaddition_correction(full_frame, coadds)

        #full_frame = full_frame/coadds
        #print(np.max(full_frame))
        #cc
        #check_image = create_image_active_region(full_frame)
        #raw_image = create_final_image(np.array(full_frame))
       # print('Max. Val. Raw = ', np.max(raw_image))
        quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']

        ##########################OFFSET REMOVAL###############################
        # Input : Full_Frame TEMPO IMAGE
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
        #to store storage summation information
        # For dark data, only offset removal is needed to compute the dark current
        # For light data, additional processing such as linearity, smear and
        #cross talk is needed

        #diff_tsoc_odd_even = calculate_median_offset_diff(full_frame)
        bias_removed_quads = perform_bias_subtraction(full_frame)

        #print('Max. Val. Offset Removed  = ', np.max(bias_removed_quads))

        text1 = telemetry_file_name+'.img' +' (Raw Data)\n Int. time:' + \
                                    str(round(int_time, 3))+ \
                           'ms, Co-adds:' +str(int(coadds))+\
                           ', FPE temp:'+ str(round(fpe_temp, 1))+'C, ' + \
                           ' FPA temp: ' + str(round(fpa_temp, 1))+'C'

        # Now let us save the raw image
        raw_image_save = os.path.join(image_save_dir, 'raw_image')
        if not os.path.exists(raw_image_save):
            os.makedirs(raw_image_save)
        plot_save_dir = raw_image_save + '/'+ data_path+'.png'
        #create_image(raw_image/coadds, text1, plot_save_dir)
        #prnu_spectrometer[prnu_spectrometer < 0.9] = 0.9
        #prnu_spectrometer[prnu_spectrometer > 1.2] = 1.2
        #print(np.min(check_image))
        #create_image(check_image, text1, plot_save_dir)
        #---------------------------------------------------------------------

       #########################NON LINEARITY CORRECTION#######################
        # Input : Bias Subtracted Active regions
        # Output : Linearized Active Region Quads
        # pass the difference instead of coadds
        if linearity:
            linearized_a = apply_linearity_correction(bias_removed_quads[0, :, :],
                                                      quads[0])
            linearized_b = apply_linearity_correction(bias_removed_quads[1, :, :],
                                                      quads[1])
            linearized_c = apply_linearity_correction(bias_removed_quads[2, :, :],
                                                      quads[2])
            linearized_d = apply_linearity_correction(bias_removed_quads[3, :, :],
                                                      quads[3])
            linearized_quads = np.array([linearized_a, linearized_b,
                                         linearized_c, linearized_d])
            #print(np.mean(linearized_quads))
        else:
            linearized_quads = bias_removed_quads
            #----------------------------------------------------------------------
            ##########################SMEAR REMOVAL################################
            # Input : linearized quads ( all quads together)
            # Output : SMEAR offset corrected Quad

        #### lets' create the masked array with outlier mask################
        # The outlier mask is array of list of 4 quads.
        #print('Max. Val. Linearized  = ', np.max(linearized_quads))
        #outlier_mask = read_outlier_mask()

        # Note : all the arrays after this step are masked arrays

        if smear:
            smear_removed_quads = perform_smear_removal(linearized_quads, int_time)
        else:
            smear_removed_quads = linearized_quads
            #----------------------------------------------------------------------

            ##########################CROSS-TALK REMOVAL###########################
            # Input : smear removed quads (all quads together)
            # Output : cross talk removed quads
        #print('Max. Val. Smeared  = ', np.max(smear_removed_quads))

        if cross_talk:
            cross_talk_removed_quads = remove_cross_talk(np.array(smear_removed_quads))
        else:
            cross_talk_removed_quads = smear_removed_quads
            #----------------------------------------------------------------------
        #print('Max. Val. Cross Talked = ', np.max(cross_talk_removed_quads))

        #gain_applied_quads = apply_electronics_gain(np.array(cross_talk_removed_quads),
                                                     #diff_tsoc_odd_even)
        processed_image = create_final_image(np.array(cross_talk_removed_quads))
        #processed_image = processed_image - dark_current

        #prnu_map = parse_prnu_file()
        #prnu_map = read_prnu_files()
        #prnu_spectrometer = create_final_image(np.array(prnu_map))
        #prnu_spectrometer[prnu_spectrometer > 1.03] = 1.02
        #prnu_spectrometer[prnu_spectrometer < 0.97] = 0.98
        #outlier_spectrometer = create_final_image(np.array(outlier_mask))
        #nx_quad, ny_quad = processed_image.shape
        #outlier_spectrometer = np.reshape(outlier_spectrometer,
         #                                 (nx_quad*ny_quad, 1))
        #outlier_detectors = np.array(np.where([outlier_spectrometer == 1]))
        #print('outliers =', outlier_detectors.shape[1])
        #outlier_spectrometer = np.reshape(outlier_spectrometer, (nx_quad, ny_quad))
        #create_image(outlier_spectrometer,'outliers = '+str(outlier_detectors.shape[1]),'b')
        #processed_image = processed_image/(prnu_spectrometer)
        #processed_image[processed_image >= 0.90*16383] = np.NAN
        #processed_image = processed_image*180/int_time
        #print(np.nanmax(processed_image))
        #processed_image[processed_image>=0.85*16383] = 16383

        text1 = telemetry_file_name+'.img, ' + 'Int. time:' + str(round(int_time, 2))+ \
                           'ms\n Co-adds:' +str(int(coadds))+\
                           ', FPE temp:'+ str(round(fpe_temp, 1))+'C, ' + \
                           ' FPA temp: ' + str(round(fpa_temp, ))+'C'
        processed_image_save = os.path.join(image_save_dir, 'processed_image')
        if not os.path.exists(processed_image_save):
            os.makedirs(processed_image_save)
        plot_save_dir = processed_image_save + '/' + data_path+'.png'

        processed_image = np.array(processed_image)
        processed_image[processed_image >= 0.85*16383] = 16383
        create_image(processed_image, text1, plot_save_dir)
        # Gain is 17.2
        processed_file_name = image_save_dir +'/'+ data_path+'.csv'
        #print(processed_file_name)

        np.savetxt(processed_file_name, processed_image, fmt='%1.3f', delimiter=",")
        count = count+1


    #np.savetxt(variance_file_name, var_all,fmt='%1.3f',  delimiter=",")
    #np.savetxt(mean_file_name, signal_all,fmt='%1.3f', delimiter=",")


if __name__ == "__main__":
    main()
