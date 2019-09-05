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
#from random import randint
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import uniform_filter
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            perform_bias_subtraction,\
                                            apply_linearity_correction,\
                                            perform_smear_removal,\
                                            remove_cross_talk,\
                                            create_final_image,\
                                            read_outlier_mask,\
                                            parse_telemetry_file,\
                                            read_prnu_files,\
                                            parse_prnu_file


from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

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

def create_image(image_data, max_row, text1, figure_name):
    """ This function creates the TEMPO Active Area Image. This is a placeholder
    to check the output of various image processing steps
    """

    plt.figure()
    fig_ax = plt.gca()
    print(image_data.shape)
    image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          origin='lower', interpolation='none')
    plt.title(text1)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show()
    #plt.pause(0.1)
    plt.show()
    text1 = 'Laser WL = 727.5 nm'




#    rows = image_data.shape[1]
#    fig = plt.figure(figsize=(12, 4))
#    fig.subplots_adjust(hspace=0.4, wspace=0.4)
#    ax1 = fig.add_subplot(1, 5, 1)
#    ax1.imshow(np.array(image_data), cmap='nipy_spectral', origin='lower', interpolation='none')
#    ax1.set_title(text1)
#    ax1.set_xlim(0, 15)
#    y1 = 1972
#    y2 = 2002
#    ax1.set_ylim(y1, y2)
#    ax1.set_ylabel('Spectral Pixel Index', fontsize=14)
#    fig.text(0.5, 0.04, 'Spatial Pixel Index', ha='center', fontsize=14)
#
#    ax2 = fig.add_subplot(1, 5, 2)
#    ax2.imshow(np.array(image_data), cmap='nipy_spectral', origin='lower', interpolation='none')
#    ax2.set_xlim(670, 685)
#    ax2.set_ylim(y1, y2)
#    ax2.set_title(text1)
#
#    ax3 = fig.add_subplot(1, 5, 3)
#    ax3.imshow(np.array(image_data), cmap='nipy_spectral', origin='lower', interpolation='none')
#    ax3.set_xlim(1110, 1125)
#    ax3.set_ylim(y1, y2)
#    ax3.set_title(text1)
#
#    ax4 = fig.add_subplot(1, 5, 4)
#    ax4.imshow(np.array(image_data), cmap='nipy_spectral', origin='lower', interpolation='none')
#    ax4.set_xlim(1500, 1515)
#    ax4.set_ylim(y1, y2)
#    ax4.set_title(text1)
#
#    ax5 = fig.add_subplot(1, 5, 5)
#    ax5.imshow(np.array(image_data), cmap='nipy_spectral', origin='lower', interpolation='none')
#    ax5.set_xlim(rows-15, rows)
#    ax5.set_ylim(y1, y2)
#    ax5.set_title(text1)
#    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
#    plt.show()
#    plt.close('all')




def create_image_both(my_image1, my_image2, my_image3):
    """Compares the raw image and processed image
    """

    fig, axes = plt.subplots(nrows=1, ncols=1)
    image_name = axes.imshow(my_image1, origin='lower')
    clim = image_name.properties()['clim']
    fig.colorbar(image_name, shrink=0.8)
    plt.title('Original Image')
    plt.show()

    fig, axes = plt.subplots(nrows=1, ncols=1)
    axes.imshow(my_image2, origin='lower', clim=clim)
    fig.colorbar(image_name, shrink=0.8)
    plt.title('Processed Image + Ball PRNU')
    plt.show()
    fig, axes = plt.subplots(nrows=1, ncols=1)
    axes.imshow(my_image3, origin='lower', clim=clim)
    fig.colorbar(image_name, shrink=0.6)
    plt.title('Processed Image + Local PRNU')
    plt.show(block=False)

def calculate_dark_current(data_dir, coadds):
    """ Dark current calculation. FInal dark current is the mean
    """
    all_dark_current = []
    dark_data_dir = os.path.join(data_dir, 'Dark_data')
    data_path_all = sorted([each for each in os.listdir(dark_data_dir)
                            if each.endswith('.sav')])

    for data_path in data_path_all:
        data_file = os.path.join(dark_data_dir, data_path)
        full_frame = read_idl_file(data_file)
        bias_removed_quads = perform_bias_subtraction(full_frame)
        dark_current_image = create_final_image(np.array(bias_removed_quads))
        all_dark_current.append(dark_current_image)
    all_dark_current = np.array(all_dark_current)
    #mean_dark_current = np.mean(np.mean(all_dark_current, axis=0)/coadds)

    return np.mean(all_dark_current, axis=0)

def read_wavelength_file(file_path, wavelen_val):
    """ Parse wavelength values from the log file
    """
    wavelen_data = pd.read_csv(os.path.join(file_path, 'wavelength_spectrum_analyzer.csv'),
                               delimiter=",")
    return wavelen_data[wavelen_val].values




def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass\UV\297.8_nm'
    image_dir = os.path.join(data_dir, 'saved_quads')
    save_dir = os.path.join(image_dir, 'processed_image')
    image_save_dir = os.path.join(image_dir, 'saved_plots')
    telemetry_dir = os.path.join(data_dir, 'processed/h5')
    wavelen_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\Spectral_Band_pass'
    wavelen_val = data_dir.strip('_')[-8:]

    linearity = 0 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    bandpass_val = []
    bandpass_val_norm = []
    normalization_factor_all = []
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
#    data_path_all = sorted([each for each in os.listdir(image_dir)
#                            if each.startswith('2017_07_30_00_24_59_38064')
#                            and  each.endswith('.sav')])
#
    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('2017_07_28_17_07_51_683839.img.sav')])
#
#    #dark_data = data_path_all[1:len(data_path_all):2]

    wavelength = read_wavelength_file(wavelen_dir, wavelen_val)

    print('Total data = ', len(data_path_all))
#    resultFile = open (os.path.join(save_dir,'file_name_all.csv'),'w')
#
#    for results in data_path_all:
#         resultFile.write(results+"\n")
#    resultFile.close()
    dark_current = calculate_dark_current(image_dir, coadds=1)
    count = 0
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
            print('FPE Temp. = ', round(fpe_temp, 2), 'FPA Temp = ',
                  round(fpa_temp, 2))
        else:

            print('Telemetry file missing')
            coadds = 10
            int_time = 3500
            fpe_temp = 43
            fpa_temp = -21

       ##############Image Processing begins here##############################
        full_frame = read_idl_file(data_file)
        check_image = create_image_active_region(full_frame)
        raw_image = create_final_image(np.array(full_frame))
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


        bias_removed_quads = perform_bias_subtraction(full_frame)

        #print('Max. Val. Offset Removed  = ', np.max(bias_removed_quads))

        text1 = telemetry_file_name+'.img' +' (Raw Data)\n Int. time:' + str(round(int_time, 1))+ \
                           'ms, Co-adds:' +str(int(coadds))+\
                           ', FPE temp:'+ str(round(fpe_temp, 1))+'C, ' + \
                           ' FPA temp: ' + str(round(fpa_temp, 1))+'C'
        # Now let us save the raw image
        raw_image_save = os.path.join(image_save_dir, 'raw_image')
        if not os.path.exists(raw_image_save):
            os.makedirs(raw_image_save)
        plot_save_dir = raw_image_save + '/'+ data_path+'.png'
        #prnu_spectrometer[prnu_spectrometer < 0.9] = 0.9
        #prnu_spectrometer[prnu_spectrometer > 1.2] = 1.2
        #create_image(check_image/coadds, text1, plot_save_dir)
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
        outlier_mask = read_outlier_mask()

        # Note : all the arrays after this step are masked arrays

        if smear:
            smear_removed_quads = perform_smear_removal(linearized_quads, int_time,
                                                        outlier_mask)
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
        processed_image = create_final_image(np.array(cross_talk_removed_quads))
        processed_image = processed_image - (dark_current)
        #prnu_map = parse_prnu_file() # BATC Map
        #prnu_spectrometer = prnu_map

        #Uncomment the lines below if my PRNU is to be used.
        #prnu_map = read_prnu_files() # LaRC map
        #prnu_spectrometer = create_final_image(np.array(prnu_map))
        #prnu_spectrometer[prnu_spectrometer > 1.03] = 1.02
        #prnu_spectrometer[prnu_spectrometer < 0.97] = 0.98
        #create_image(prnu_map, 'TEMPO PRNU Map', 'a')
        #create_image(prnu_spectrometer, 'TEMPO PRNU Map', 'a')
        outlier_spectrometer = create_final_image(np.array(outlier_mask))

        #create_image(outlier_spectrometer,'Outliers = 504', 'b')
        #cc
        nx_quad, ny_quad = processed_image.shape
        outlier_spectrometer = np.reshape(outlier_spectrometer,
                                          (nx_quad*ny_quad, 1))
#        outlier_detectors = np.array(np.where([outlier_spectrometer == 1]))
        #print('outliers =', outlier_detectors.shape[1])
        outlier_spectrometer = np.reshape(outlier_spectrometer, (nx_quad, ny_quad))
        #create_image(outlier_spectrometer,'outliers = '+str(outlier_detectors.shape[1]),'b')
        #processed_image = processed_image/(prnu_spectrometer)
        #processed_image = processed_image/(coadds*prnu_spectrometer)
        processed_image[processed_image >= 0.9*16383] = 16383
#        text1 = 'WL = ' + str(wavelength[count])+ 'nm'
        processed_image_save = os.path.join(image_save_dir, 'processed_image')
        if not os.path.exists(processed_image_save):
            os.makedirs(processed_image_save)
#        plot_save_dir = processed_image_save + '/' +str(wavelength[count])+ \
#                       '_' + str(count)+'_nm' +'.png'

        #print('min_val=', np.min(processed_image))
        #processed_image = np.reshape(processed_image, (nx_quad*ny_quad, 1))
        #processed_image[outlier_spectrometer] = 16383

        #col_index = np.sort([randint(100, 2000) for p in range(0, 10)])
        row_mean = np.mean(processed_image, axis=1)
        max_row_each = np.where(row_mean == max(np.mean(processed_image, axis=1)))[0]
        max_row_each = max_row_each[0]
        max_row = max_row_each
        print(max_row)
        max_row = 29       
        create_image(processed_image, max_row, text1, plot_save_dir)
#        print(max_row_each)
#        cc

        max_row = 29 # for 297.8 nm
        #max_row = 1546 #for 640nm
        #max_row = 1369 #for 605nm
        #max_row = 1050 #for 605nm
       # max_row = 142 # for 320 nm
       # max_row = 319  # for 355 nm
       # max_row = 496  # for 390 nm
       # max_row = 673  # for 425 nm
        #max_row = 849  # for 460 nm
        #max_row = 992  # for 488.2 nm
        #max_row = 192 # for 330 nm
        #max_row = 91 # for 310 nm
        #max_row = 1724 # 675 nm
        #max_row = 2032 # for 736 nm
        #max_row = 1989 # for 727.5 nm
        #check_val = 1392

        #subset_image = processed_image[max_row,  :]
        subset_image_normalized = processed_image[max_row-18 : max_row+18, :]
        normalization_factor = np.sum(subset_image_normalized[subset_image_normalized < 16383])
        normalization_factor_all.append(normalization_factor)
        #normalization_factor_subset = np.sum(subset_image_normalized \
                                              #[subset_image_normalized<16383])
        #normalization_factor_all.append(normalization_factor_subset)
        bandpass_val.append(processed_image[max_row, :])
        bandpass_val_norm.append(processed_image[max_row, :]/normalization_factor)
        count = count+1
        #print(processed_image.shape)
       # cc


#        plt.plot(processed_image[max_row-4, :], 'purple', label='4 lines down')
#        plt.plot(processed_image[max_row-3, :], 'y-', label='3 lines down')
#        plt.plot(processed_image[max_row-2, :], 'g-', label='2 lines down')
#        plt.plot(processed_image[max_row-1, :], 'r-', label='1 line down')
#        plt.plot(processed_image[max_row, :], 'k-', label='max (296.54 nm), '\
#                  + 'row:' +str(max_row))
#        plt.plot(processed_image[max_row+1, :], 'm-', label='1 line up')
##        plt.plot(processed_image[max_row+2, :], 'b-', label='2 lines up')
##        plt.plot(processed_image[max_row+3, :], 'orange', label='3 lines up')
##        plt.plot(processed_image[max_row+4, :], 'cyan', label='4 lines up')
#        plt.title('Spatial Profile along illuminated rows (Laser WL: 727.5 nm)')
#        plt.ylabel('Counts (DN)')
#        plt.xlabel('Spatial Pixel Index')
#        plt.grid(True, linestyle=':')
#
#        #plt.yscale('log')
#        #plt.xlim(0, 2000)
#
#
#        #plt.legend(loc='best')
#        plt.show()
#        cc
#        #plt.ylim(2000, 40000)
#        plt.show(block=False)
#        #cc
#        plt.figure(3)
#        for i in col_index:
#
#            plt.plot(processed_image[:, i],
#                     'o--', label='Col. ' + str(i))
#            plt.legend(loc='best', ncol=4)
#        plt.grid(True, linestyle=':')
#        plt.title('Profile along Spectral direction (Spectral Bandpass 297 nm)')
#        plt.xlabel('Spectral Pixel Index')
#        plt.ylabel('Counts (DN)')
#        #plt.xlim(655, 675)
#        #plt.yscale('log')
#
#        plt.show()
#        cc

    save_bandpass = save_dir +'/'+'bandpass_radiance_second_iteration.csv'
    save_bandpass_norm = save_dir +'/'+'bandpass_radiance_normalized_second_iteration.csv'
    norm_factor = save_dir +'/'+'integrated_radiance_second_iteration.csv'

    np.savetxt(save_bandpass, np.array(bandpass_val), delimiter=",")
    np.savetxt(save_bandpass_norm, np.array(bandpass_val_norm), delimiter=",")
    np.savetxt(norm_factor, np.array(normalization_factor_all), delimiter=",")
    #cc
if __name__ == "__main__":
    main()
