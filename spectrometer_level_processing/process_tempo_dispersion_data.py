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
import pandas as pd
#from scipy.ndimage import uniform_filter
#from random import randint
from scipy.optimize import curve_fit
#from scipy.interpolate import interp1d
#pylint: disable= E1101



import matplotlib.pyplot as plt
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
    print(figure_name)
    plt.figure()
    fig_ax = plt.gca()
    print(image_data.shape)
    image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          origin='lower', interpolation='none')
    plt.title(title)
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



def calculate_dark_current(data_dir, coadds):

    """ Function to calculate the dark current. FInal dark current is the average
    """
    all_dark_current = []
    print(coadds)
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



def moving_average(data, num_points=10):
    """ Calculate the moving average of sets of data
    """
    ret = np.cumsum(data, dtype=float)
    ret[num_points:] = ret[num_points:] - ret[:-num_points]
    return ret[num_points - 1:] / num_points

def gauss_function(x_val, a_val, x0_val, sigma):
    """ Fit Gaussian function
    """
    return a_val*np.exp(-(x_val-x0_val)**2/(2*sigma**2))


def full_width_half_max(wavelen, data):
    """ Calculated FWHM
    """
    half_max = max(data) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - data[])
    inflection_point = np.sign(half_max - np.array(data[0:-1])) - \
                              np.sign(half_max - np.array(data[1:]))
    #plt.plot(wavelen[0:3839], inflection_point) #if you are interested
    #plt.show()
    #cc
    #find the left and right most indexes
    left_idx = np.where(inflection_point > 0)[0]
    #print(wavelen[left_idx])
    right_idx = np.where(inflection_point < 0)[-1]
    if wavelen[right_idx] - wavelen[left_idx] is None:
        fwhm = 0
    else:
        fwhm = wavelen[right_idx]-wavelen[left_idx]
        fwhm = fwhm[0]
    return fwhm


def find_nearest(array, value):
    """ Find the nearest neighbor of the given list element
    """
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def fit_gaussian_func(processed_image, wavelen):
    """ Function to fit the Gaussian fit to the wavelength Vs. Amplitude Data
    """
    num_data = processed_image.shape[1] # not all columns are iluuminated by laser.
    # Hence loop starts at 4 to total cols-7

    row_mean = np.mean(processed_image, axis=1)
    max_row = np.where(row_mean == max(np.mean(processed_image, axis=1)))[0]
    max_row = max_row[0]
    #print(max_row)
    # Now fit Gaussian
    start_index = max_row-5
    end_index = max_row+5

    estimates_gauss = []
    all_fwhm = []
    dframe_gauss = []

    for i in range(4, num_data-7): # there are 5 pin holes
        radiance_data = processed_image[start_index:end_index, i]
        #print(radiance_data)
        #radiance_data= uniform_filter(radiance_data, size=1, mode='constant')
        #cols = radiance_data.shape[0]
        pixel_loc = np.arange(start_index, end_index)

        # sometime curve can peak outisde the illuminated rows because of outliers.
        #Hence outlier filter has to be removed
        max_val = np.where(radiance_data == max(radiance_data))[0]
        max_val = max_val[0]
        if max_val not in [4, 5, 6, 7]:
            radiance_data[max_val] = np.min(max_val)

        radiance_data = radiance_data / np.max(radiance_data)
        #n = len(pixel_loc)
       # print(i)
       # mean_val = sum(radiance_data * pixel_loc)/n
        #print(mean_val)
        #sigma_val =  np.sqrt(sum(radiance_data * (pixel_loc - mean_val)**2) / n)
        mean_val = np.mean(pixel_loc)
        sigma_val = np.std(pixel_loc)


        best_vals, pcov = curve_fit(gauss_function, pixel_loc, radiance_data,
                                    p0=[1, mean_val, sigma_val])
        print(pcov)
        estimates_gauss.append(best_vals)
        sampled_pixel_loc = np.arange(min(pixel_loc), max(pixel_loc), 0.01)
        fitted_data = gauss_function(sampled_pixel_loc, *best_vals)
        fitted_data = fitted_data/np.max(fitted_data)
        fwhm = full_width_half_max(sampled_pixel_loc, fitted_data)
        #print(FWHM)
        all_fwhm.append(fwhm)

        plt.plot(pixel_loc, radiance_data, 'ro--', label='Original Data')
        plt.plot(sampled_pixel_loc, fitted_data, 'b', label='Fitted Data')
        plt.title('Example of Gaussian Fit\n Spatial Index: '+ str(i+1)+
                  ', Laswer WL = ' + str(wavelen)+' nm')
        plt.grid(True, linestyle=':')
        plt.xlabel('Spectral Pixel Index')
        plt.ylabel('Normalized Counts')
        plt.legend(loc='best')
       # print(best_vals)
        plt.text(max_row-5.2, 0.8, 'Peak Location = '+ str(round(best_vals[1], 2)),
                 fontsize=10, fontweight='bold', style='italic',
                 bbox={'facecolor':'Orange', 'alpha':0.6, 'pad':4})
        plt.show()

    all_fwhm = np.array(all_fwhm)
    estimates_gauss = np.array(estimates_gauss)
    dframe_gauss = pd.DataFrame({'Peak Amp.' : estimates_gauss[:, 0],
                                 'Peak Loc.' : estimates_gauss[:, 1],
                                 'Sigma1' : np.abs(estimates_gauss[:, 2]),
                                 'FWHM' : all_fwhm
                                })

    #print(all_peak_gauss)
    return dframe_gauss



def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Dispersion'
    image_dir = os.path.join(data_dir, 'saved_quads')
    save_dir = os.path.join(image_dir, 'processed_image')
    image_save_dir = os.path.join(image_dir, 'saved_plots')
    telemetry_dir = os.path.join(data_dir, 'processed/h5')
    linearity = 0 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    #bandpass_val = []
    #bandpass_val_norm = []
    #normalization_factor_all = []
    dark_current = 0
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
    if dark_current:
        dark_current = calculate_dark_current(image_dir, coadds=1)
    else:
        dark_current = 0

    count = 0
    #peak_loc_all = []
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
            int_time = 6000
            fpe_temp = 43
            fpa_temp = -21.5

       ##############Image Processing begins here##############################
        full_frame = read_idl_file(data_file)
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
        #print(np.min(check_image))
        #create_image(check_image, text1, plot_save_dir)
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
        processed_image = processed_image - dark_current

        #prnu_map = parse_prnu_file()
        prnu_map = read_prnu_files()
        prnu_spectrometer = create_final_image(np.array(prnu_map))
        prnu_spectrometer[prnu_spectrometer > 1.03] = 1.02
        prnu_spectrometer[prnu_spectrometer < 0.97] = 0.98

        outlier_spectrometer = create_final_image(np.array(outlier_mask))

        nx_quad, ny_quad = processed_image.shape
        outlier_spectrometer = np.reshape(outlier_spectrometer,
                                          (nx_quad*ny_quad, 1))
        #outlier_detectors = np.array(np.where([outlier_spectrometer == 1]))
        #print('outliers =', outlier_detectors.shape[1])
        outlier_spectrometer = np.reshape(outlier_spectrometer, (nx_quad, ny_quad))
        #create_image(outlier_spectrometer,'outliers = '+str(outlier_detectors.shape[1]),'b')
        processed_image = processed_image/(coadds*prnu_spectrometer)
        #print(np.min(processed_image))

        # Read the laser wavelength
        wavelen = np.loadtxt(os.path.join(image_dir, 'Wavelength.csv'), delimiter=',')

        # Read the spectral index
#        spec_index = np.loadtxt(os.path.join(image_dir, 'Spectral_column.csv'), delimiter=',')
        wavelen = wavelen[count]
#        spec_pix = int(spec_index[count])
        #print('Test Data= ',  wavelen, spec_pix)


                #processed_image[processed_image>=0.85*16383] = 16383
        text1 = telemetry_file_name+'.img' +' (Laser WL:' + \
                           str(wavelen) +' nm)\n Int. time:' + str(coadds)+ \
                           'ms, Co-adds:' +str(int(coadds))+\
                           ', FPE temp:'+ str(round(fpe_temp, 1))+'C, ' + \
                           ' FPA temp: ' + str(round(fpa_temp, ))+'C'
        processed_image_save = os.path.join(image_save_dir, 'processed_image')
        if not os.path.exists(processed_image_save):
            os.makedirs(processed_image_save)
        plot_save_dir = processed_image_save + '/' + data_path+'.png'
        #processed_image = uniform_filter(processed_image, size=(5, 3), mode='mirror')
        create_image(processed_image, text1, plot_save_dir)


        # Now let's find the peak for each monochromatic laser wavelengths
        dframe_gauss = fit_gaussian_func(processed_image, wavelen)
        print(dframe_gauss.head())
        #file_name = os.path.join(save_dir, 'Gaussian_parameters_' + \
                                  #str(count)+'_'+ str(wavelen)+'.csv')
        #dframe_gauss.to_csv(file_name)
        count = count+1

if __name__ == "__main__":
    main()
