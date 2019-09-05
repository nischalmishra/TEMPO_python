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
import pandas as pd


import matplotlib.pyplot as plt
from Process_TEMPO_Spectrometer_Data import read_idl_file,\
                                            perform_coaddition_correction,\
                                            perform_bias_subtraction,\
                                            perform_temp_correction,\
                                            apply_linearity_correction,\
                                            perform_smear_removal,\
                                            remove_cross_talk,\
                                            create_final_image,\
                                            read_outlier_mask,\
                                            parse_telemetry_file,\
                                            read_prnu_files,\
                                            apply_electronics_gain
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
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\FPA_Gain_vs_Temp'
    telemetry_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2016.09.19'    
    image_save_dir = r'C:\Users\nmishra\Workspace\TEMPO\TEMPO_DC_FPA_TEMP_SENSTIVITY_NEW'
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')
    linearity = 1 # option to turn on/off linearity correction
    smear = 1
    cross_talk = 1
    temp_correction = 1

    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)
    
    temperature_files = [each for each in os.listdir(data_dir) \
                        if each.endswith('_PT_Dark')]
    for k in range(1, len(temperature_files)):
        
        image_data_files = os.path.join(data_dir, temperature_files[k],
                                        'Script_Data', 'saved_quads')
        
        data_path_all = [each for each in os.listdir(image_data_files)
                         if not each.endswith('_118000.dat.sav')
                         if each.endswith('.dat.sav')
                         ]
        
        data_path_op_int = [each for each in os.listdir(image_data_files) \
                         if each.endswith('118000.dat.sav')][-1]
        
        data_path_all.append(data_path_op_int)
        save_dir = os.path.join(image_save_dir, temperature_files[k])
        all_quad_D_odd = []
        all_quad_D_even = []
        all_quad_C_odd = []
        all_quad_C_even = []
        all_quad_B_odd = []
        all_quad_B_even = []  
        all_quad_A_odd = [] 
        all_quad_A_even = []
        all_int_time = []
        
      
        for data_path in data_path_all:
            data_path_name_split = data_path.split('_')
            int_time = round(int(data_path_name_split[-1].split('.')[0]))
            int_time = int(int_time)/1000
            #print(int_time)
            all_int_time.append(int_time)
                          
            data_file = os.path.join(image_data_files, data_path)
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
                print('Integ. Time =' , int_time)
                print(coadds)
                
            else:
                #print('Telemetry file missing')
                coadds = 100
                #int_time = 6000
                fpe_temp = 25
                fpa_temp = 25
            
           ##############Image Processing begins here##############################
        
            full_frame = read_idl_file(data_file)
            
            if np.array(full_frame).ndim == 4:
                full_frame = np.mean(full_frame, axis=0)
            else:            
                full_frame = perform_coaddition_correction(full_frame, coadds)
            
            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
    
            ##########################OFFSET REMOVAL###############################
            # Input : Full_Frame TEMPO IMAGE
            # Otput : Bias Subtracted Active Region. The active region dimensions are
            #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
            #to store storage summation information
            bias_removed_quads = perform_bias_subtraction(full_frame)         
            #---------------------------------------------------------------------
            #########################NON LINEARITY CORRECTION#######################
            if temp_correction:
                temp_corrected_quads =   perform_temp_correction(bias_removed_quads, fpe_temp)
            else:
                temp_corrected_quads = bias_removed_quads            
    
            
           #########################NON LINEARITY CORRECTION#######################
            # Input : Bias Subtracted Active regions
            # Output : Linearized Active Region Quads
            # pass the difference instead of coadds
            if linearity:
                linearized_a = apply_linearity_correction(temp_corrected_quads [0, :, :],
                                                          quads[0])
                linearized_b = apply_linearity_correction(temp_corrected_quads [1, :, :],
                                                          quads[1])
                linearized_c = apply_linearity_correction(temp_corrected_quads [2, :, :],
                                                          quads[2])
                linearized_d = apply_linearity_correction(temp_corrected_quads[3, :, :],
                                                          quads[3])
                linearized_quads = np.array([linearized_a, linearized_b,
                                             linearized_c, linearized_d])
                #print(np.mean(linearized_quads))
            else:
                linearized_quads = temp_corrected_quads       
    
    
            if cross_talk:
                cross_talk_removed_quads = remove_cross_talk(linearized_quads)
            else:
                cross_talk_removed_quads = linearized_quads
                
            if smear:
                smear_removed_quads = perform_smear_removal(np.array(cross_talk_removed_quads),
                                                            int_time)
            else:
                smear_removed_quads = cross_talk_removed_quads
            
            processed_image = create_final_image(np.array(smear_removed_quads))
            quad_d = processed_image[0:1028, 0:1024]        
            quad_d_odd = quad_d[:, ::2]
            quad_d_even = quad_d[:, 1::2]
            all_quad_D_odd.append(np.mean(filter_outlier_median(quad_d_odd[300:900, 200:400])))
            all_quad_D_even.append(np.mean(filter_outlier_median(quad_d_even[300:900, 200:400])))
        
            quad_c = processed_image[0:1028, 1024:]
            quad_c_odd = quad_c[:, ::2]
            quad_c_even = quad_c[:, 1::2]
            all_quad_C_odd.append(np.mean(filter_outlier_median(quad_c_odd[300:900, 200:400])))
            all_quad_C_even.append(np.mean(filter_outlier_median(quad_c_even[300:900, 200:400])))
        
            quad_a = processed_image[1028:, 0:1024]
            quad_a_odd = quad_a[:, ::2]
            
            quad_a_even = quad_a[:, 1::2]
            all_quad_A_odd.append(np.mean(filter_outlier_median(quad_a_odd[300:900, 200:400])))
            all_quad_A_even.append(np.mean(filter_outlier_median(quad_a_even[300:900, 200:400])))
        
            quad_b = processed_image[1028:, 1024:]
            quad_b_odd = quad_b[:, ::2]
            quad_b_even = quad_b[:, 1::2]
            all_quad_B_odd.append(np.mean(filter_outlier_median(quad_b_odd[300:900, 200:400])))
            all_quad_B_even.append(np.mean(filter_outlier_median(quad_b_even[300:900, 200:400])))
            
        dframe1 = pd.DataFrame({'Int_time.' : all_int_time,
                                'Avg_Quad_A_odd' : all_quad_A_odd,
                                'Avg_Quad_A_even' : all_quad_A_even,
                                'Avg_Quad_B_odd' : all_quad_B_odd,
                                'Avg_Quad_B_even' : all_quad_B_even,
                                'Avg_Quad_C_odd' : all_quad_C_odd,
                                'Avg_Quad_C_even' : all_quad_C_even,
                                'Avg_Quad_D_odd' : all_quad_D_odd,
                                'Avg_Quad_D_even' : all_quad_D_even                                
                               })
        processed_file_dir = os.path.join(save_dir, temperature_files[k])
        if not os.path.exists(processed_file_dir):
            os.makedirs(processed_file_dir)
        
        processed_file_name = processed_file_dir + '/'+ temperature_files[k]+ '_Photon_transfer_data_all_FPA.csv'
        dframe1.to_csv(processed_file_name)
            
       
    
    #np.savetxt(variance_file_name, var_all,fmt='%1.3f',  delimiter=",")
    #np.savetxt(mean_file_name, signal_all,fmt='%1.3f', delimiter=",")
    

if __name__ == "__main__":
    main()
