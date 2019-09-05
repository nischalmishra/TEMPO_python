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
from scipy.ndimage import uniform_filter
#from random import randint
from scipy.optimize import curve_fit
#from scipy.interpolate import interp1d
#pylint: disable= E1101

from Process_TEMPO_Spectrometer_Data import read_idl_file, parse_telemetry_file

import matplotlib.pyplot as plt



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


def extract_trailing_overclocks(raw_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads_odd = []
    all_quads_even = []
    #num_quads = raw_quads.shape[0]
    num_quads=2
    for quads in range(1, num_quads):        
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]        
        even_detector_bias = trailing_overclocks[ :, ::2]
        #even_detector_bias = trailing_overclocks
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        #even_detector_bias = even_detector_bias[:, 4:]  
#        plt.imshow(even_detector_bias)
#        plt.show()
#        cc
        odd_detector_bias = trailing_overclocks[:, 1::2]
        #odd_samples = odd_detector_bias[:, 4:]
        #rows, cols = odd_samples.shape
        #odd_detector_bias = odd_samples[:, 0:cols-1]
        all_quads_even.append(even_detector_bias)
        all_quads_odd.append(odd_detector_bias)
        
    return np.array(even_detector_bias), np.array(odd_detector_bias)





def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Dark_data\2017.07.10'
    telemetry_dir = r'E:\Image Data\2018.07.10'
    image_dir = os.path.join(data_dir, 'saved_quads')
    save_dir = os.path.join(image_dir, 'processed_image')
    image_save_dir = os.path.join(image_dir, 'saved_plots')
    telemetry_dir = os.path.join(telemetry_dir, 'processed/h5')    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
#    data_path_all = sorted([each for each in os.listdir(image_dir)
#                            if each.startswith('2017_07_30_00_24_59_38064')
#                            and  each.endswith('.sav')])
#
    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.sav')])
    
    all_quads_odd = []
    all_quads_even = []
    print('Total data = ', len(data_path_all))

    for data_path in data_path_all:
        data_file = os.path.join(image_dir, data_path)
        print(data_path)
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
            print('Integ. Time =' , int_time)
            
        else:
            print('Telemetry file missing')
            coadds = 10
            int_time = 6000
            fpe_temp = 43
            fpa_temp = -21.5
        
       ##############Image Processing begins here##############################
        full_frame = read_idl_file(data_file)
        full_frame = full_frame/coadds
        num_quads = full_frame.shape[0]
        #print(num_quads)
        
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


        even_quad, odd_quad = extract_trailing_overclocks(full_frame)
        all_quads_even.append(even_quad)
        all_quads_odd.append(odd_quad)

        
        
    all_quads_odd = np.array(all_quads_odd)
    all_quads_even = np.array(all_quads_even)
    x,y,z = np.shape(all_quads_odd)
    all_quads_odd = np.reshape(all_quads_odd,(x*y*z, 1))
    
    
    x,y,z = np.shape(all_quads_even)
    all_quads_even = np.reshape(all_quads_even,(x*y*z, 1))
    
    
    
    #variance_all_odd = np.var(all_quads_odd, axis=0)
    #rows, cols = variance_all_odd.shape
    
    
    
    #variance_all = np.reshape(variance_all_odd, (rows*cols,1))
    #mean_variance = np.mean(variance_all)
    #text = 'Mean Var = '+ str(round(mean_variance, 2)) 
   # text1 = '\nRN =' + str(round(np.sqrt(mean_variance), 2)) + ' DN'
    #label1 = text+text1
    
    plt.figure()
    plt.hist(all_quads_odd, 30, facecolor='blue', alpha=0.5)
    #plt.ylim(0, 380)
    #plt.xlim(0, 9)
    plt.grid(linestyle=':')
    plt.xlabel('DN')
    plt.ylabel('Frequency')
    plt.title('Histogram of Trailing Serial Overclocks (Quad B Odd)')
    #plt.legend(loc='best')
    #plt.show()
   
    plt.savefig(r'C:\Users\nmishra\Desktop\RN\2017_07_10_FPE_39.69_FPA_-21.18\Overclocks\Outlier_included' +'/'+'QuadB_Odd.png')
    plt.close('all')    
    #variance_all_even = np.var(all_quads_odd, axis=0)
    #rows, cols = variance_all_even.shape
    
    #variance_all = np.reshape(variance_all_even, (rows*cols,1))
    #mean_variance = np.mean(variance_all)
    #text = 'Mean Var = '+ str(round(mean_variance, 2)) 
    #text1 = '\nRN =' + str(round(np.sqrt(mean_variance), 2)) + ' DN'
    #label1 = text+text1
    
    plt.figure()
    plt.hist(all_quads_even, 30, facecolor='red', alpha=0.5)
    #plt.ylim(0, 380)
    #plt.xlim(0, 9)
    plt.grid(linestyle=':')
    plt.xlabel('DN')
    plt.ylabel('Frequency')
    plt.title('Histogram of Trailing Serial Overclocks (Quad B Even)')
    #plt.legend(loc='best')
    plt.savefig(r'C:\Users\nmishra\Desktop\RN\2017_07_10_FPE_39.69_FPA_-21.18\Overclocks\Outlier_included' +'/'+'QuadB_Even.png')
    


if __name__ == "__main__":
    main()
