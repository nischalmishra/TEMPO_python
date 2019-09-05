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
import h5py
from Process_TEMPO_Spectrometer_Data import read_idl_file

import matplotlib.pyplot as plt

def parse_telemetry_file(file_path, telemetry_file_name):
    """ Read the variables required from telemetry files to process the spectro
    meter image. Values required are integration time, coadds and temperatures
    """
    h5file_name = telemetry_file_name + '.h5'
    hdf_name = os.path.join(file_path, h5file_name)
    #print(h5file_name)
    with h5py.File(hdf_name, 'r') as hdf:
        group1 = hdf.get('telemetry')
        group2 = group1.get('calculated_values')
        coadds = group2.attrs['ncoadds']
        int_time = group2.attrs['tint']
        ccd_temp = group2.attrs['ccd_temp']
        fpe_temp = group2.attrs['fpe_temp']
        # sanity check
        #print(coadds, int_time, ccd_temp, fpe_temp)

    return coadds, int_time, fpe_temp, ccd_temp

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


def plot_tsoc(raw_quads, int_time):
    quad = ['Quad A','Quad B', 'Quad C','Quad D']
    
    num_quads = raw_quads.shape[0]
    for quads in range(3, num_quads):
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]
        odd_detector_bias = trailing_overclocks[:, ::2]
        odd_detector_bias = odd_detector_bias[:, 4:]
        even_detector_bias = trailing_overclocks[:, 1::2]
        even_detector_bias = even_detector_bias[:, 4:-1]
        labels =["Sixth","Seventh","Eighth" ,"Ninth", "Tenth"]
        fig = plt.figure()
        ax1 = fig.add_subplot(211)        
        ax1.set_title('Trailing Serial Overclock Profile of '+ quad[quads]  + ', Integ, Time = '+ str(round(int_time, 0))+' ms')
        ax1.set_ylabel('TSOC, Odd Lines (DN)')
        ax1.plot(odd_detector_bias)       
        ax1.grid(True, linestyle=':', )
        #ax1.legend(loc='best')
        
        ax2 = fig.add_subplot(212)
        ax2.plot(even_detector_bias)       
        ax2.grid(True, linestyle=':', )
        ax2.set_ylabel('TSOC, Even Lines (DN)')
        ax2.set_xlabel('Pixel #')
        #ax2.legend(loc='best')
         
        plt.show()
        cc
        
        

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
    num_quads = raw_quads.shape[0]
    
    #num_quads=2
    for quads in range(0, num_quads):        
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]        
        odd_detector_bias = trailing_overclocks[ :, ::2]
        
                
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        odd_detector_bias = odd_detector_bias[:, 4:-1]  
        
        even_detector_bias = trailing_overclocks[:, 1::2]
        even_detector_bias= even_detector_bias[:, 4:-1]        
        all_quads_even.append(np.mean(even_detector_bias, axis=1))        
        all_quads_odd.append(np.mean(odd_detector_bias, axis=1))
        plt.plot(np.mean(even_detector_bias, axis=1).T,'k.')
        plt.grid(True, linestyle=':')
        plt.title('Trailing Overclock Plot, Quad D, int.time = 8 ms')
        plt.xlabel('----> Read out sequence')
        plt.ylabel('TSOC (DN)')
        
        plt.show()
        
    
    return np.array(all_quads_even), np.array(all_quads_odd)





def extract_active_region(raw_quads):
     # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads_odd = []
    all_quads_even = []
    num_quads = raw_quads.shape[0]
    #saturation = 14744
    #num_quads=2
    
    for quads in range(0, num_quads):        
        active_quad = raw_quads[quads, 2:1030, 7:1024].astype('float')
       
        active_quad[active_quad >= 0.90*16383] = np.NAN          
        active_quad = active_quad.astype(np.float)
        odd_detector_active = (active_quad[ :, ::2]) 
        #odd_detector_active[odd_detector_active>=saturation] = np.nan
        even_detector_active = active_quad[:, 1::2] 
       # even_detector_active[even_detector_active>=saturation] = np.nan
        all_quads_even.append(np.nanmean(even_detector_active, axis=1))        
        all_quads_odd.append(np.nanmean(odd_detector_active, axis=1))
        
    return np.array(all_quads_even), np.array(all_quads_odd)
        
    


def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC'
    telemetry_dir = os.path.join(r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.28\processed\h5') 
    image_dir = os.path.join(data_dir, 'saved_quads')
   
    image_save_dir = os.path.join(image_dir, 'saved_plots', 'SDO/updated_05_30')
       
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)

    data_path_all = sorted([each for each in os.listdir(image_dir)
                            if each.endswith('.sav')])
    
    all_active_A_odd = []
    all_active_B_odd = []
    all_active_C_odd = []
    all_active_D_odd = []
    
    all_active_A_even = []
    all_active_B_even = []
    all_active_C_even = []
    all_active_D_even = []   
    
    all_tsoc_A_odd = []
    all_tsoc_B_odd = []
    all_tsoc_C_odd = []
    all_tsoc_D_odd = []
    
    all_tsoc_A_even = []
    all_tsoc_B_even = []
    all_tsoc_C_even = []
    all_tsoc_D_even = []
    all_int_time = [ ]
    
    print('Total data = ', len(data_path_all))

    for data_path in range(0, len(data_path_all), 100):
        data_path = data_path_all[data_path]
        #print(data_path)
        data_file = os.path.join(image_dir, data_path)
       #################### Telemetry Statistics###############################
        # Now lets parse the telemetry file to get the required information###
        #Note; The function can return the CCD temp and FPE temp. But these are
        #not needed on the image processing algorithm.. atleast for the time being
        telemetry_file_name = data_path.split('.')[0]
        telemetry_file = os.path.join(telemetry_dir, telemetry_file_name+'.h5')
        

        if os.path.exists(telemetry_file):
            coadds, int_time, fpe_temp, fpa_temp = parse_telemetry_file(telemetry_dir,
                                                                        telemetry_file_name)
#            print('FPE Temp. = ', round(fpe_temp, 2), 'FPA Temp = ',
#                  round(fpa_temp, 2))
            print('Integ. Time =' , int_time)
            #if float(int_time)>=40:
               # break
            
        else:
            continue
            print('Telemetry file missing')
#            coadds = 10
#            int_time = 6000
#            fpe_temp = 43
#            fpa_temp = -21.5
        
       ##############Image Processing begins here##############################
        full_frame = read_idl_file(data_file)
        #full_frame = full_frame/coadds
        num_quads = full_frame.shape[0]
        #print(num_quads)
        
        #check_image = create_image_active_region(full_frame)
        #raw_image = create_final_image(np.array(full_frame))
       # print('Max. Val. Raw = ', np.max(raw_image))
        quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
        #plot_tsoc(full_frame, int_time)
        all_int_time.append(int_time)
        ##########################OFFSET REMOVAL###############################
        # Input : Full_Frame TEMPO IMAGE
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
        #to store storage summation information
        # For dark data, only offset removal is needed to compute the dark current
        # For light data, additional processing such as linearity, smear and
        #cross talk is needed

        even_quad_tsoc, odd_quad_tsoc = extract_trailing_overclocks(full_frame)
        all_tsoc_A_odd.append(odd_quad_tsoc[0, :])        
        all_tsoc_B_odd.append(odd_quad_tsoc[1, :])
        all_tsoc_C_odd.append(odd_quad_tsoc[2, :])
        all_tsoc_D_odd.append(odd_quad_tsoc[3, :])
        
        all_tsoc_A_even.append(even_quad_tsoc[0, :])        
        all_tsoc_B_even.append(even_quad_tsoc[1, :])
        all_tsoc_C_even.append(even_quad_tsoc[2, :])
        all_tsoc_D_even.append(even_quad_tsoc[3, :])
        
    
        
        even_quad_active, odd_quad_active = extract_active_region(full_frame)
        


        all_active_A_odd.append(odd_quad_active[0, :])        
        all_active_B_odd.append(odd_quad_active[1, :])
        all_active_C_odd.append(odd_quad_active[2, :])
        all_active_D_odd.append(odd_quad_active[3, :])
        
        all_active_A_even.append(even_quad_active[0, :])        
        all_active_B_even.append(even_quad_active[1, :])
        all_active_C_even.append(even_quad_active[2, :])
        all_active_D_even.append(even_quad_active[3, :])
        
       

        
    all_tsoc_A_odd = np.array(all_tsoc_A_odd)     
    all_tsoc_B_odd = np.array(all_tsoc_B_odd) 
    all_tsoc_C_odd= np.array(all_tsoc_C_odd)
    all_tsoc_D_odd= np.array(all_tsoc_D_odd)
        
    all_tsoc_A_even = np.array(all_tsoc_A_even) 
    all_tsoc_B_even = np.array(all_tsoc_B_even) 
    all_tsoc_C_even= np.array(all_tsoc_C_even)
    all_tsoc_D_even= np.array(all_tsoc_D_even)
    
    all_active_A_odd = np.array(all_active_A_odd)   
    all_active_B_odd = np.array(all_active_B_odd) 
    all_active_C_odd= np.array(all_active_C_odd)
    all_active_D_odd= np.array(all_active_D_odd)
        
    all_active_A_even = np.array(all_active_A_even) 
    all_active_B_even = np.array(all_active_B_even) 
    all_active_C_even= np.array(all_active_C_even)
    all_active_D_even= np.array(all_active_D_even)
    
  
    
    
    
    all_int_time = np.array(all_int_time)
    
    int_time_save = os.path.join(image_save_dir,'integration_time.csv')
    quad_A_active_odd = os.path.join(image_save_dir,' quad_A_active_odd.csv')
    quad_B_active_odd = os.path.join(image_save_dir,' quad_B_active_odd.csv')
    quad_C_active_odd = os.path.join(image_save_dir,' quad_C_active_odd.csv')
    quad_D_active_odd = os.path.join(image_save_dir,' quad_D_active_odd.csv')
    
    quad_A_active_even = os.path.join(image_save_dir,' quad_A_active_even.csv')
    quad_B_active_even = os.path.join(image_save_dir,' quad_B_active_even.csv')
    quad_C_active_even = os.path.join(image_save_dir,' quad_C_active_even.csv')
    quad_D_active_even = os.path.join(image_save_dir,' quad_D_active_even.csv')
    
    
    
    
    
    quad_A_tsoc_odd = os.path.join(image_save_dir,' quad_A_tsoc_odd.csv')
    quad_B_tsoc_odd = os.path.join(image_save_dir,' quad_B_tsoc_odd.csv')
    quad_C_tsoc_odd = os.path.join(image_save_dir,' quad_C_tsoc_odd.csv')
    quad_D_tsoc_odd = os.path.join(image_save_dir,' quad_D_tsoc_odd.csv')
    
    quad_A_tsoc_even = os.path.join(image_save_dir,' quad_A_tsoc_even.csv')
    quad_B_tsoc_even = os.path.join(image_save_dir,' quad_B_tsoc_even.csv')
    quad_C_tsoc_even = os.path.join(image_save_dir,' quad_C_tsoc_even.csv')
    quad_D_tsoc_even = os.path.join(image_save_dir,' quad_D_tsoc_even.csv')
    
    
    
    np.savetxt(quad_A_active_odd, all_active_A_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_B_active_odd, all_active_B_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_C_active_odd, all_active_C_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_D_active_odd, all_active_D_odd, fmt='%1.3f',  delimiter=",")
    
    np.savetxt(quad_A_active_even, all_active_A_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_B_active_even, all_active_B_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_C_active_even, all_active_C_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_D_active_even, all_active_D_even, fmt='%1.3f',  delimiter=",")
   
    
    np.savetxt(quad_A_tsoc_odd, all_tsoc_A_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_B_tsoc_odd, all_tsoc_B_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_C_tsoc_odd, all_tsoc_C_odd, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_D_tsoc_odd, all_tsoc_D_odd, fmt='%1.3f',  delimiter=",")
    
    np.savetxt(quad_A_tsoc_even, all_tsoc_A_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_B_tsoc_even, all_tsoc_B_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_C_tsoc_even, all_tsoc_C_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(quad_D_tsoc_even, all_tsoc_D_even, fmt='%1.3f',  delimiter=",")
    np.savetxt(int_time_save, all_int_time, fmt='%1.3f',  delimiter="," )
   
    
    
    


if __name__ == "__main__":
    main()
