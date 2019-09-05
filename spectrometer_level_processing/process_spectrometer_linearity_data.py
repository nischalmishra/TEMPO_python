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
#from random import randint

#from scipy.interpolate import interp1d
#pylint: disable= E1101

from Process_TEMPO_Spectrometer_Data import read_idl_file, perform_bias_subtraction,\
                                     parse_telemetry_file


def extract_active_region(bias_subtracted_quads):
     # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads_odd = []
    all_quads_even = []
    num_quads = bias_subtracted_quads.shape[0]
    #print(bias_subtracted_quads.shape)

    #num_quads=2
    for quads in range(0, num_quads):
        active_quad = bias_subtracted_quads[quads, :, 7:1024].astype('float')
        #active_quad[active_quad >= 0.90*16383] = np.NAN
        odd_detector_active = active_quad[ :, ::2]
        even_detector_active = active_quad[:, 1::2]
        
        all_quads_even.append(np.nanmean(even_detector_active))
        all_quads_odd.append(np.nanmean(odd_detector_active))
    return all_quads_even, all_quads_odd



def extract_uncertainty(bias_subtracted_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads_odd_var = []
    all_quads_even_var = []
    num_quads = bias_subtracted_quads.shape[0]
    #saturation = 14744
    #num_quads=2
    for quads in range(0, num_quads):
        active_quad = bias_subtracted_quads[quads :, 7:1024].astype('float') 
        active_quad[active_quad >= 0.90*16383] = np.NAN          
        active_quad = active_quad.astype(np.float)
        odd_detector_active = (active_quad[ :, ::2]) 
        #odd_detector_active[odd_detector_active>=saturation] = np.nan
        even_detector_active = active_quad[:, 1::2] 
       # even_detector_active[even_detector_active>=saturation] = np.nan
        all_quads_even_var.append(np.nanvar(even_detector_active))        
        all_quads_odd_var.append(np.nanvar(odd_detector_active))
        
    return np.array(all_quads_even_var), np.array(all_quads_odd_var)


def main():
    """
    This is the main function that does all the processing. It calls the required
    functions to process the raw image. Raw images are saved in .sav (IDL)
    format and read through this script.
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC'
    telemetry_dir = os.path.join(r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.28\processed\h5')
    image_dir = os.path.join(data_dir, 'saved_quads')
    image_save_dir = os.path.join(image_dir, 'saved_plots', 'Linearity','Mean_Stats\Variance_inlcuded')
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
    
    all_active_A_odd_var = []
    all_active_B_odd_var = []
    all_active_C_odd_var = []
    all_active_D_odd_var = []
    
    all_active_A_even_var = []
    all_active_B_even_var = []
    all_active_C_even_var = []
    all_active_D_even_var = []
    
    
    int_time_all = []

    print('Total data = ', len(data_path_all))
    for data_path in range(0, len(data_path_all)):
        #print(data_path)
        data_path = data_path_all[data_path]
        print(data_path)

        data_file = os.path.join(image_dir, data_path)
       
       ################## Telemetry Statistics###############################
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
        
        int_time_all.append(int_time)
        
       #############Image Processing begins here##############################
        full_frame = read_idl_file(data_file)
        #full_frame = full_frame/coadds
        #print(num_quads)

        #check_image = create_image_active_region(full_frame)
        #raw_image = create_final_image(np.array(full_frame))
       # print('Max. Val. Raw = ', np.max(raw_image))


        ##########################OFFSET REMOVAL###############################
        # Input : Full_Frame TEMPO IMAGE
        # Otput : Bias Subtracted Active Region. The active region dimensions are
        #now 1028*1024. 2 pixel lines from SMEAR overclocks are now used for
        #to store storage summation information
        # For dark data, only offset removal is needed to compute the dark current
        # For light data, additional processing such as linearity, smear and
        #cross talk is needed
        bias_removed_quads = perform_bias_subtraction(full_frame)
        even_quad_active, odd_quad_active = extract_active_region(bias_removed_quads)

        all_active_A_odd.append(odd_quad_active[0])
        all_active_B_odd.append(odd_quad_active[1])
        all_active_C_odd.append(odd_quad_active[2])
        all_active_D_odd.append(odd_quad_active[3])

        all_active_A_even.append(even_quad_active[0])
        all_active_B_even.append(even_quad_active[1])
        all_active_C_even.append(even_quad_active[2])
        all_active_D_even.append(even_quad_active[3])
        
        
        even_quad_var, odd_quad_var = extract_uncertainty(bias_removed_quads)
       
        all_active_A_odd_var.append(odd_quad_var[0])        
        all_active_B_odd_var.append(odd_quad_var[1])
        all_active_C_odd_var.append(odd_quad_var[2])
        all_active_D_odd_var.append(odd_quad_var[3])
        
        all_active_A_even_var.append(even_quad_var[0])        
        all_active_B_even_var.append(even_quad_var[1])
        all_active_C_even_var.append(even_quad_var[2])
        all_active_D_even_var.append(even_quad_var[3])

    all_active_A_odd = np.array(all_active_A_odd)
    all_active_B_odd = np.array(all_active_B_odd)
    all_active_C_odd = np.array(all_active_C_odd)
    all_active_D_odd = np.array(all_active_D_odd)

    all_active_A_even = np.array(all_active_A_even)
    all_active_B_even = np.array(all_active_B_even)
    all_active_C_even = np.array(all_active_C_even)
    all_active_D_even = np.array(all_active_D_even)
    
    all_active_A_odd_var = np.array(all_active_A_odd_var)   
    all_active_B_odd_var = np.array(all_active_B_odd_var) 
    all_active_C_odd_var= np.array(all_active_C_odd_var)
    all_active_D_odd_var= np.array(all_active_D_odd_var)
        
    all_active_A_even_var = np.array(all_active_A_even_var) 
    all_active_B_even_var = np.array(all_active_B_even_var) 
    all_active_C_even_var= np.array(all_active_C_even_var)
    all_active_D_even_var= np.array(all_active_D_even_var)
    
    

    quad_A_active_odd = os.path.join(image_save_dir, 'quad_A_active_odd.csv')
    quad_B_active_odd = os.path.join(image_save_dir, 'quad_B_active_odd.csv')
    quad_C_active_odd = os.path.join(image_save_dir, 'quad_C_active_odd.csv')
    quad_D_active_odd = os.path.join(image_save_dir, 'quad_D_active_odd.csv')

    quad_A_active_even = os.path.join(image_save_dir, 'quad_A_active_even.csv')
    quad_B_active_even = os.path.join(image_save_dir, 'quad_B_active_even.csv')
    quad_C_active_even = os.path.join(image_save_dir, 'quad_C_active_even.csv')
    quad_D_active_even = os.path.join(image_save_dir, 'quad_D_active_even.csv')
    
    
    quad_A_active_odd_var = os.path.join(image_save_dir, 'quad_A_active_odd_var.csv')
    quad_B_active_odd_var = os.path.join(image_save_dir, 'quad_B_active_odd_var.csv')
    quad_C_active_odd_var = os.path.join(image_save_dir, 'quad_C_active_odd_var.csv')
    quad_D_active_odd_var = os.path.join(image_save_dir, 'quad_D_active_odd_var.csv')

    quad_A_active_even_var = os.path.join(image_save_dir, 'quad_A_active_even_var.csv')
    quad_B_active_even_var = os.path.join(image_save_dir, 'quad_B_active_even_var.csv')
    quad_C_active_even_var = os.path.join(image_save_dir, 'quad_C_active_even_var.csv')
    quad_D_active_even_var = os.path.join(image_save_dir, 'quad_D_active_even_var.csv')

    integration_time = os.path.join(image_save_dir, 'int_time.csv')
    
    
    np.savetxt(quad_A_active_odd, all_active_A_odd, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_B_active_odd, all_active_B_odd, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_C_active_odd, all_active_C_odd, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_D_active_odd, all_active_D_odd, fmt='%1.3f', delimiter=",")

    np.savetxt(quad_A_active_even, all_active_A_even, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_B_active_even, all_active_B_even, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_C_active_even, all_active_C_even, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_D_active_even, all_active_D_even, fmt='%1.3f', delimiter=",")
    
    
    np.savetxt(quad_A_active_odd_var, all_active_A_odd_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_B_active_odd_var, all_active_B_odd_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_C_active_odd_var, all_active_C_odd_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_D_active_odd_var, all_active_D_odd_var, fmt='%1.3f', delimiter=",")

    np.savetxt(quad_A_active_even_var, all_active_A_even_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_B_active_even_var, all_active_B_even_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_C_active_even_var, all_active_C_even_var, fmt='%1.3f', delimiter=",")
    np.savetxt(quad_D_active_even_var, all_active_D_even_var, fmt='%1.3f', delimiter=",")
    
    
    
    np.savetxt(integration_time, int_time_all, fmt='%1.3f', delimiter=",")


if __name__ == "__main__":
    main()
