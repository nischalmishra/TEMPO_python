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
import h5py
#import pandas as pd

#from random import randint

#from scipy.interpolate import interp1d
#pylint: disable= E1101



import matplotlib.pyplot as plt


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
    image = fig_ax.imshow(image_data[0:1028, :], cmap='nipy_spectral',
                          origin='lower', interpolation='none')
    #image = fig_ax.imshow(np.array(image_data), cmap='nipy_spectral',
                          #origin='lower', interpolation='none')
    plt.title(title)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show()
    #plt.pause(0.1)
    #plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    plt.close('all')



def calculate_mean_dark(data_dir):
    """Calculate the mean of all dark processed_data
    """

    data = ([each for each in os.listdir(data_dir)
             if each.endswith('.h5')])
    
    all_data = []
    for num_data in data:
        #print(num_data)
        processed_data = os.path.join(data_dir, num_data)
        file = h5py.File(processed_data, 'r')        
        data = file.get('Processed_data')
        all_data.append(data)
        #print

    all_data = np.array(all_data)
    all_data = np.mean(all_data, axis=0)
    return all_data
    


def calculate_mean(data_dir):
    """ Calculate the mean of Light data with straylight correction.
    Stray light correction version 1 is to use shadow pixels ( first 3 and last 5)
    """
    data = ([each for each in os.listdir(data_dir)
             if each.endswith('.h5')])
    all_data = []
    for num_data in data:
        processed_data = os.path.join(data_dir, num_data)
        file = h5py.File(processed_data, 'r')        
        data = file.get('Processed_data')       
        all_data.append(data)
    all_data = np.array(all_data)
    all_data = np.mean(all_data, axis=0)
    return all_data
         


def find_nearest(array, value):
    """ Find the nearest neighbor of the given list element
    """
    idx = (np.abs(array-value)).argmin()
    return array[idx]




def main():
    """
    This is the main function that does all the processing. It does all the analysis portion

    """
    #diffuser_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Diffuser_Irradiance\saved_quads\3845_ms\saved_plots_modified'
    #diffuser_light_data = os.path.join(diffuser_data_dir,'Light_data')
    #diffuser_dark_data = os.path.join(diffuser_data_dir, 'Dark_data')
    #print(diffuser_data_dir)
    #cc
    #int_time_diffuser = 3845.0

#    radiance_data_dir_UV = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_UV_Lamp\saved_quads\saved_plots_modified'
#    radiance_light_data_UV = os.path.join(radiance_data_dir_UV,'DSS-Y')
#    radiance_dark_data_UV = os.path.join(radiance_data_dir_UV, 'Dark_data')
##    #int_time_radiance = 93.0
##    print(radiance_data_dir_UV)
##
    radiance_data_dir_VIS = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\processed_h5'
    radiance_light_data_VISminusY = os.path.join(radiance_data_dir_VIS,'DSS-Y')
    radiance_light_data_VIS_Center = os.path.join(radiance_data_dir_VIS,'DSS_Center')
    radiance_light_data_VISplusY = os.path.join(radiance_data_dir_VIS,'DSS+Y')    
    radiance_dark_data_VIS = os.path.join(radiance_data_dir_VIS, 'Dark_Data')



   # mean_diffuser_data = calculate_mean(diffuser_light_data)
    #mean_diffuser_dark_data = calculate_mean(diffuser_dark_data)
#
#
#     #Let's correct for dark current and work in signal rates unit
#    diffuser_dc_corrected = (mean_diffuser_data - mean_diffuser_dark_data)
#    diffuser_dc_corrected = np.round(diffuser_dc_corrected, 2)
#    diffuser_dc_corrected[diffuser_dc_corrected <0] = 0
#    mean_save_dir = os.path.join(diffuser_data_dir,'processed_average_data')
#    mean_save_irradiance = os.path.join(mean_save_dir, 'mean_irradiance_3845ms.csv')
#    np.savetxt(mean_save_irradiance, diffuser_dc_corrected, fmt='%1.3f', delimiter=",")
   



#    mean_radiance_data_UV = calculate_mean(radiance_light_data_UV)
#    mean_radiance_dark_data_UV = calculate_mean(radiance_dark_data_UV)
#    radiance_dc_corrected_UV = (mean_radiance_data_UV - mean_radiance_dark_data_UV)
#    radiance_dc_corrected_UV = np.round(radiance_dc_corrected_UV, 2)
#    radiance_dc_corrected_UV[radiance_dc_corrected_UV < 0] = 0
###
    mean_radiance_dark_data_VIS = calculate_mean_dark(radiance_dark_data_VIS)
    
    
    # Correct for Dark current
    mean_radiance_data_VISminusY = calculate_mean(radiance_light_data_VISminusY) - mean_radiance_dark_data_VIS
    mean_radiance_data_VIS_Center = calculate_mean(radiance_light_data_VIS_Center)- mean_radiance_dark_data_VIS
    mean_radiance_data_VISplusY = calculate_mean(radiance_light_data_VISplusY) - mean_radiance_dark_data_VIS
 
#    mean_save_dir_UV = os.path.join(radiance_data_dir_UV,'processed_average_data')
    mean_save_dir_VIS = os.path.join(radiance_data_dir_VIS,'Mean_Processed_Data')
    if not os.path.exists(mean_save_dir_VIS):
        os.makedirs(mean_save_dir_VIS)
    #
#    mean_save_radiance_UV = os.path.join(mean_save_dir_UV, 'mean_radiance_DSSminus_UV.csv')
#    mean_save_radiance_VIS = os.path.join(mean_save_dir_VIS, 'mean_radiance_DSSminus_VIS.csv')
#    #
#   
#    np.savetxt(mean_save_radiance_UV, radiance_dc_corrected_UV, fmt='%1.3f', delimiter=",")
#    np.savetxt(mean_save_radiance_VIS, radiance_dc_corrected_VIS, fmt='%1.3f', delimiter=",")
   

    #Write into h5file
    hf_name_center = os.path.join(mean_save_dir_VIS,'DSS_Center.h5')
    hf1 = h5py.File(hf_name_center,'w')
    hf1.create_dataset('Processed_data', data=mean_radiance_data_VIS_Center)
    hf1.close()

    hf_name_plus_Y = os.path.join(mean_save_dir_VIS,'DSS_plusY.h5')
    hf2 = h5py.File(hf_name_plus_Y,'w')
    hf2.create_dataset('Processed_data', data=mean_radiance_data_VISplusY )
    hf2.close()

    hf_name_minus_Y = os.path.join(mean_save_dir_VIS,'DSS_minusY.h5')
    hf3 = h5py.File( hf_name_minus_Y ,'w')
    hf3.create_dataset('Processed_data', data=mean_radiance_data_VISminusY)
    hf3.close()
    print('DONE')

if __name__ == "__main__":
    main()
