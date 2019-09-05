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
#import h5py
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
     """Calculate the mean of all processed dark data
    """

     data = ([each for each in os.listdir(data_dir)
             if each.endswith('.csv')])
     all_data = []
     for num_data in data:
        #print(num_data)
        processed_data = os.path.join(data_dir, num_data)
        data = np.genfromtxt(processed_data, delimiter=',')       
        all_data.append(data)
        #print

     all_data = np.array(all_data)    
     return np.mean(all_data, axis=0)



def calculate_mean(data_dir, dark_data):
    """Calculate the mean of all processed_data
    """

    data = ([each for each in os.listdir(data_dir)
             if each.endswith('.csv')])
    all_data = []
    for num_data in data:
        #print(num_data)
        processed_data = os.path.join(data_dir, num_data)
        data = np.genfromtxt(processed_data, delimiter=',')
        data = data-dark_data
        all_data.append(data)
        #print

    all_data = np.array(all_data)    
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
    diffuser_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Diffuser_Irradiance\saved_quads\970_ms\saved_plots'
    diffuser_dark_data = os.path.join(diffuser_data_dir, 'Dark_data')
    
    mean_diffuser_dark_data = calculate_mean_dark(diffuser_dark_data)
    diffuser_dc_corrected = calculate_mean(diffuser_data_dir, mean_diffuser_dark_data)

     #Let's correct for dark current and work in signal rates unit
    
    diffuser_dc_corrected_mean = np.round(np.mean(diffuser_dc_corrected, axis=0), 2)
    diffuser_dc_corrected_var = np.round(np.var(diffuser_dc_corrected, axis=0), 2)
    


    mean_save_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Signal_To_Noise\Diffuser_970ms'
    mean_save_irradiance = os.path.join(mean_save_dir, 'mean_signal_970ms.csv')
    var_save_irradiance = os.path.join(mean_save_dir, 'var_signal_970ms.csv')
    np.savetxt(mean_save_irradiance, diffuser_dc_corrected_mean[:, 6:2041], fmt='%1.3f', delimiter=",")
    np.savetxt(var_save_irradiance, diffuser_dc_corrected_var[:, 6:2041], fmt='%1.3f', delimiter=",")
  






if __name__ == "__main__":
    main()
