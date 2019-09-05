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







def main():
    """
    This is the main function that does all the processing. It does all the analysis portion

    """
    dark_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_UV_Lamp\saved_quads\saved_plots_modified\Dark_data'
    data = ([each for each in os.listdir(dark_data_dir)
             if each.endswith('.csv')])
    
    all_quad_A_even = []
    all_quad_A_odd = []
    all_quad_B_even = []
    all_quad_B_odd = []
    all_quad_C_even = []
    all_quad_C_odd = []
    all_quad_D_even = []
    all_quad_D_odd = []
    dframe1 = pd.DataFrame()
    for num_data in data[1:]:
        #print(num_data)
        processed_data = os.path.join(dark_data_dir, num_data)
        data = np.genfromtxt(processed_data, delimiter=',')
        quad_d = data[0:1028, 0:1024]        
        quad_d_odd = quad_d[:,::2]
        quad_d_even = quad_d[:, 1::2]
        all_quad_D_odd.append(np.mean(filter_outlier_median(quad_d_odd[300:900, 200:400])))
        all_quad_D_even.append(np.mean(filter_outlier_median(quad_d_even[300:900, 200:400])))
        
        quad_c = data[0:1028, 1024:]
        quad_c_odd = quad_c[:, ::2]
        quad_c_even = quad_c[:, 1::2]
        all_quad_C_odd.append(np.mean(filter_outlier_median(quad_c_odd[300:900, 200:400] )))
        all_quad_C_even.append(np.mean(filter_outlier_median(quad_c_even[300:900, 200:400] )))
        
        quad_a = data[1028:, 0:1024]
        quad_a_odd = quad_a[:, ::2]
        quad_a_even = quad_a[:, 1::2]
        all_quad_A_odd.append(np.mean(filter_outlier_median(quad_a_odd[300:900, 200:400] )))
        all_quad_A_even.append(np.mean(filter_outlier_median(quad_a_even[300:900, 200:400] )))
        
        quad_b = data[1028:, 1024:]
        quad_b_odd = quad_b[:, ::2]
        quad_b_even = quad_b[:, 1::2]
        all_quad_B_odd.append(np.mean(filter_outlier_median(quad_b_odd[300:900, 200:400] )))
        all_quad_B_even.append(np.mean(filter_outlier_median(quad_b_even[300:900, 200:400] )))
        
    dframe1 = pd.DataFrame({'Avg_Quad_A_odd' : all_quad_A_odd,
                                'Avg_Quad_A_even' : all_quad_A_even,
                                'Avg_Quad_B_odd' : all_quad_B_odd,
                                'Avg_Quad_B_even' : all_quad_B_even,
                                'Avg_Quad_C_odd' : all_quad_C_odd,
                                'Avg_Quad_C_even' : all_quad_C_even,
                                'Avg_Quad_D_odd' : all_quad_D_odd,
                                'Avg_Quad_D_even' : all_quad_D_even
                             })   
        
    dframe1.to_csv(r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Radiance_Cal_VIS_Lamp\saved_quads\Dark_current_data\processed_dark_data'+ '/'+'dark_current_all.csv')


if __name__ == "__main__":
    main()
