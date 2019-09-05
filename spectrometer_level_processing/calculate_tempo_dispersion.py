# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 10:36:14 2018

@author: nmishra
"""

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
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 3})
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
import matplotlib.ticker as mtick

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






def moving_average(data, n=10):
    ret = np.cumsum(data, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def gauss_function(x, a, x0, sigma):
    """ Fit Gaussian function
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))




def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def calculate_dispersion_all (save_dir, data_path_all, wavelen):
    

    dispersion_all = []
    dispersion_whole = []

    csfont = {'fontname':'Comic Sans MS'}
    hsfont = {'fontname':'Arial'}
    count = 0
    for i in wavelen:
        
        print(wavelen[count])
   
        count1=0
        for data in range(0, len(wavelen)-1):
            
            print(wavelen[count1+1])
            
            
        
            dframe_ref = pd.read_csv(save_dir+'/Iteration1/'+data_path_all[count], delimiter=',')
            dframe_comp = pd.read_csv(save_dir+'/Iteration1/'+data_path_all[count1+1], delimiter=',')
            #print(dframe_ref.head())
            #print(dframe_comp.head())
            #print(wavelen)
            delta_wavelen = float(wavelen[count]) -float( wavelen[count1+1])
            #print(delta_wavelen)


            delta_pixel = (dframe_ref['Peak Loc.'] -  dframe_comp['Peak Loc.']).values
            #print(delta_pixel)
            dispersion = np.true_divide(delta_wavelen, delta_pixel)
            dispersion_whole.append(dispersion)
            count1 = count1+1
        
        count=count+1
        dispersion_all.append(dispersion_whole)
   
    dispersion_all = np.array(dispersion_all)
    n, rows, cols = dispersion_all.shape
    #print(dispersion_all.shape)
    dispersion_all = np.reshape(dispersion_all, (n*rows*cols,1))
    unct = 100* np.std(dispersion_all)/np.mean(dispersion_all)
    text1 = 'Mean : ' + str(round(np.mean(dispersion_all), 4)) + \
           '\n Unct. : ' + str(round(unct, 2)) +'%'
    plt.hist(dispersion_all,color='blue', bins=80, label=text1)
    #plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))


    plt.xlim(0.197, 0.200)
    plt.grid(True, linestyle=':')
    plt.title('Histogram of TEMPO Spectral Sampling Estimates', fontweight='bold', **csfont)
    plt.ylabel('Frequency', fontsize=12, **hsfont)
    plt.xlabel('Spectral Sampling Estimate (nm/pixel)', fontsize=12, **hsfont)
    L = plt.legend()
    plt.setp(L.texts, family='Times New Roman')
    plt.legend(loc='best', fontsize=12)
    plt.show()
    cc

    print(np.mean(dispersion_all))
    print(np.std(dispersion_all))

def calculate_dispeersion_nearest(save_dir, data_path_all, wavelen):
    dispersion_all = []
    wavelen_all = []
    Visible = data_path_all[0:10]
    UV = data_path_all[10:]
    channels = [Visible, UV]
    #print(Visible)
    for bands in channels:


        for i in range(0, len(bands)-1):

            dframe_ref = pd.read_csv(save_dir+'/Iteration1/'+bands[i], delimiter=',')
            wavelen_ref = float(bands[i].split('_')[-1].split('.csv')[0])
            dframe_comp = pd.read_csv(save_dir+'/Iteration1/'+bands[i+1], delimiter=',')
            wavelen_comp =  float(bands[i+1].split('_')[-1].split('.csv')[0])
            print(wavelen_ref, wavelen_comp)


            delta_wavelen = wavelen_ref - wavelen_comp

            delta_pixel = (dframe_ref['Peak Loc.'] -  dframe_comp['Peak Loc.']).values
            #print(delta_pixel)
            #print(delta_wavelen)

            dispersion = np.true_divide(delta_wavelen, delta_pixel)
            #print(dispersion)
            Mean_val = 'Mean = ' +str(np.round(np.mean(dispersion), 4))
            stdev = 100*np.std(dispersion)/np.mean(dispersion)
            unct = '\n Unct = ' + str(np.round(stdev,2)) +'%'
            label = Mean_val + unct


            plt.plot(dispersion,'b.', label=label)
            dispersion_all.append(dispersion)
            wavelen_all.append(wavelen_comp)
            plt.ylim(0.196, 0.2)
            plt.title('Spectral Sampling Estimates for '
                      + str(wavelen_comp) + 'nm\n'
                       + '(Reference WL:' +str(wavelen_ref)+ 'nm)', fontweight='bold')
            legend = plt.legend(loc='best', shadow=True)
            legend.get_frame().set_facecolor('orange')
            #plt.legend(loc='best')
            plt.xlabel('Spatial Pixel Indices', fontsize=12)
            plt.ylabel('Spectral Sampling (nm/pixel)', fontsize=12)
            plt.grid(True, linestyle=':')
            #plt.show()
            plt.close('all')

    mean_dispersion = np.mean(np.array(dispersion_all), axis=1)
    Mean_val = np.mean(mean_dispersion[:-1])
    Unct = 100*np.std(mean_dispersion[:-1])/Mean_val
    text = 'Mean (Spectral Dir.) = ' + str(round(Mean_val, 3)) + '\nUnct (Spectral Dir.) = ' + str(round(Unct,3)) +'%'

    std_dispersion = 100*np.std(np.array(dispersion_all), axis=1)/mean_dispersion
    plt.errorbar(wavelen_all[:-1], mean_dispersion[:-1], yerr = std_dispersion[:-1],
                 mfc='red',mec='red', fmt='s', ms=7, elinewidth=2, ecolor='b', 
                 label=text)
#    plt.errorbar(wavelen_all[:-1], mean_dispersion[:-1], yerr = 0,
#                 mfc='red',mec='red', fmt='s--', ms=7, elinewidth=2, ecolor='b', 
#                 linecolor ='red', label=text)
    plt.annotate(text, xy=(178, 250), xycoords='axes points',
                 size=10, ha='right', va='top',
                 bbox=dict(boxstyle='round', fc='w'))
    plt.grid(True, linestyle=':')
    plt.title('Spectral Sampling Estimates Vs. Wavelength', fontweight='bold')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Spectral Sampling (nm/pixel)', fontsize=12)
    plt.show()






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

    data_path_all = sorted([each for each in os.listdir(os.path.join(save_dir, 'Iteration1'))
                            if each.startswith('Gaussian_parameters')], reverse=True)
    wavelen = ([data_path.split('_')[-1].split('.csv')[0] for data_path in data_path_all])
    #calculate_dispersion_all (save_dir, data_path_all, wavelen)
    calculate_dispeersion_nearest(save_dir, data_path_all, wavelen)














if __name__ == "__main__":
    main()
