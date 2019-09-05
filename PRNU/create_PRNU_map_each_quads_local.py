# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

The purpose of this code is to generate PRNU map for each TEMPO quads.
Over the course of time this script will be modified as we learn more about the
TEMP CCD behavior. This script takes a 3 by 3 area and slides it along the image
to calculate PRNU as Det.DN/Local Mean

"""
import os
import numpy as np
from scipy.ndimage.filters import uniform_filter, median_filter
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import pandas as pd


#*****************************************************************************

def perform_bias_subtraction_ave(active_quad, trailing_overclocks):

    """This function calculates the offset with trailing overclocks and subtracts
    the offset from the active area
    """
    spec_pix, spat_pix = active_quad.shape
    bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
    even_detector_bias = trailing_overclocks[ :, ::2]
    # Toss out the outliers
    even_detector_bias = even_detector_bias[:, 4:]
    avg_bias_even = np.mean(even_detector_bias, axis=1)
    odd_detector_bias = trailing_overclocks[:, 1::2]
    # Toss out the outliers
    odd_samples = odd_detector_bias[:, 4:]
    #rows, cols = odd_samples.shape
    #print(rows, cols)
    cols = odd_samples.shape[1]
    #print (cols)
    odd_detector_bias = odd_samples[:, 0:cols-1]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]
    odd_detector_active_quad = active_quad[:, 1::2]
    #Subtract off the offset from active region
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (spec_pix, spat_pix))
    # Reform the quad
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad

def perform_linearity(active_quad, quad_name):
    """This function does the non-linearity correction basedupon the look up
    table. The active quads needs to be broken up in even and odd as linearity table is
    created per octants
    """
    quad = quad_name.split(' ')[1]
    active_quad[active_quad < 0] = 0
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file = 'Linearity_look_up_table_final_DN.csv'
    linearity_file = pd.read_csv(path+'/'+ linearity_file)
    # for even quad
    lut_even = linearity_file[['Input Signal (DN)', quad+str(1)]].dropna()
    lut_even['DN'] = lut_even.pop('Input Signal (DN)').astype(int)
    lut_even['Value'] = lut_even.pop(quad+str(1)).astype(float)
    lut_even = lut_even.set_index('DN')
    even_detector_active_quad = active_quad[:, ::2]
    nx_half_quad, ny_half_quad = even_detector_active_quad.shape
    even_detector_active_quad = np.reshape(np.array(even_detector_active_quad),
                                           (nx_half_quad*ny_half_quad, 1))
    dframe = pd.DataFrame(data=even_detector_active_quad)
    dframe.columns = ['DN']
    datin_even = dframe['DN'].astype(int)
    linearized_quad_even = datin_even.apply(lambda x: lut_even.ix[x])
    even_pixels_quad = np.reshape(linearized_quad_even.values, (nx_half_quad,
                                                                ny_half_quad))
    dframe = None
    #for odd
    lut_odd = linearity_file[['Input Signal (DN)', quad+str(2)]].dropna()
    lut_odd['DN'] = lut_odd.pop('Input Signal (DN)').astype(int)
    lut_odd['Value'] = lut_odd.pop(quad+str(2)).astype(float)
    lut_odd = lut_odd.set_index('DN')
    odd_detector_active_quad = active_quad[:, 1::2]

    nx_half_quad, ny_half_quad = odd_detector_active_quad.shape
    odd_detector_active_quad = np.reshape(np.array(odd_detector_active_quad),
                                          (nx_half_quad*ny_half_quad, 1))
    dframe = pd.DataFrame(data=odd_detector_active_quad)
    dframe.columns = ['DN']
    datin_odd = dframe['DN'].astype(int)
    linearized_quad_odd = datin_odd.apply(lambda x: lut_odd.ix[x])
    odd_pixels_quad = np.reshape(linearized_quad_odd.values,
                                 (nx_half_quad, ny_half_quad))
    # now let's put quad in place
    nx_quad, ny_quad = active_quad.shape
    linearized_quad = np.array([[0]*ny_quad]*nx_quad)
    linearized_quad = np.reshape(linearized_quad, (nx_quad, ny_quad))
    linearized_quad[:, ::2] = even_pixels_quad
    linearized_quad[:, 1::2] = odd_pixels_quad
    return linearized_quad

def perform_smear_subtraction(active_quad, int_time):
    """the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    """
    frame_transfer = 8*10**(3)
    smear_factor = (frame_transfer / (int_time+ frame_transfer))* np.mean(active_quad, axis=0)
    #print(active_quad.shape)
    #print(smear_factor.shape)
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad

def calculate_dark_current(image, i, int_time):
    """ Calculate the dark current based off the dark data
    It takes the filename of Light data and searches for the mathing integration
    time   in the dark data directory, find the required dark data, subtracts off
    the offset and smear and computes the dark current
    """
    dark_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark'
    data_path_name_split = image.split('_')
    #print(data_path_name_split)
    all_int_files = [each for each in os.listdir(dark_data_dir) \
                         if each.endswith('_'+data_path_name_split[-1])]
    print(all_int_files)

    dark_data_file = os.path.join(dark_data_dir, all_int_files[0])
    IDL_variable = readsav(dark_data_file)
    all_full_frame = IDL_variable.q
    quad = all_full_frame[:, i, :, :]
    active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)
    tsoc = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
    bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc)
    smear_subtracted_quad = perform_smear_subtraction(bias_subtracted_quad[10:1000, :],
                                                      int_time)
    return smear_subtracted_quad

def calculate_window_mean(image):
    """
    This function creates a mean mask. Mean is a 3*3 moving average mean
    """
    #rows, cols = image.shape
    mean_mask = []
    image = image
    mean_mask = median_filter(image, size=(5, 5), mode='reflect')     
    return mean_mask


def fill_vignetted_pixels (subset_quad):
    print(subset_quad.shape)
    rows, cols = subset_quad.shape
    final_image = np.zeros((1028, 1024))
    final_image = np.array(final_image)
    final_image[0:10, :] = np.flipud(subset_quad[0:10, :])
    final_image[10:1000, :] = subset_quad
    final_image[1000:, :] = np.flipud(subset_quad[962:, :])
    return final_image
    

def interpolate_vignetted_pixels(subset_quad):
    # Let's interpolate forward and backward to fill the vginetted rows
    
    #active_quad_A[15:990, :] slice taken for PRNU
    #active_quad_A[10:1000, :]
    rows, cols = subset_quad.shape
    print(rows, cols)
    x1_rows = list(range(10, rows+10))    
    new_points_end = list(range(rows+10-1, 1028))   
    new_points_begin = list(range(0, 11))    
    y2_all_end = []
    y2_all_begin = []
    for i in range(150, cols):
        fit = np.polyfit(x1_rows[900:], subset_quad[900:, i], 5)   
        line = np.poly1d(fit)
        val = np.polyval(fit, x1_rows[900:])     
        desired_val_end = line(new_points_end)
        desired_val_begin = line(new_points_begin)
        plt.figure()
        plt.plot(x1_rows[900:], val,'r.', label='5th order polyfit')
        plt.plot(x1_rows[900:], subset_quad[900:, i],'go--', label='actual data')
        #plt.plot(new_points_begin, desired_val_begin,'bo--', label='extrapolated (beginning)')
        plt.plot(new_points_end, desired_val_end,'mo--', label='extrapolated (end)')
        plt.grid(True, linestyle=':')
        plt.legend(loc='best', ncol=1)
        plt.title('Example of extrapolation to fill vignetted pixels')
        plt.ylabel('Raw Image - Offset - SMEAR (DN)')
        plt.xlabel('Spectral Indices (#)')
        plt.show()
        
        plt.figure()
        plt.plot(100*(val-subset_quad[:, i])/subset_quad[:, i],'k.')
        plt.ylabel('%Diff')
        plt.xlabel('Spectral Indices (#)')
        plt.title('Residue of model and the data')
        plt.ylim([-1, 1])
        plt.grid(True, linestyle=':')
        plt.show()
#        cc
        y2_all_end.append(desired_val_end)
        y2_all_begin.append(desired_val_begin)
    y2_all_end = np.array(y2_all_end).T
    y2_all_begin = np.array(y2_all_begin).T  
    #print(y2_all_begin.shape)
    final_mask = np.zeros((1028, 1024))
    final_mask[0:16, :] = y2_all_begin
    final_mask[16:rows+16, :] = subset_quad
    final_mask[rows+14:, :] = y2_all_end  
    #print(final_mask.shape)
#    plt.plot(final_mask[:, 100],'r')
#    plt.plot(final_mask[:, 500], 'g')
#    plt.show()
    #cc
    return final_mask


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
    """ Create image as a sanity check step
    """
    plt.figure()
    axes_ = plt.gca()
    image = axes_.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(axes_)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    #plt.show()
    #cc
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')

def create_hist(prnu_map, title, figure_name):
    """ Create histogram of PRNU to understand the distribution
    """
    nx_quad, ny_quad = prnu_map.shape
    label = 'Mean = '+ str(round(np.mean(prnu_map), 3)) + \
            '\n Median = '+ str(round(np.median(prnu_map), 3))
            #'\n Std. = '+ str(round(np.std(PRNU_map), 3))+ \
            #'\n Max = '+ str(round(np.max(PRNU_map), 3)) + \
            #'\n Min = '+ str(round(np.min(PRNU_map), 3))
    plt.hist(np.reshape(prnu_map, (nx_quad* ny_quad, 1)), 400, facecolor='red',
             label=label)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                        prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.xlim(0.95, 1.05)
    plt.ylabel('Frequency (# of pixels)', fontsize=12, fontweight="bold")
    plt.xlabel('PRNU', fontsize=12, fontweight="bold")
    plt.title(title)
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')


def plot_few_PRNU_vals_spectral(PRNU_map, figure_name, title):
    """ Plot PRNU in spectral direction for few samples
    """
    plt.figure()
    plt.plot(PRNU_map[:, 10], 'r', label='Spatial Index 1')
    plt.plot(PRNU_map[:, 55], 'k', label='Spatial Index 55')
    plt.plot(PRNU_map[:, 100], 'g', label='Spatial Index 100')
    plt.plot(PRNU_map[:, 381], color='orange', label='Spatial Index 381')
    plt.plot(PRNU_map[:, 500], 'b', label='Spatial Index 500')
    plt.plot(PRNU_map[:, 1000], 'm', label='Spatial Index 1000')
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.title(title)
    plt.ylabel('PRNU Value (Ratio)')
    #plt.ylim([0.90, 1.05])
    plt.xlabel('Pixel Indices (#)')
    #plt.show()
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')
    #plt.show()

def plot_few_PRNU_vals_spatial(PRNU_map, figure_name, title):
    """ Plot PRNU in spatiall direction for few samples
    """
    plt.figure()
    plt.plot(PRNU_map[:, 10], 'r', label='Spectral Index 1')
    plt.plot(PRNU_map[:, 55], 'k', label='Spectral Index 55')
    plt.plot(PRNU_map[:, 100], 'g', label='Spectral Index 100')
    plt.plot(PRNU_map[:, 381], color='orange', label='Spectral Index 381')
    plt.plot(PRNU_map[:, 500], 'b', label='Spectral Index 500')
    plt.plot(PRNU_map[:, 900], 'm', label='Spectral Index 1000')
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.title(title)
    plt.ylabel('PRNU Value (Ratio)')
    #plt.ylim([0.90, 1.05])
    plt.xlabel('Pixel Indices (#)')
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')
    #plt.show()

def main():
    """
    Tte main function
    """
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_map_Median\5_by_5'
    if 'Integration_Sweep' in file_path1:
        not_used_collects = ['FT6_SHORT_INT_25048.dat.sav', 'FT6_SHORT_INT_0.dat.sav',
                             'FT6_SHORT_INT_10038.dat.sav', 'FT6_SHORT_INT_19981.dat.sav',
                             'FT6_SHORT_INT_50000.dat.sav', 'FT6_SHORT_INT_54971.dat.sav',
                             'FT6_SHORT_INT_34990.dat.sav', 'FT6_SHORT_INT_30019.dat.sav',
                             'FT6_SHORT_INT_45028.dat.sav', 'FT6_SHORT_INT_4971.dat.sav',
                             'FT6_SHORT_INT_39961.dat.sav', 'FT6_SHORT_INT_54971.dat.sav',
                             'FT6_SHORT_INT_15009.dat.sav', 'FT6_SHORT_INT_10038.dat.sav'
                             'FT6_SHORT_INT_19981.dat.sav', 'FT6_SHORT_INT_15009dat.sav',
                             'FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                             'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                             'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                             'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                             'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                             'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                             'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                             'FT6_LONG_INT_199999.dat.sav', 'FT6_LONG_INT_125047.dat.sav']

    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
    nominal_int_files = [items for items in all_int_files if items not in not_used_collects]
    # let's work with each quad. PRNU map will be generated for each quad
    for i in range(0, 4):
        all_prnu_map = []
        linearity=0
        for data_files in nominal_int_files:
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            print(data_file)
            IDL_variable = readsav(data_file)
            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                string1 = 'VA_'
                string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '

            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q
            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)
            tsoc_A = np.mean(quad[:, 2:1028, 1034:1056], axis=0)
            # Let's remove the offset
            print('Ok, Now remove Offset using trailing overclocks')
            active_quad_A = perform_bias_subtraction_ave(quad_A[2:1028, 10:1034],
                                                         tsoc_A)
#            plt.imshow(active_quad_A, cmap='bwr', origin='lower')
#            plt.show()
            # removed the vignetted zones and perform linearity
            print('Ok, Now apply non-linearity correction')
            if linearity:
                linearized_quad_A = perform_linearity(active_quad_A[10:1000, :], quads[i])
            else:
                linearized_quad_A = active_quad_A[10:1000, :]
            print('Ok, Now apply SMEAR Removal')
            active_quad_A = perform_smear_subtraction(linearized_quad_A, int_time)
            print('OK, Now remove the dark current')
            dark_current = calculate_dark_current(data_files, i, int_time)
            active_quad_A = active_quad_A - dark_current
            #active_quad_A[active_quad_A<0] = 0
            active_quad_A = fill_vignetted_pixels(active_quad_A)
           
            # PRNU map generation begins from here
            # Window size is input. If it's 3. it will create a 3*3 moving average filter
            mean_mask = calculate_window_mean(active_quad_A)
            prnu_single_image = active_quad_A/mean_mask   
           # prnu_single_image[prnu_single_image>1.2] = 1.2
            #prnu_single_image[prnu_single_image<0.8] = 0.8 
            #plt.plot(prnu_single_image[:, 10],'k--')
            #prnu_single_image =  fill_vignetted_pixels(prnu_single_image)
            #create_image(prnu_single_image,'a','b')
               
            #******************************************************************
            quad_dir = quads[i]
            orig_image_dir = 'saved_quads'
            plot_dir = os.path.join(save_dir, quad_dir, orig_image_dir)
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)

           #Let's plot the image of the active quad
            figure_name = plot_dir + '/' + string1 + str(int_time) +'_new_small.png'
            title = 'Active Pixels in ' + quads[i]+', ' + string2 + str(int_time)
            create_image(active_quad_A, title, figure_name)

            orig_image_dir = 'saved_mean_mask'
            plot_dir = os.path.join(save_dir, quad_dir, orig_image_dir)
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)

           #Let's plot the image of the mean mask
            figure_name = plot_dir + '/' + string1 + str(int_time) +'_new_small.png'
            title = 'Example of Median Filter Mask\n ' + quads[i]+', ' + string2 + str(int_time)
            create_image(np.array(mean_mask), title, figure_name)

            PRNU_plot_name = 'PRNU_plots_spectral_dir'
            PRNU_dir = os.path.join(save_dir, quad_dir, PRNU_plot_name)
            if not os.path.exists(PRNU_dir):
                os.makedirs(PRNU_dir)
            figure_name = PRNU_dir + '/' + string1 + str(int_time) +'_PRNU_plot'
            title = 'PRNU values along spectral direction\n '+ quads[i]+', ' + \
                    string2 + str(int_time)+ r" $\mu$" +'secs'
            plot_few_PRNU_vals_spectral(prnu_single_image, figure_name, title)

            PRNU_plot_name = 'PRNU_plots_spatial_dir'
            PRNU_dir = os.path.join(save_dir, quad_dir, PRNU_plot_name)
            if not os.path.exists(PRNU_dir):
                os.makedirs(PRNU_dir)
            figure_name = PRNU_dir + '/' + string1 + str(int_time) +'_PRNU_plot'
            title = 'PRNU values along satial direction\n '+ quads[i] +', ' + \
                     string2 + str(int_time)+ r" $\mu$" +'secs'
            plot_few_PRNU_vals_spatial(prnu_single_image.T, figure_name, title)

            # lets' plot the PRNU map
            # Let's create a separate directory
            PRNU_image_dir = 'PRNU_map'
            PRNU_hist_dir = 'PRNU_hist'
            image_dir = os.path.join(save_dir, quad_dir, PRNU_image_dir)
            if not os.path.exists(image_dir):
                os.makedirs(image_dir)
            hist_dir = os.path.join(save_dir, quad_dir, PRNU_hist_dir)
            if not os.path.exists(hist_dir):
                os.makedirs(hist_dir)

            #plot the PRNU map
            figure_name = image_dir + '/' + string1 + str(int_time) + '_new_small.png'
            title = 'PRNU Map, '+ quads[i]+', ' + string2 + str(int_time) + \
                     r" $\mu$" +'secs'
            create_image(prnu_single_image, title, figure_name)
            csv_file_name = image_dir+ '/' + quads[i] + string1 + \
                            str(int_time) + '_new_small.csv'
            np.savetxt(csv_file_name, np.array(prnu_single_image), delimiter=",")

            # plot the histograms of PRNU map
            figure_name = hist_dir + '/' + string1 + str(int_time) + '._new_small.png'
            title = 'Histogram of PRNU Map, ' + quads[i] + ', ' + string2 + \
                    str(int_time) + r" $\mu$" +'secs'
            create_hist(prnu_single_image, title, figure_name)
            all_prnu_map.append(prnu_single_image)
            

        # Ok, take average of all the PRNU maps and compute a final mask
        
        final_PRNU = np.mean(np.array(all_prnu_map), axis=0)        
        
        final_PRNU_std = np.std(np.array(all_prnu_map), axis=0)

        PRNU_ave_dir = 'Final_PRNU/Full_quad_PRNU'
        quad_dir = quads[i]
        final_image_dir = os.path.join(save_dir, quad_dir, PRNU_ave_dir)
        if not os.path.exists(final_image_dir):
            os.makedirs(final_image_dir)

        title = 'Average PRNU Map, ' + quads[i]
        figure_name = final_image_dir+'/'+ 'final_mask_image'+'.png'
        create_image(final_PRNU, title, figure_name)

        title = ' Uncertainty associated with average PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'unct_final_mask_image'+'.png'
        create_image(100* final_PRNU_std/final_PRNU, title, figure_name)

        title = ' Histogram of Average PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'final_mask_hist'+'.png'
        create_hist(final_PRNU, title, figure_name)

        title = ' Histogram of Uncertainty associated with final PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'unct_final_mask_hist'+'.png'
        #create_hist(100* final_PRNU_std/final_PRNU, title, figure_name)

        csv_file_name = final_image_dir+'/'+ quads[i]+'_Final_PRNU.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU), delimiter=",")
        csv_file_name = final_image_dir+'/'+ quads[i]+'_Final_PRNU_Std.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU_std), delimiter=",")


if __name__ == "__main__":
    main()
