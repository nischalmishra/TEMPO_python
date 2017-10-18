# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""
# This script uses the PRNU map derived from integration time sweep data to some
# of the images acquired during intesnity sweep.

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from scipy.io.idl import readsav
from create_PRNU_map import filter_outlier_median,\
                             perform_bias_subtraction_ave, create_image,\
                             create_striping_metric, create_striping_metric_plot

def create_final_PRNU_mask():
    """ Function to create a final PRNU mask. This bascially fills the pixels
        not used during PRNU mask to one and creates a full frame image"""
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Intensity_Sweep\Script_Data\saved_quads\PRNU\Final_PRNU'
    map_file = 'Final_PRNU.csv'
    PRNU_file = os.path.join(file_path, map_file)
    PRNU = genfromtxt(PRNU_file, delimiter=',')

    active_quad_A = np.ones((1024, 1024))
    active_quad_B = np.ones((1024, 1024))
    active_quad_C = np.ones((1024, 1024))
    active_quad_D = np.ones((1024, 1024))
    # for the pixels not use in the generation of mask use fill value of 1
     #remember, this is what I did
     #lower_quads = np.concatenate((active_quad_A[0:900, :], np.fliplr(active_quad_B[0:900, :])), axis=1)
     #upper_quads = np.concatenate((np.flipud(active_quad_D[0:900, :]), np.rot90(active_quad_C[0:900, :], 2)), axis=1)

    active_quad_A[0:1000, 0:1024] = PRNU[0:1000, 0:1024] # quad A
    active_quad_B[0:1000, 0:1024] = PRNU[0:1000, 1024:2048] # quad B
    active_quad_C[24:, 0:1024] = PRNU[1000:, 1024:2048] # quad C
    active_quad_D[24:, 0:1024] = PRNU[1000:, 0:1024] # quad D
    

    lower_quads = np.concatenate((active_quad_A, active_quad_B), axis=1)
    upper_quads = np.concatenate((active_quad_D, active_quad_C), axis=1)
    final_PRNU_mask = np.concatenate((lower_quads, upper_quads), axis=0)
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(final_PRNU_mask, cmap='nipy_spectral', origin='lower')
    plt.title('TEMPO PRNU Mask')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.savefig(file_path+'/'+'Final_prnu_mask_filled.png', dpi=100, bbox_inches="tight")
    plt.close('all')
    csv_file_name = file_path+'/'+ 'Full_Frame_PRNU_Intensity_sweep.csv'
    np.savetxt(csv_file_name, np.array(final_PRNU_mask), delimiter=",")
    print(np.array(final_PRNU_mask).shape)
    return csv_file_name

def main():
	# This is the main function
    #mask = create_final_PRNU_mask ()
    #cc
    # read PRNU mask file
    PRNU_path = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU\Final_PRNU'
    PRNU_file_name = 'Full_Frame_PRNU_Integ_sweep.csv'
    PRNU_file = os.path.join(PRNU_path, PRNU_file_name)
    PRNU_mask = genfromtxt(PRNU_file, delimiter=',')

    # ok, now let us apply PRNU Mask to the intensity sweep data
    file_path = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    if 'Integration_Sweep' in file_path:
        saturated_collects = ['FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']

    elif 'Intensity_Sweep' in file_path:
        saturated_collects = ['162_OP_INT_118000.dat.sav', '164_OP_INT_118000.dat.sav',
                              '166_OP_INT_118000.dat.sav', '168_OP_INT_118000.dat.sav',
                              '170_OP_INT_118000.dat.sav', '172_OP_INT_118000.dat.sav',
                              '174_OP_INT_118000.dat.sav', '176_OP_INT_118000.dat.sav',
                              '178_OP_INT_118000.dat.sav', '180_OP_INT_118000.dat.sav',
                              '182_OP_INT_118000.dat.sav', '184_OP_INT_118000.dat.sav',
                              '186_OP_INT_118000.dat.sav', '188_OP_INT_118000.dat.sav',
                              '190_OP_INT_118000.dat.sav', '192_OP_INT_118000.dat.sav',
                              '194_OP_INT_118000.dat.sav', '196_OP_INT_118000.dat.sav',
                              '198_OP_INT_118000.dat.sav', '200_OP_INT_118000.dat.sav',
                              '202_OP_INT_118000.dat.sav', '254_OP_INT_118000.dat.sav',
                              '252_OP_INT_118000.dat.sav', '250_OP_INT_118000.dat.sav',
                              '248_OP_INT_118000.dat.sav', '246_OP_INT_118000.dat.sav',
                              '244_OP_INT_118000.dat.sav', '242_OP_INT_118000.dat.sav',]


    all_int_files = [each for each in os.listdir(file_path) \
                     if each.endswith('.dat.sav')]

    nominal_int_files = [items for items in all_int_files if items not in saturated_collects]

    save_dir = 'PRNU/PRNU_Validation_Intensity_Sweep'
    orig_image_dir = 'saved_quads'
    plot_dir = os.path.join(file_path, save_dir, orig_image_dir)
    if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    for data_files in nominal_int_files:

            data_path_name_split = data_files.split('_')
            print(data_files)

            if 'Intensity_Sweep' in file_path:
                 int_time = data_path_name_split[0]
                 string1 = 'VA_'
                 string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '

            data_file = os.path.join(file_path, data_files)
            IDL_variable = readsav(data_file)
            all_full_frame = IDL_variable.q
            # identify the regions to create PRNU map and create a sample image
            quad_A = np.mean(all_full_frame[:, 0, 4:1028, 10:1034], axis=0)
            tsoc_A = np.mean(all_full_frame[:, 0, 4:1028, 1034:1056], axis=0)
            active_quad_A = perform_bias_subtraction_ave(quad_A, tsoc_A)

            quad_B = np.mean(all_full_frame[:, 1, 4:1028, 10:1034], axis=0)
            tsoc_B = np.mean(all_full_frame[:, 1, 4:1028, 1034:1056], axis=0)
            active_quad_B = perform_bias_subtraction_ave(quad_B, tsoc_B)

            quad_C = np.mean(all_full_frame[:, 2, 4:1028, 10:1034], axis=0)
            tsoc_C = np.mean(all_full_frame[:, 2, 4:1028, 1034:1056], axis=0)
            active_quad_C = perform_bias_subtraction_ave(quad_C, tsoc_C)

            quad_D = np.mean(all_full_frame[:, 3, 4:1028, 10:1034], axis=0)
            tsoc_D = np.mean(all_full_frame[:, 3, 4:1028, 1034:1056], axis=0)
            active_quad_D = perform_bias_subtraction_ave(quad_D, tsoc_D)

            lower_quads = np.concatenate((active_quad_A, active_quad_B), axis=1)
            upper_quads = np.concatenate((active_quad_D, active_quad_C), axis=1)
            active_quads = np.concatenate((lower_quads, upper_quads), axis=0)

            # The following line applied PRNU correction to the active quads
            active_quads = np.true_divide(active_quads, PRNU_mask)

            # Ok, now let's generate a PRNU corrected image
            figure_name = plot_dir+'/'+ string1 + str(int_time)+'.png'
            title = 'Active Quad Image (PRNU Applied), ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
            create_image(active_quads, title, figure_name)

            # Ok, now lets us generate a striping metric of the PRNU corrected image
            active_quads_filtered = filter_outlier_median(active_quads)
            striping_directory = 'striping_metric_hist'
            striping_met_directory = os.path.join(file_path, save_dir, striping_directory)
            if not os.path.exists(striping_met_directory):
                os.makedirs(striping_met_directory)
            figure_name = striping_met_directory +'/'+ string1+ str(int_time)+'.png'
            title = 'Striping metric (PRNU Applied), ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
            striping_metric_all_rows, striping_metric_all_cols = create_striping_metric(active_quads_filtered, title, figure_name)

                # now the actual plot
            striping_plot_directory = 'striping_metric_plots'
            striping_met_plot_directory = os.path.join(file_path, save_dir, striping_plot_directory)
            if not os.path.exists(striping_met_plot_directory):
                os.makedirs(striping_met_plot_directory)
            figure_name = striping_met_plot_directory +'/'+ string1+ str(int_time)+'.png'
            title = 'Striping metric plot, ' + string2+ str(int_time) + ', image size = '+ str(active_quads.shape)
            create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name)
            #cc
if __name__ == "__main__":
    main()
    