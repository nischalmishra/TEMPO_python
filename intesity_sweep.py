# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 09:56:05 2019

@author: nmishra
"""
import os
from scipy.io.idl import readsav
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def perform_bias_subtraction(active_quad, trailing_overclocks):
    """
    Subtract out the bias
    """
    # sepearate out even and odd detectors
    spec_pix, spat_pix = active_quad.shape
    bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
    odd_detector_bias = trailing_overclocks[:, ::2]
    odd_detector_bias = odd_detector_bias[:, 4:]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    # Now repeat the process for even lines
    even_detector_bias = trailing_overclocks[:, 1::2]
    even_detector_bias = even_detector_bias[:, 4:]
    cols = even_detector_bias.shape[1]
    even_detector_bias = even_detector_bias [:, 0:cols-1]
    avg_bias_even = np.mean(even_detector_bias, axis=1)
    # Now subtract the bias

    odd_detector_active_quad = active_quad[:, ::2]
    even_detector_active_quad = active_quad[:, 1::2]
    # None keyword helps to subtract 1-D array from 2D array row wise or columnwise
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (spec_pix, spat_pix))
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_odd
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_even
    return bias_subtracted_quad



def create_image(image_data):
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

    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
    #plt.show()
    #plt.pause(0.1)
    #plt.show()
#    plt.draw()
#    plt.pause(0.001)
#    print('OK, Move forward')
    #plt.show(block=False)
    #plt.close('all')



def main():
    """"
    The main function
    """
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Intensity_Sweep\Light\saved_quads'
    all_int_files = [each for each in os.listdir(file_path1) \
                         if each.endswith('248_OP_INT_118000.dat.sav')]
    #print(all_int_files)
    if 'Integration_Sweep' in file_path1:
        saturated_collects = ['FT6_LONG_INT_130018.dat.sav',#'FT6_SHORT_INT_0.dat.sav',
                              'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']
    elif 'Intensity_Sweep' in file_path1:
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
                              '202_OP_INT_118000.dat.sav']
    nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files]
    #print(nominal_int_files)

    for data_files in nominal_int_files:
        print(data_files)
        data_file = os.path.join(file_path1, data_files)
        idl_variable = readsav(data_file)
        all_full_frame = idl_variable.q
        all_quads = []
        for i in range(0, 4):
            quad = np.mean(all_full_frame[:, i, :, :], axis=0)
            tsoc_all = quad[2:1030, 1034:1056]
            active_quad = quad[2:1030, 10:1034]
            bias_subtracted = perform_bias_subtraction(active_quad, tsoc_all)
            all_quads.append(bias_subtracted)
        lower_quads = np.concatenate((all_quads[0], np.fliplr(all_quads[1])),
                                     axis=1)
        upper_quads = np.concatenate((np.flipud(all_quads[3]),
                                      np.rot90(all_quads[2], 2)), axis=1)
        active_quads = np.concatenate((lower_quads, upper_quads), axis=0)
        print(active_quads.shape)
        spatial_dir = np.arange(0, active_quads.shape[1])
        plt.plot(spatial_dir, active_quads[200, :], 'b.', linewidth=0.3)
        plt.grid(True, linestyle=':')
        plt.title('Signal Counts Along Spatial Direction')
        plt.ylabel('Bias Subtracted Counts (DN)')
        plt.xlabel('Spatial Pixel Locations')

        plt.show()
        #plt.close('all')
        #create_image(active_quads)

if __name__ == "__main__":
    main()
