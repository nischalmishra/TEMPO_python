# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 10:34:21 2017

@author: nmishra

This routine is used to process TEMPO Spectrometer Data. While the function
are not in order, follow following steps to process the image

1. Mask the outliers in the active region
2. Perform Offset subtraction
3. Perform Non-Linearity_correction
4. Subtract out the SMEAR
5. Apply PRNU Map
6. Apply the Cross Talk Correction


"""

import os
import numpy as np
#import numpy.ma as ma
from scipy.io.idl import readsav
import pandas as pd
import h5py
#import matplotlib.pyplot as plt
#pylint: disable= E1101

# To Do: Add Dark current removal correction. Waiting for file structures



def read_idl_file(data_file):
    """Read IDL variables. Dave's IDL code to read spectrometer data runs
     much faster than my Python script. Also if there are .sav files,
     I can use them in  future rather than reading the whole file again
     """
    read_variable = readsav(data_file)
    try:
        full_frame_image = read_variable.quads
    except KeyError:
        full_frame_image = read_variable.q
    return full_frame_image


def parse_telemetry_file(file_path, telemetry_file_name):
    """ Read the variables required from telemetry files to process the spectro
    meter image. Values required are integration time, coadds and temperatures
    """
    h5file_name = telemetry_file_name + '.h5'
    hdf_name = os.path.join(file_path, h5file_name)
    with h5py.File(hdf_name, 'r') as hdf:
        group1 = hdf.get('telemetry')
        group2 = group1.get('calculated_values')
        coadds = group2.attrs['ncoadds']
        int_time = group2.attrs['tint']
        ccd_temp = group2.attrs['ccd_temp']
        fpe_temp = group2.attrs['fpe_temp']
        # sanity check
        #print(coadds, int_time, ccd_temp, fpe_temp)
        #cc

    return coadds, int_time, fpe_temp, ccd_temp

def perform_coaddition_correction(file_name, num_coadds):
    """ This module dooes the coaddition correction. The coaddtion info is read
    from the telemetry files"""
    corrected_data = file_name/num_coadds
    return corrected_data


def ping_pong_identification_fps():
    """ During FPS it was observed that even lines of the trailing overclocks
    where greated than odd lines.The same check needs to be done for all the
    test images to ensure correct LUTs like linearity and electronics gains
    are applied to the image
    """
    # The first line is going to be odd line. So for all FPS images, odd-even
    #is less than one. Hence return all negatives for each quads
    sign_fps = [-1, -1, -1, -1]
    return sign_fps


def calculate_mean_offset_diff(raw_quads):
    """ This function calculates the median offset of each quads based on trailing
    overclock pixels. This is to check if the ping/pong changed from FPS
    """
    num_quads = raw_quads.shape[0]
    sign_spectro = []
    for quads in range(0, num_quads):
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]
        odd_detector_bias = trailing_overclocks[ :, ::2]
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        odd_detector_bias = odd_detector_bias[:, 4:]
        avg_bias_odd = np.mean(odd_detector_bias, axis=1)

        # Now repeat the process for even lines
        even_detector_bias = trailing_overclocks[:, 1::2]
        even_detector_bias = even_detector_bias[:, 4:]
        cols = even_detector_bias.shape[1]
        even_detector_bias = even_detector_bias [:, 0:cols-1]
        avg_bias_even = np.mean(even_detector_bias, axis=1)
        if np.mean(avg_bias_odd) - np.mean(avg_bias_even) < 0:
            sign = -1
        else:
            sign = 1
        sign_spectro.append(sign)

    return sign_spectro

def read_outlier_mask():
    """ The outlier mask is created based on dark data. The idea is to create
    masked array of active quads and not use these in further calculations
    """

    file_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask'+\
                r'Outlier_median\save_masks_logical_OR\to_SAO'
    mask_a = np.genfromtxt(file_path + '/' + 'quad_A_outlier_mask.csv',
                           delimiter=',')
    mask_b = np.genfromtxt(file_path + '/' + 'quad_B_outlier_mask.csv',
                           delimiter=',')
    mask_c = np.genfromtxt(file_path + '/' + 'quad_C_outlier_mask.csv',
                           delimiter=',')
    mask_d = np.genfromtxt(file_path + '/' + 'quad_D_outlier_mask.csv',
                           delimiter=',')
    outlier_mask = [mask_a, mask_b, mask_c, mask_d]
    return np.array(outlier_mask)


def perform_bias_subtraction(raw_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads = []
    num_quads = raw_quads.shape[0]

    for quads in range(0, num_quads):
        active_quad = raw_quads[quads, 2:1030, 10:1034]
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]
        odd_detector_bias = trailing_overclocks[:, ::2]
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        odd_detector_bias = odd_detector_bias[:, 4:]
        cols = odd_detector_bias.shape[1]
        odd_detector_bias = odd_detector_bias [:, 0:cols-1]
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
        bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
        bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]

        active_quad[:, ::2] = bias_subtracted_quad_odd
        active_quad[:, 1::2] = bias_subtracted_quad_even
        all_quads.append(active_quad)



    return np.array(all_quads)

def perform_temp_correction(active_quads, fpe_temp):
    """
    Apply the FPE temp sensitivity correction based on exponential model
    The temperature will be scaled to 42C, the temp at which the radiance cal.
    data was acquiredSo the numerator will be fixed foe each octant. The demoninator
    will chage as the FPE temp of the image acquired

    """
    a = [17.250, 17.220, 17.410, 17.450,
         17.210, 17.130, 18.170, 18.110]
    b = [1.042E-04, 1.167E-04, 1.942E-04, 1.059E-04,
         2.011E-04,	2.443E-04,	-2.559E-04, -2.258E-04]
    c = [-4.755E-02, -5.159E-02, -6.725E-03, -7.535E-02,
         -1.330E-01,	-1.132E-01, -9.538E-01,	-9.326E-01]
    d = [-6.404E-02,	-6.422E-02, -1.479E-01,	-5.194E-02,
         -4.342E-02,	-5.010E-02, -1.790E-02,	-1.904E-02]
    numerator = [17.322, 17.301, 17.553,	17.519,
                 17.335, 17.293, 17.526,	17.520]

    all_quads = []
    num_quads = active_quads.shape[0]

    for i in range(0, num_quads):
        active_quad = active_quads[i]
        active_quad_odd = active_quad[:, ::2]
        active_quad_even = active_quad[:, 1::2]
        correction_factor_odd = numerator[2*i] / (a[2*i]*np.exp(b[2*i]*fpe_temp)+
                                                  c[2*i]*np.exp(d[2*i]*fpe_temp))
        correction_factor_even = numerator[2*i+1] / (a[2*i+1]*np.exp(b[2*i+1]*fpe_temp)+
                                                     c[2*i+1]*np.exp(d[2*i+1]*fpe_temp))
        active_quad[:, ::2] = active_quad_odd * correction_factor_odd
        active_quad[:, 1::2] = active_quad_even * correction_factor_even
        all_quads.append(active_quad)

    return np.array(all_quads)

def apply_linearity_correction(active_quad, quad_name):
    """ Ok Read the look up table to perfrom non-linearity correction
     based on the look up table"""

    quad = quad_name.split(' ')[1]
    active_quad[active_quad < 0] = 0 # Replace negative values by 0
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file = 'Linearity_look_up_table_final_DN.csv'
    linearity_file = pd.read_csv(path+'/'+ linearity_file)
    lut_even = linearity_file[['Input Signal (DN)', quad+str(1)]].dropna()
    lut_even['DN'] = lut_even.pop('Input Signal (DN)').astype(int)
    lut_even = lut_even.set_index('DN')
    even_detector_active_quad = active_quad[:, 1::2]
    nx_half_quad, ny_half_quad = even_detector_active_quad.shape
    even_detector_active_quad = np.reshape(np.array(even_detector_active_quad),
                                           (nx_half_quad*ny_half_quad, 1))
    dframe = pd.DataFrame(data=even_detector_active_quad)
    dframe.columns = ['DN']
    datin_even = (dframe['DN']).astype(int)
    linearized_quad_even = datin_even.apply(lambda x: lut_even.ix[x])
    even_pixels_quad = np.reshape(linearized_quad_even.values, (nx_half_quad, ny_half_quad))
    dframe = None
    lut_odd = linearity_file[['Input Signal (DN)', quad + str(2)]].dropna()
    lut_odd['DN'] = lut_odd.pop('Input Signal (DN)').astype(int)
    lut_odd['Value'] = lut_odd.pop(quad+str(2)).astype(float)
    lut_odd = lut_odd.set_index('DN')
    odd_detector_active_quad = active_quad[:, ::2]
    nx_half_quad, ny_half_quad = odd_detector_active_quad.shape
    odd_detector_active_quad = np.reshape(np.array(odd_detector_active_quad),
                                          (nx_half_quad*ny_half_quad, 1))
    dframe = pd.DataFrame(data=odd_detector_active_quad)
    dframe.columns = ['DN']
    datin_odd = (dframe['DN']).astype(int)
    linearized_quad_odd = datin_odd.apply(lambda x: lut_odd.ix[x])
    odd_pixels_quad = np.reshape(linearized_quad_odd.values, (nx_half_quad,
                                                              ny_half_quad))
    # now let's put quad in place
    active_quad[:, 1::2] = even_pixels_quad
    active_quad[:, ::2] = odd_pixels_quad
    return  active_quad

def perform_smear_removal(active_quads, int_time):
    """Perform the SMEAR subtraction using Active Quad Method.
    The underlying assumption in smear subtraction is that the dark current
    in the storage region is really small and hence neglected from the analysis.
    typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor*tRO)
    tft = 8ms
    Note all corrections are applied per octant. Please update this.
    """
    all_quads = []
    #frame_transfer = 8.3333ms
    tft = 8.3333
    num_quads = active_quads.shape[0]
    for quads in range(0, num_quads):

        active_quad = active_quads[quads]
        active_quad_odd = active_quad[:, ::2]
        active_quad_even = active_quad[:, 1::2]
        smear_factor_odd = tft / (int_time + tft)* active_quad_odd.mean(axis=0)
        smear_factor_even = tft / (int_time + tft)* active_quad_even.mean(axis=0)
        smear_subtracted_quad_odd = active_quad_odd - smear_factor_odd[None, :]
        smear_subtracted_quad_even = active_quad_even - smear_factor_even[None, :]
        active_quad[:, ::2] = smear_subtracted_quad_odd
        active_quad[:, 1::2] = smear_subtracted_quad_even
        all_quads.append(active_quad)

    return all_quads


def parse_prnu_file():
    """ Read the PRNU hdf file provided by BATC. It takes the spectrometer
        orientation including the overclocks.
        """
    hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\
                batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
    file = h5py.File(hdf_name, 'r')
    prnu = file.get('prnu')
    prnu = np.array(prnu).transpose()
    quad_d = prnu[2:1030, 10:1034]
    quad_c = prnu[2:1030, 1078:2102]
    quad_a = prnu[1062:2090, 10:1034]
    quad_b = prnu[1062:2090, 1078:2102]
    prnu_map_lower = np.concatenate((quad_d, quad_c), axis=1)
    prnu_map_upper = np.concatenate((quad_a, quad_b), axis=1)
    prnu_map = np.concatenate((prnu_map_lower, prnu_map_upper), axis=0)
    return prnu_map


def read_prnu_files():
    """ The PRNU correction is based on each quad. These were derived using
    very bright and "homogenous" images acquired during the photon transfer
    integration time sweep
    """
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_map_Median'
    prnu_mask_a = np.genfromtxt(file_path +'/' + 'Quad A_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_b = np.genfromtxt(file_path +'/' + 'Quad B_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_c = np.genfromtxt(file_path +'/' + 'Quad C_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask_d = np.genfromtxt(file_path +'/' + 'Quad D_Final_PRNU.csv',
                                delimiter=',')
    prnu_mask = [prnu_mask_a, prnu_mask_b, prnu_mask_c, prnu_mask_d]
    return prnu_mask

def remove_cross_talk(full_frame):
    """ Apply the cross talk correction. The correction factor was provided by
    is 0.0015 for quads within the same CCD and 0 outside the CCD
    """
    quad_a = full_frame[0, :, :] - 0.0015*full_frame[1, :, :]
    quad_b = full_frame[1, :, :] - 0.0016*full_frame[0, :, :]
    quad_c = full_frame[2, :, :] - 0.00144*full_frame[3, :, :]
    quad_d = full_frame[3, :, :] - 0.00143*full_frame[2, :, :]
    corrected_full_frame = [quad_a, quad_b, quad_c, quad_d]
    return corrected_full_frame

def apply_electronics_gain(full_frame, difference):
    """ This function applied the electronics gain on octant basis. However, these
    numbers were derived using FPS testing where TSOC magnitude of even lines were
    greater than odd lines. The even/odd lines read out AD may change.
    Therefore, need to check the magnitude of even and odd TSOC for each image
    to apply the gains correctly.
    """
    electronics_gain_odd = [0.0601, 0.0596, 0.0604, 0.0605]
    electronics_gain_even = [0.0602, 0.0599, 0.0605, 0.0608]
    all_quads = []
    num_quads = full_frame.shape[0]
    for quads in range(0, num_quads):
        active_quad = full_frame[quads, :, :]
        if difference[quads] < 0: # Note: Difference is odd-even
            gain_even = 1/electronics_gain_even[quads]
            gain_odd = 1/electronics_gain_odd[quads]
        elif  difference[quads] > 0:
            gain_even = 1/electronics_gain_odd[quads]
            gain_odd = 1/electronics_gain_even[quads]
        gain_even = 1/electronics_gain_even[quads]
        gain_odd = 1/electronics_gain_odd[quads]
        even_detector_active_quad = gain_even*active_quad[:, ::2]
        odd_detector_active_quad = gain_odd*active_quad[:, 1::2]
        active_quad[:, ::2] = even_detector_active_quad
        active_quad[:, 1::2] = odd_detector_active_quad
        all_quads.append(active_quad)
    return np.array(all_quads)

def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image.

     740nm =================================================================
         |                              |                                |
         |                              |                                |
         |        Amplifier A           |         Amplifier B            |
         |                              |                                |
         |                              |                                |
         |------------------------------|--------------------------------|
         |                              |                                |
         |                              |                                |
         |        Amplifier C           |         Amplifier D            |
         |                              |                                |
         |    *(Index 1,1)              |                                |
   290nm =================================================================
          NORTH EDGE                                         SOUTH EDGE

    """
    quad_a = full_frame[0, :, :]
    quad_b = full_frame[1, :, :]
    quad_c = full_frame[2, :, :]
    quad_d = full_frame[3, :, :]
    uv_ccd = np.concatenate((quad_d, np.fliplr(quad_c)),
                            axis=1)
    visible_ccd = np.concatenate((np.flipud(quad_a), np.rot90(quad_b, 2)),
                                 axis=1)
    processed_image = np.concatenate((uv_ccd, visible_ccd), axis=0)
    return processed_image
