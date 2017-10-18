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
from scipy.io.idl import readsav
import scipy.io as sio
import pandas as pd
import h5py


def read_mat_file(filename):
    """ Read .mat files provided by BALL to process Spectrometer Data. May
    need this for quick validation
    """
    mat_contents = sio.loadmat(filename)
    image_data = mat_contents['image']['img'][0, 0]
    return image_data


def read_idl_file(data_file):
    """Read IDL variables. Dave's IDL code to read spectrometer data runs
     much faster than my Python script. Also if there are .sav files,
     I can use them in    future rather than reading the whole file again
     """
    read_variable = readsav(data_file)
    full_frame_image = read_variable.quads
    return full_frame_image


def read_outlier_mask():
    """ The outlier mask is created based on dark data. The idea is to create
    masked array of active quads and not use these in further calculations
    """
    #TO DO : Fix the quad orientation. Also, try and correlate the outliers
    # from FPS test to spectrometer data

    file_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask'
    file_name = 'final_outlier_mask.csv'
    outlier_mask = np.genfromtxt(file_path +'/'+ file_name, delimiter=',')
    rows, cols = outlier_mask.shape
    mask_a = outlier_mask[0: rows/2, 0:cols/2]
    mask_b = outlier_mask[rows/2:rows, 0:cols/2]
    mask_c = outlier_mask[rows/2:rows, cols/2:cols]
    mask_d = outlier_mask[0:rows/2, cols/2:cols]
    outlier_mask = [mask_a, mask_b, mask_c, mask_d]
    return outlier_mask


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

def perform_bias_subtraction(raw_quads):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    all_quads = []
    num_quads, nx_quad, ny_quad = raw_quads.shape
    for quads in range(0, num_quads):
        active_quad = raw_quads[quads, 2:1030, 10:1034]
        trailing_overclocks = raw_quads[quads, 2:1030, 1034:1056]
        spec_pix, spat_pix = active_quad.shape
        bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
        even_detector_bias = trailing_overclocks[ :, ::2]
        # remove outliers
        # First 4 hot lines in even and odd
        # last odd lne in odd
        even_detector_bias = even_detector_bias[:, 4:]
        avg_bias_even = np.mean(even_detector_bias, axis=1)
        odd_detector_bias = trailing_overclocks[:, 1::2]
        odd_samples = odd_detector_bias[:, 4:]
        rows, cols = odd_samples.shape
        odd_detector_bias = odd_samples[:, 0:cols-1]
        avg_bias_odd = np.mean(odd_detector_bias, axis=1)
        even_detector_active_quad = active_quad[:, ::2]
        odd_detector_active_quad = active_quad[:, 1::2]
        bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None]
        bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None]
        bias_subtracted_quad = np.reshape(bias_subtracted_quad, (spec_pix, spat_pix))
        bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
        bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
        all_quads.append(bias_subtracted_quad)

    return np.array(all_quads)

def apply_linearity_correction(active_quad, quad_name, num_coadds):
    """ Ok Read the look up table to perfrom non-linearity correction
     based on the look up table"""
    quad = quad_name.split(' ')[1]
    active_quad[active_quad < 0] = 0 # Replace negative values by 0
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file = 'Linearity_look_up_table_final_DN.csv'
    linearity_file = pd.read_csv(path+'/'+ linearity_file)
    lut_even = linearity_file[['Input Signal (DN)', quad+str(1)]].dropna()
    lut_even['DN'] = lut_even.pop('Input Signal (DN)').astype(int)
#    #sanity check
#    print(lut_even['DN'].head())
#    print(lut_even['Value'].head())
#    cc

    lut_even = lut_even.set_index('DN')
    even_detector_active_quad = active_quad[:, ::2]
    nx_half_quad, ny_half_quad = even_detector_active_quad.shape
    even_detector_active_quad = np.reshape(np.array(even_detector_active_quad),
                                           (nx_half_quad*ny_half_quad, 1))
    dframe = pd.DataFrame(data=even_detector_active_quad)
    dframe.columns = ['DN']
    datin_even = (dframe['DN']/num_coadds).astype(int)
    linearized_quad_even = datin_even.apply(lambda x: lut_even.ix[x])
    even_pixels_quad = np.reshape(linearized_quad_even.values, (nx_half_quad, ny_half_quad))
    dframe = None
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
    datin_odd = (dframe['DN']/num_coadds).astype(int)
    linearized_quad_odd = datin_odd.apply(lambda x: lut_odd.ix[x])
    odd_pixels_quad = np.reshape(linearized_quad_odd.values, (nx_half_quad,
                                                              ny_half_quad))
    # now let's put quad in place
    nx_quad, ny_quad = active_quad.shape
    linearized_quad = np.array([[0]*ny_quad]*nx_quad)
    linearized_quad = np.reshape(linearized_quad, (nx_quad, ny_quad))
    linearized_quad[:, ::2] = even_pixels_quad
    linearized_quad[:, 1::2] = odd_pixels_quad
    return linearized_quad

def perform_smear_removal(active_quad, int_time):
    """Perform the SMEAR subtraction using Active Quad Method.
    The underlying assumption in smear subtraction is that the dark current
    in the storage region is really small and hence neglected from the analysis.
    typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor*tRO)
    tft = 8ms
    """
    all_quads = []
    frame_transfer = 8
    num_quads, nx_quad, ny_quad = active_quad.shape
    for quads in range(0, num_quads):
        smear_factor = (frame_transfer / (int_time + frame_transfer))* \
                       np.mean(active_quad[quads], axis=0)
        smear_subtracted_quad = active_quad - smear_factor[None, :]
        all_quads.append(smear_subtracted_quad)
    return smear_subtracted_quad

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
        #ccd_temp = group2.attrs['ccd_temp']
        #fpe_temp = group2.attrs['fpe_temp']
    # sanity check
    #print(coadds, int_time, ccd_temp, fpe_temp)

    return coadds, int_time


def parse_prnu_file():
    """ Read the PRNU hdf file provided by BATC. It takes the spectrometer
        orientation including the overclocks.
        """
    hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
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


def apply_prnu_correction(full_frame):
    """ The PRNU correction is based on each quad. These were derived using
    very bright and "homogenous" images acquired during the photon transfer
    integration time sweep
    """
    #prnu_all_quads = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\Final_Map\PRNU_mask_active_regions.csv'
    prnu_mask =np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\spectro_prnu_interpol.csv', delimiter=',')
    #prnu_mask = np.genfromtxt(prnu_all_quads, delimiter=',')
    full_frame = full_frame/prnu_mask
    return full_frame

def remove_cross_talk(full_frame):
    """ Apply the cross talk correction. The correction factor was provided by
    is 0.0015 for quads within the same CCD and 0 outside the CCD

    """
    quad_a = full_frame[0, :, :] - 0.0015*full_frame[1, :, :]
    quad_b = full_frame[1, :, :] - 0.0015*full_frame[0, :, :]
    quad_c = full_frame[2, :, :] - 0.0015*full_frame[3, :, :]
    quad_d = full_frame[3, :, :] - 0.0015*full_frame[2, :, :]
    corrected_full_frame = [quad_a, quad_b, quad_c, quad_d]
    return corrected_full_frame

def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image

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
