# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 10:34:21 2017

@author: nmishra

 This function contains all the routines (in order) required to process the
TEMPO RAW data.Each of the image processing steps is optional  meaning the main
program has the switch option whether to call them or not. While the code has
some redundancy to it priority has been given to code it in such a way that anyone
can understand what is going on in each stages.This code assumes that the raw
image files are stored in hdf format in the order of read outs.
quadA = image[0, :, :]
quadB = image[1, :, :]
quadC = image[2, :, :]
quadD = image[3, :, :]
Although, most of the image processing is mosstly based off octants, each image
processing steps the inputs are in the form of quadrants.
"""

# TO DO : Add the gain conversion, conversion to rates, include FPA temp correction

import os
import numpy as np
#import numpy.ma as ma
from scipy.io.idl import readsav
import h5py
#import matplotlib.pyplot as plt
#pylint: disable= E1101

def read_h5_file(data_file):
    """
    This code assumes that input files are hdf files where data are
    saved as QuadA,quadB, quadC, quadD
    """
   # quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
    file = h5py.File(data_file, 'r')
    quada = file.get('Quad A')
    quadb = file.get('Quad B')
    quadc = file.get('Quad C')
    quadd = file.get('Quad D')
    all_quads = [quada, quadb, quadc, quadd]
    return np.array(all_quads)

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
    meter image. Values returned are integration time, coadds and temperatures
    of electronics and detector array (FPE and FPA)
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
        return coadds, int_time, fpe_temp, ccd_temp

def perform_coaddition_correction(file_name, num_coadds):
    """ This module dooes the coaddition correction. This code reads the coaddition
    info from the image telemetry files that BATC had delivered along with the
    raw image files"""
    corrected_data = file_name/num_coadds
    return corrected_data

def identify_ping_pong_fps():
    """ This function takes mean trailing overclock values from even and odd lines
    and indetifies whoch one has larger magnitude. This information is passed
    along to the main function to check if the ping and pong swtiched states.
    """
    mean_tsoc_fps = np.array([[811.27, 834.42], [809.5, 838.16],
                              [809.22, 835.25], [810.63, 832.7]])
    mean_tsoc_diff_fps = mean_tsoc_fps[:, 0] - mean_tsoc_fps[:, 1]
    phase_tsoc_fps = [0 if i < 0 else 1 for i in mean_tsoc_diff_fps]
    return phase_tsoc_fps

def identify_ping_pong_image(trailing_overclocks):
    """ This function calculates the median offset of each quads based on trailing
    overclock pixels. This is to check if the ping/pong swtiched states from FPS
    level testing. Any  correction algorithm that requires a database
    item needs this check. Examples are non-linearity corrections and
    temperature related corrections and gain
    """
    num_quads = trailing_overclocks.shape[0]
    avg_offset_tsoc = []
    for quads in range(0, num_quads):
        odd_detector_bias = trailing_overclocks[quads, :, ::2]
        even_detector_bias = trailing_overclocks[quads, :, 1::2]
        # remove outliers first 4 lines and last lines
        avg_bias = np.mean(odd_detector_bias[:, 4:-1]) - np.mean(even_detector_bias[:, 4:-1])
        avg_offset_tsoc.append(avg_bias)
    phase_tsoc_image = [0 if i < 0 else 1 for i in avg_offset_tsoc]
    return phase_tsoc_image

def perform_bias_subtraction(active_quads, trailing_overclocks):
    """ Remove offset from active quads.For each octants, toss the outliers
        average the trailing overclock pixels and subtract the offset on a row-by
        row basis.
        """
    num_quads = active_quads.shape[0]
    #active_quads = np.array([[[0.0]*1024]*1028]*4)
    for quads in range(0, num_quads):
        # Break the TSOC into octants
        odd_detector_bias = trailing_overclocks[quads, :, ::2]
        even_detector_bias = trailing_overclocks[quads, :, 1::2]
        avg_bias_odd = np.mean(odd_detector_bias[:, 4:-1], axis=1)
        avg_bias_even = np.mean(even_detector_bias[:, 4:-1], axis=1)
        #active_quad_odd = active_quads[quads, :, ::2]
        #active_quad_even = active
        bias_subtracted_quad_odd = active_quads[quads, :, ::2] - avg_bias_odd[:, None]
        bias_subtracted_quad_even = active_quads[quads, :, 1::2] - avg_bias_even[:, None]
        active_quads[quads, :, ::2] = bias_subtracted_quad_odd
        active_quads[quads, :, 1::2] = bias_subtracted_quad_even
        #print(active_quads.shape)
        #cc
    return active_quads

def apply_non_linearity_correction(active_quads, ping_pong_phase_image,
                                   ping_pong_phase_fps):
    """ Ok Read the look up table to perfrom non-linearity correction
     based on the look up table"""

    base_dir = r'C:\Users\nmishra\Workspace\TEMPO'
    file_dir = r'linearity_testing\Ping_pong_included\Integration_sweep_latest'
    linearity_file_name = os.path.join(base_dir, file_dir,
                                       'Linearity_look_up_table_final_DN.csv')
    linearity_file = np.genfromtxt(linearity_file_name, delimiter=",",
                                   skip_header=0, usecols=list(range(1, 9)))
    num_quads, nx_rows, ny_cols = active_quads.shape
    for i in range(0, num_quads):
        if  ping_pong_phase_image == ping_pong_phase_fps:# check the phase of ping pong
            linearity_file_odd = linearity_file[:, 2*i]
            linearity_file_even = linearity_file[:, 2*i+1]
        else:
            linearity_file_odd = linearity_file[:, 2*i+1]
            linearity_file_even = linearity_file[:, 2*i]
        odd_detector_active_quad = np.reshape(active_quads[i, :, ::2].astype(int),
                                              [nx_rows*int(ny_cols/2), 1])
        odd_detector_active_quad = [linearity_file_odd[a] for a in odd_detector_active_quad]
        even_detector_active_quad = np.reshape(active_quads[i, :, 1::2],
                                               [nx_rows*int(ny_cols/2), 1]).astype(int)
        even_detector_active_quad = [linearity_file_even[a] for a in even_detector_active_quad]
        active_quads[i, :, ::2] = np.reshape(odd_detector_active_quad, [nx_rows, int(ny_cols/2)])
        active_quads[i, :, 1::2] = np.reshape(even_detector_active_quad, [nx_rows, int(ny_cols/2)])
    return active_quads

def perform_temp_correction(active_quads, fpe_temp, ping_pong_phase_image,
                            ping_pong_phase_fps):
    """
    Apply the FPE temp sensitivity correction based on exponential model
    The temperature will be scaled to 42C, the temp at which the radiance cal.
    data was acquired. So the numerator will be fixed foe each octant. The demoninator
    will chage as the FPE temp of the image acquired

   NOTE!! UPDATE THESE COEFFICIENTS PLEASE
    """
    a_val = [17.250, 17.220, 17.410, 17.450,
             17.210, 17.130, 18.170, 18.110]

    b_val = [1.042E-04, 1.167E-04, 1.942E-04, 1.059E-04,
             2.011E-04,	2.443E-04,	-2.559E-04, -2.258E-04]

    c_val = [-4.755E-02, -5.159E-02, -6.725E-03, -7.535E-02,
             -1.330E-01, -1.132E-01, -9.538E-01, -9.326E-01]

    d_val = [-6.404E-02, -6.422E-02, -1.479E-01, -5.194E-02,
             -4.342E-02, -5.010E-02, -1.790E-02, -1.904E-02]

    numerator = [17.322, 17.301, 17.553, 17.519,
                 17.335, 17.293, 17.526, 17.520]
    num_quads = active_quads.shape[0]
    for i in range(0, num_quads):
        correction_factor_odd = numerator[2*i] / (a_val[2*i]*np.exp(b_val[2*i]*fpe_temp)+
                                                  c_val[2*i]*np.exp(d_val[2*i]*fpe_temp))

        correction_factor_even = numerator[2*i+1] / (a_val[2*i+1]*np.exp(b_val[2*i+1]*fpe_temp)+
                                                     c_val[2*i+1]*np.exp(d_val[2*i+1]*fpe_temp))
        if ping_pong_phase_image == ping_pong_phase_fps:#check the phase of ping and pong
            active_quads[i, :, ::2] = active_quads[i, :, ::2] * correction_factor_odd
            active_quads[i, :, 1::2] = active_quads[i, :, 1::2] * correction_factor_even
        else:
            active_quads[i, :, ::2] = active_quads[i, :, ::2] * correction_factor_even
            active_quads[i, :, 1::2] = active_quads[i, :, 1::2] * correction_factor_odd
    return active_quads

def remove_cross_talk(full_frame):
    """ Apply the cross talk correction. The correction factor was provided by
    is 0.0015 for quads within the same CCD and 0 outside the CCD. LaRC estimates
    were slightly different.
    """
    quad_a = full_frame[0, :, :] - 0.0015 * full_frame[1, :, :]
    quad_b = full_frame[1, :, :] - 0.0016 * full_frame[0, :, :]
    quad_c = full_frame[2, :, :] - 0.00144 * full_frame[3, :, :]
    quad_d = full_frame[3, :, :] - 0.00143 * full_frame[2, :, :]
    all_quads = [quad_a, quad_b, quad_c, quad_d]
    return all_quads

def perform_smear_removal(active_quads, int_time):
    """Perform the SMEAR subtraction using Active Quad Method.
    Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor*tRO)
    The underlying assumption in smear subtraction is that the dark current
    in the storage region is really small and hence not used in the smear calculation.
    typically,  tft = 8.333 ms. Note all corrections are applied per octant.
    """
    tft = 8.3333 # frame transfer time
    num_quads = active_quads.shape[0]
   # processed_quads = np.array([[[0.0]*1024]*1028]*4)
    for quads in range(0, num_quads):
        active_quad = active_quads[quads, :, :]
        active_quad_odd = active_quad[:, ::2]
        active_quad_even = active_quad[:, 1::2]
        smear_factor_odd = tft / (int_time + tft)* active_quad_odd.mean(axis=0)
        smear_factor_even = tft / (int_time + tft)* active_quad_even.mean(axis=0)
        smear_subtracted_quad_odd = active_quad_odd - smear_factor_odd[None, :]
        smear_subtracted_quad_even = active_quad_even - smear_factor_even[None, :]
        active_quads[quads, :, ::2] = smear_subtracted_quad_odd
        active_quads[quads, :, 1::2] = smear_subtracted_quad_even
    return active_quads


def perform_gain_correction(active_quads, ping_pong_phase_image,
                            ping_pong_phase_fps):
    """ Apply electronics gain to converts into electrons units.
    """
    larc_gain =	[16.859, 16.766, 17.027, 16.947,
                 17.038, 16.957, 16.861, 16.785]
    num_quads = active_quads.shape[0]
    for i in range(0, num_quads):
        if ping_pong_phase_image == ping_pong_phase_fps:#check the phase of ping and pong
            active_quads[i, :, ::2] = active_quads[i, :, ::2] *larc_gain[2*i]
            active_quads[i, :, 1::2] = active_quads[i, :, 1::2] *larc_gain[2*i+1]
        else:
            active_quads[i, :, ::2] = active_quads[i, :, ::2] *larc_gain[2*i+1]
            active_quads[i, :, 1::2] = active_quads[i, :, 1::2] *larc_gain[2*i]
    return active_quads


def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image.
     290nm =================================================================
         |                              |                                |
         |                              |                                |
         |        Amplifier D           |         Amplifier C            |
         |                              |                                |
         |                              |                                |
         |------------------------------|--------------------------------|
         |                              |                                |
         |                              |                                |
         |        Amplifier A           |         Amplifier B            |
         |                              |                                |
         |    *(Index 1,1)              |                                |
   740nm =================================================================
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
    return np.flipud(processed_image)
