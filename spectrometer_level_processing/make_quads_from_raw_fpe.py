# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 10:59:09 2017

@author: nmishra
This function contains the list of functions to organize various quads.
Includes function to organize the 4 quads from binary files, stack 4 quads on
top of each other, seperate out the overclock pixels, active regions, buffer
regions etc. The dimensions for these reconstruction have been taken directly
from Ball document  @ https://nx.larc.nasa.gov/dsweb/View/Collection-76892.

Please note that the row and column indices used in the python script look
different than in the document because of how python indexes rows and columns.
For instance python index starts from 0 (same as in IDL). Ball scripts provided
at the end of the document uses matlab script and matlab index starts at one.

    Example:
    a = np.arange(1,10)
    array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    a[0:5]
    [0:5]
    array([1, 2, 3, 4, 5])

"""

import os
import numpy as np
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from astropy.io import fits

#******************************************************************************

def make_quads_from_fits_file(file_path, data_path):
    data_path_name_split = data_path.split('_')
    int_time = data_path_name_split[4].split('-')[1]
    temp = data_path_name_split[-2].split('-')[1]
    frames = data_path
    data_file = os.path.join(file_path, data_path)
    hdulist = fits.open(data_file)
    full_frame = hdulist[0].data
    return full_frame, int_time, temp, frames

def make_quads_from_IDL_sav_file(file_path, data_path, names_to_check):
    """ This function takes Dave Flittner's IDL variables are forms quads with
    it
    """
    data_path_name_split = data_path.split('_')
    #print(data_path_name_split) 
    #cc
    
    integration_type = [x for x in names_to_check if x in data_path_name_split]
    collection_type = integration_type[0]
    frames_int = data_path_name_split[-1]
    frames = data_path
    integration_time = round(int(frames_int.split('.')[0])) #for integ_sweep
    #integration_time = round(int(data_path_name_split[0])) # for intensity sweep
    #print(integration_time)
    #cc
    #temp = data_path_name_split[5]
    #print(temp)
    #cc
    
    
    
    #frames = data_path
    #print('Integration Time->', integration_time)
    data_file = os.path.join(file_path, data_path)
    #print(data_file)   
    #cc
    IDL_variable = readsav(data_file)    
    #full_frame = IDL_variable.quads
    full_frame = IDL_variable.q    
    #return full_frame, integration_time
    return full_frame, collection_type, frames,integration_time
 



def make_full_frame_from_raw_fpe(raw_bits, nx_quad, ny_quad):
    """ This function makes the full frame image including the overclocks from
    the raw binary data.

    ******Inputs***********
    1. Raw binary stream
    2. nx_quad : num of pixels in x-direction
    3. ny_quad : num of pixels in y-direction

    *******Outputs********
    full_frame[x,y]
    full frame glued together such that the spatial direction, or field angle
    direction, or position along the slit is the column #, or x direction
    AND the dispersion direction or spectral direction is the row #, or the
    y direction. This is like Fig. 5-2 of BATC SER2450865, TEMPO Processing
    Algorithm for Image Reconstruction

    Note : IDL displays columns and rows of matrix. Python diplays rows and
    colums. This function has been made longer purposely  to ensure that I
    understand how quads A, B, C and D are organized in a full frame image.

    TO DO : make the code shorter and optimize the numpy array to make it run
    faster. Not an immediate priority though.

    """
    j = np.arange(0, nx_quad*ny_quad*4, 4)
    quad_a = np.reshape(raw_bits[j], (ny_quad, nx_quad))
    quad_b = np.fliplr(np.reshape(raw_bits[j+1], (ny_quad, nx_quad)))
    quad_d = np.flipud(np.reshape(raw_bits[j+3], (ny_quad, nx_quad)))
    quad_c = np.rot90(np.resize(raw_bits[j+2], (ny_quad, nx_quad)), 2)
    lower = np.concatenate((quad_a, quad_b), axis=1)
    upper = np.concatenate((quad_d, quad_c), axis=1)
    full_frame = np.concatenate((lower, upper), axis=0)
   # print('Size of full frame is ', full_frame.shape)
    return full_frame

#*****************************************************************************
def full_frame_to_quads(full_frame, nx_quad, ny_quad):
    """
    This function basically does the same stuff as make_full_frame_from_raw_fpe
    function but organizes the quads differently. This may be handy when each
    quad needs to be analyzed individualy. The orientation of all the quads
    have been made consistent.
    """
    quad_a = full_frame[0:ny_quad, 0:nx_quad]
    quad_b = np.fliplr(full_frame[0:ny_quad, nx_quad:])
    quad_c = np.rot90(full_frame[ny_quad:, nx_quad:], 2)
    quad_d = np.flipud(full_frame[ny_quad:, 0:nx_quad])
    quads = [quad_a, quad_b, quad_c, quad_d]
    return quads

#*****************************************************************************
def find_active_region_each_quads(quads):
    """
    This function grabs the active regions within each quads and puts them in
    array. For rows and columns of active areas in each quads, refer to
    Table 5.1 and 5.2 of TEMPO Processing Algorithm for Image Reconstruction.
    The quads have been oriented, for convinience purpose, such that the starting
    and rows and columns for each quads are same.
    Please note that the python array indexing starts at zero. For example
    A = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    A[0,1:] prints array([2, 3])
    A[1:3, 1:2] prints array([[5],
                        [8]])
    Input : the 4 quads arrange in an array
    Output : the active areas of 4 quads arranged in an array
    """
    active_area_quads = []
    for k in list(range(0, 4)):
        quadi = np.asarray(quads[k])
        active_area = quadi[4:1028, 10:1034]
        active_area_quads.append(active_area)
    active_area_quads = np.asarray(active_area_quads)
    #print('The size of active region is ', active_area_quads.shape)
    return active_area_quads

#*****************************************************************************
def find_trailing_overclock_pixels(quads):
    """This function grabs the trailing pixels within each quads & puts them in
    array. For rows and columns of trailing pixels in each quads, refer to
    Table 5.1 and 5.2 of TEMPO Processing Algorithm for Image Reconstruction.
    The quads have been oriented for convinience purpose such that the starting
    and rows and columns for each quad are same.
    Please note that the python array indexing starts at zero. For example
    A = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    A[0,1:] prints array([2, 3])
    A[1:3, 1:2] prints array([[5],
                             [8]])
    """
    trailing_overclock_quads = []
    for k in list(range(0, 4)):
        quadi = np.asarray(quads[k])
        trailing_area = quadi[0:1046, 1034:1056]
        trailing_overclock_quads.append(trailing_area)
    trailing_overclock_quads = np.asarray(trailing_overclock_quads)
    print('The size of trailing overclock in each quad is ',
          trailing_overclock_quads.shape)
    return trailing_overclock_quads

#*****************************************************************************
def find_leading_overclock_pixels(quads):
    """This function grabs the trailing pixels within each quads & puts them in
    array. For rows and columns of trailing pixels in each quads, refer to
    Table 5.1 and 5.2 of TEMPO Processing Algorithm for Image Reconstruction.
    The quads have been oriented for convinience purpose such that the starting
    and rows and columns for each quad are same.
    Please note that the python array indexing starts at zero. For example
    A = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    A[0,1:] prints array([2, 3])
    A[1:3, 1:2] prints array([[5],
                        [8]])
    """
    leading_buffer_quads = []
    for k in list(range(0, 4)):
        quadi = np.asarray(quads[k])
        leading_buffer_area = quadi[0:1046, 0:10]
        leading_buffer_quads.append(leading_buffer_area)
    leading_buffer_quads = np.asarray(leading_buffer_quads)
    print('The shape of leading buffer in each quads is ',
          leading_buffer_quads.shape)
    return leading_buffer_quads
#*****************************************************************************

def find_smear_pixels(quads):
    """This function grabs the trailing pixels within each quads & puts them in
    array. For rows and columns of trailing pixels in each quads, refer to
    Table 5.1 and 5.2 of TEMPO Processing Algorithm for Image Reconstruction.
    The quads have been oriented for convinience purpose such that the starting
    and rows and columns for each quad are same.
    Please note that the python array indexing starts at zero. For example
    A = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    A[0, 1:] prints array([2, 3])
    A[1:3, 1:2] prints array([[5],
                        [8]])
    """
    smear_pixel_quads = []
    for k in list(range(0, 4)):
        quadi = np.asarray(quads[k])
        smear_pixel_area = quadi[1028:1046, 10:1034]
        smear_pixel_quads.append(smear_pixel_area)
    smear_pixel_quads = np.asarray(smear_pixel_quads)
    print(' The shape of the smear pixel in each quads is ',
          smear_pixel_quads.shape)
    return smear_pixel_quads
#*****************************************************************************

def bias_corrected_quads(active_region, trailing_overclocks):
    """ Performs the bias correction using trailing overclock pixels.
    We still need to decide whether it's the correct approach or not.
    """
    bias_subtracted_quads = []
    for k in list(range(0, 4)):
        active_quadi = np.asarray(active_region[k])
        bias = np.mean(trailing_overclocks[k], axis=1)
        bias = np.reshape(bias, (1046, 1))
        bias_subtract = active_quadi - bias
        bias_subtracted_quads.append(bias_subtract)
    bias_subtracted_quads = np.asarray(bias_subtracted_quads)
    print('the size of bias corrected quad is ', bias_subtracted_quads.shape)
    return bias_subtracted_quads
