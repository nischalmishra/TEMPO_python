# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:45:14 2017
This code was written by Nischal Mishra to grab the  TEMPO Raw data as delivered by BATC.
For details ragarding ccsds and grddp , please folloe the BATC deliverable.

@author: nmishra
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
#pylint: disable= E1101
from mpl_toolkits.axes_grid1 import make_axes_locatable


def parse_ccsds_header(ccsds_header_buf):
    """
    Check the ccsds headers
    """

    ccsds_header = {}
    ccsds_header['primary'] = {}
    ccsds_header['secondary'] = {}
    #let us populate the dictionary
    ccsds_header['primary']['ccsds_packet_size'] = np.bitwise_or(np.uint32(ccsds_header_buf[4] <<8), np.uint32(ccsds_header_buf[5]))
    ccsds_header['primary']['ccsds_version'] = ccsds_header_buf[0] & int('11100000', 2)
    ccsds_header['primary']['ccsds_type'] = ccsds_header_buf[0] & int('00010000', 2)
    ccsds_header['primary']['second_header_flag'] = ccsds_header_buf[0] & int('00001000', 2)
    ccsds_header['primary']['APID'] = np.bitwise_or(np.uint32(ccsds_header_buf[0]<<8), np.uint32(ccsds_header_buf[1])) & int('0000011111111111', 2)
    ccsds_header['primary']['CCSDS_seq_flag'] = ccsds_header_buf[2] & int('11000000', 2)
    ccsds_header['primary']['ccsds_packet_seq_count'] = np.bitwise_or(np.uint32(ccsds_header_buf[2]<<8), np.uint32(ccsds_header_buf[3])) & int('0011111111111111', 2)
    # now secondary
    ccsds_header['secondary']['img_time_day'] = np.bitwise_or(np.uint32(ccsds_header_buf[6] <<8), np.uint32(ccsds_header_buf[7]))
    ccsds_header['secondary']['img_time_msecs'] = np.uint32(ccsds_header_buf[8]<< 8*3) | np.uint32(ccsds_header_buf[9]<<8*2) | np.uint32(ccsds_header_buf[10]<<8) | np.uint32(ccsds_header_buf[11])
    ccsds_header['secondary']['img_time_usecs'] = np.uint32(ccsds_header_buf[12]<<8) | np.uint32(ccsds_header_buf[13])
    ccsds_header['secondary']['img_sync_word'] = np.uint32(ccsds_header_buf[14]<< 8*3) | np.uint32(ccsds_header_buf[15]<<8*2) | np.uint32(ccsds_header_buf[16]<<8) | np.uint32(ccsds_header_buf[17])
    return ccsds_header

def ccsds_packet_parse(ccsds_packet_buf):
    """ Parse the ccsds packets
    """
    ccsds_header_length = 18 # bytes  (includes primary and secondary)
    #ccsds_secondary_header_length = 12
    CCSDS_HEADER = parse_ccsds_header(ccsds_packet_buf[0:ccsds_header_length])
    science_image_padding = 2 # bytes
    img_bytes = 10560 # ; (2L*2112L*20L)/8L  ; 2 rows * 2112 pix/row * 20 bits/pix  * 1byte/8bits
    crc_bytes = 4
    #zeropad_bytes = 1
    out = {}
    out['pixel_data'] = {}
    out['ccsds_header'] = {}
    i1 = ccsds_header_length + science_image_padding
    i2 = i1 + img_bytes - 1
    i3 = i2 + 1
    i4 = i3 + crc_bytes - 1
    i5 = i4 + 1
    img_data = ccsds_packet_buf[i1:i2+1]
    img_crc = ccsds_packet_buf[i3:i4+1]
    zeropad = ccsds_packet_buf[i5]
    out['ccsds_header'] = CCSDS_HEADER
    out['pixel_data']['img_data'] = img_data
    out['pixel_data']['zeropad'] = zeropad
    out['pixel_data']['img_crc'] = img_crc
    return out


def parse_pixels_from_ccsds_packet(ccsds_packet_data):
    """ Parse pixels from ccsds packets
    """
    pixel_all = []
    npix = int(len(ccsds_packet_data)*8/20/2)
    for i in range(0, npix):
        i1 = i*5
        pixel1 = ((np.int64(ccsds_packet_data[i1])<< 8*2) | (np.int64(ccsds_packet_data[i1+1])<< 8) | np.int64(ccsds_packet_data[i1+2])) >> 4
        pixel_all.append(pixel1)
        pixel2 = (ccsds_packet_data[i1+2] & int('0000000F', 16)) <<8*2 | np.int64(ccsds_packet_data[i1+3]<<8) | np.int64(ccsds_packet_data[i1+4])
        pixel_all.append(pixel2)
    return pixel_all


def make_fpa_from_bit_streams(pixel_data):
    """ Now arrange the pixel streams into quads
    """
    nx_quad = 1056
    ny_quad = 1046
    full_fpa = []
    ncols = len(pixel_data)
    #nrows = len(pixel_data[0])*8/20/2
    for i in range(0, ncols):
        all_pixels = parse_pixels_from_ccsds_packet(pixel_data[i])
        #print(all_pixels)
        full_fpa.append(all_pixels)
    full_fpa = np.array(full_fpa)
    full_fpa = np.reshape(full_fpa, (full_fpa.shape[0]*full_fpa.shape[1], 1))
    j = np.arange(0, nx_quad*ny_quad*4, 4)
    quad_a = np.reshape(full_fpa[j], (ny_quad, nx_quad))
    quad_b = np.reshape(full_fpa[j+1], (ny_quad, nx_quad))
    quad_d = np.reshape(full_fpa[j+3], (ny_quad, nx_quad))
    quad_c = np.resize(full_fpa[j+2], (ny_quad, nx_quad))
    return np.array([quad_a, quad_b, quad_c, quad_d])
 
def create_image(full_frame):
    """ Ok, lets look at image
    """
    lower = np.concatenate((full_frame[3,:,:], np.fliplr(full_frame[2, :, :])),
                           axis=1)
    upper = np.concatenate((np.flipud(full_frame[0, :,:]), 
                            np.rot90(full_frame[1,:,:], 2)), axis=1)
    full_frame = np.concatenate((lower, upper), axis=0)
    fig_ax = plt.gca()    
    image = fig_ax.imshow(full_frame, cmap='jet',
                          origin='lower', interpolation='none') 
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
   # plt.grid(False)   
    #plt.show()
    #cc
    
      

def main(input_image_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Spectral_reg_spatial_disc'):
    """
    Main function
    """
    #input_image_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\2018.06.26'
    image_save_dir = os.path.join(input_image_dir, 'raw_hdf_files')
    if not os.path.exists(image_save_dir):
        os.makedirs(image_save_dir)
    data_path_all = sorted([each for each in os.listdir(input_image_dir)
                            if each.endswith('.img')])

    for data_path in data_path_all:
        input_filename = os.path.join(input_image_dir, data_path)
        file_info = os.stat(input_filename)
        file_size = file_info.st_size
        all_image = []
        img = {}
        img['header'] = {}
        img['pixel_data'] = {}
        num_grddp_packets = -1
        with open(input_filename, mode='rb') as binary_file:
            for chunk in iter(lambda: binary_file.read(21179), b''):
                grddp_header = chunk[0:7]
                grddp_header = list(map(int, grddp_header))
                #print(grddp_header)
                grddp_packet_length = np.int64(grddp_header[4])<<8 | np.int64(grddp_header[5])
                grddp_packet_plus_crc = chunk[8 :]
                grddp_packet_plus_crc = list(map(int, grddp_packet_plus_crc))
                #crc_byte = grddp_packet_plus_crc[-1]
                ccsds_packet_within_grddp = int(grddp_packet_length/2)
                ccsds_packet1 = ccsds_packet_parse(grddp_packet_plus_crc[0:ccsds_packet_within_grddp])
                if num_grddp_packets < 0:
                    num_grddp_packets = int(file_size/(8+1+grddp_packet_length))
                    img['header'] = ccsds_packet1['ccsds_header']
                    img['pixel_data']['img_data'] = ccsds_packet1['pixel_data']['img_data']*num_grddp_packets*2
                    img['pixel_data']['zeropad'] = ccsds_packet1['pixel_data']['zeropad']*num_grddp_packets*2
                    img['pixel_data']['img_crc'] = ccsds_packet1['pixel_data']['img_crc']*num_grddp_packets*2
                    img['pixel_data']['img_data'] = bytes(0)
                all_image.append(ccsds_packet1['pixel_data']['img_data'])
                ccsds_packet1 = ccsds_packet_parse(grddp_packet_plus_crc[ccsds_packet_within_grddp:grddp_packet_length])
                all_image.append(ccsds_packet1['pixel_data']['img_data'])

        #Now lets save them as hdf file
        full_frame = make_fpa_from_bit_streams(all_image)
       
        QUADS = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
        HF_FILE_NAME = os.path.join(image_save_dir, data_path+'.h5')
        num_quads = np.array(full_frame).shape[0]
        with h5py.File(HF_FILE_NAME, 'w') as hf:
            for i in range(0, num_quads):
                file_header = QUADS[i]
                hf.create_dataset(file_header, data=full_frame[i, :, :])

if __name__ == "__main__":
    main()
