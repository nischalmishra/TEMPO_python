# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:45:14 2017

@author: nmishra
"""

import os
import numpy as np
import matplotlib.pyplot as plt

#grddp_header = bytearray(8)
#ccsds_header = bytearray(18)

def parse_ccsds_header(ccsds_header_buf):
    #break-out ccsds_header items and return structure based on C&T HANDBOOK, version ?
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
    pixel_all = []
    npix = int(len(ccsds_packet_data)*8/20/2)
    for i in range(0, npix):
        i1 = i*5
        pixel1 = (np.int32(ccsds_packet_data[i1]<< 8*2) | np.int32(ccsds_packet_data[i1+1]<< 8) | np.uint32(ccsds_packet_data[i1+2]))>> 4
        pixel_all.append(pixel1)
        pixel2 = np.int32(ccsds_packet_data[i1+2] & int('0000000F', 16) <<16) | np.int32(ccsds_packet_data[i1+3]<<8) |np.uint32(ccsds_packet_data[i1+4])
        pixel_all.append(pixel2)
    return pixel_all




def make_fpa_from_bit_streams(pixel_data):
    nx_quad = 1056
    ny_quad = 1046
    full_fpa = []
    ncols = len(pixel_data)
    #nrows = len(pixel_data[0])*8/20/2
    for i in range(0, ncols):
        all_pixels = parse_pixels_from_ccsds_packet(pixel_data[i])
        full_fpa.append(all_pixels)
    full_fpa = np.array(full_fpa)
    full_fpa = np.reshape(full_fpa, (full_fpa.shape[0]*full_fpa.shape[1], 1))
    j = np.arange(0, nx_quad*ny_quad*4, 4)
    quad_a = np.reshape(full_fpa[j], (ny_quad, nx_quad))
    quad_b = np.reshape(full_fpa[j+1], (ny_quad, nx_quad))
    quad_d = np.reshape(full_fpa[j+3], (ny_quad, nx_quad))
    quad_c = np.resize(full_fpa[j+2], (ny_quad, nx_quad))
    lower = np.concatenate((quad_d, np.fliplr(quad_c)), axis=1)
    upper = np.concatenate((np.flipud(quad_a), np.rot90(quad_b, 2)), axis=1)
    full_frame = np.concatenate((lower, upper), axis=0)
    return full_frame

def main():
    input_filename = r'C:\Users\nmishra\Workspace\TEMPO\Data\Spectrometer\2017.06.28\2017_06_28_18_15_59_235208.img'
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
          grddp_packet_length = np.bitwise_or(np.uint32(grddp_header[4]<<8), np.uint32(grddp_header[5]))
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

    full_frame = make_fpa_from_bit_streams(all_image)
    plt.figure()
    plt.imshow(full_frame, origin='lower', cmap='nipy_spectral', aspect='auto')
    plt.grid(False)
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main()




