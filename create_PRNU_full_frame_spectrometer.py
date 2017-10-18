# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 10:38:29 2017

@author: nmishra
"""
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import h5py



def parse_PRNU_file():
    """ Read the PRNU hdf file provided by BATC. It takes the spectrometer 
        orientation including the overclocks.
        """
    hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
    file = h5py.File(hdf_name, 'r')
    prnu = file.get('prnu')
    prnu = np.array(prnu).transpose()
    quad_D = prnu[2:1030, 10:1034]
    quad_C = prnu[2:1030, 1078:2102]
    quad_A = prnu[1062:2090, 10:1034]
    quad_B = prnu[1062:2090, 1078:2102]
    prnu_map_lower = np.concatenate((quad_D, quad_C), axis=1)
    prnu_map_upper = np.concatenate((quad_A, quad_B), axis=1)
    prnu_map = np.concatenate((prnu_map_lower, prnu_map_upper), axis=0)
    return prnu_map

def create_final_image(full_frame):
    """ Arrange the quads to create the final TEMPO Image
    
    """
    quad_A = full_frame[0, :, :] 
    quad_B = full_frame[1, :, :] 
    quad_C = full_frame[2, :, :] 
    quad_D = full_frame[3, :, :] 
    UV_CCD = np.concatenate((quad_D, np.fliplr(quad_C)), 
                                axis=1)
    Visible_CCD = np.concatenate((np.flipud(quad_A), np.rot90(quad_B, 2)),
                                     axis=1)
    processed_image = np.concatenate((UV_CCD, Visible_CCD), axis=0)
    
    return processed_image

def create_spectrometer_image(image):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image, cmap='bwr', origin='lower')
    plt.title('TEMPO PRNU Mask')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
   

def main():
    
    prnu_quad_A = genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Analysis_second_order\Quad A\Final_PRNU\Filtered_PRNU\Quad A_Final_PRNU.csv', delimiter=',')
    prnu_quad_B = genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Analysis_second_order\Quad B\Final_PRNU\Filtered_PRNU\Quad B_Final_PRNU.csv', delimiter=',')
    prnu_quad_C = genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Analysis_second_order\Quad C\Final_PRNU\Filtered_PRNU\Quad C_Final_PRNU.csv', delimiter=',')
    prnu_quad_D = genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Analysis_second_order\Quad D\Final_PRNU\Filtered_PRNU\Quad D_Final_PRNU.csv', delimiter=',')
   
    # Note than the final image shape in tempo spectrometer is 1028 by 1024.
    #PRNU started from 100:900 in spectral direction
    # In order to remove vignetting, the egdge pixels were not used in the PRNU
    # So we will extrapolate the values to create a full size PRNU

    prnu = [prnu_quad_A, prnu_quad_B, prnu_quad_C, prnu_quad_D]
    processed_image = create_final_image(np.array(prnu)) 
    create_spectrometer_image(processed_image)    
    np.savetxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\spectro_prnu_interpol.csv', np.array(processed_image), delimiter=",")
    
    prnu_quad = genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\spectro_prnu_interpol.csv', delimiter=',')
    print(prnu_quad.shape)
    cc
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(prnu_quad, cmap='nipy_spectral', origin='lower')
    plt.title('TEMPO PRNU Mask')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.show()
    
#    quads_name = ['QuadA','QuadB','QuadC','QuadD',]
#    
#        # find the 
#    rows, cols = prnu_quad_A.shape
#    print(rows, cols)
#   
#   
#    x1_rows = list(range(100, 900))    
#    new_points_end = list(range(900, 1028)) 
#    new_point_begin = list(range(0, 100))
#    
#   
#    prnu_all = []
#    for quads in range(0, 4):
#        y2_all_end = []
#        y2_all_begin = []
#        measured_val = prnu[quads]
#        for i in range(cols):
#            
#           fit = np.polyfit(x1_rows, measured_val[:, i], 1)
#           #print(len(measured_val[:,i]))
#           #cc
#           line= np.poly1d(fit)
#           desired_val_end = line(new_points_end)
#           desired_val_begin = line(new_point_begin)
#           y2_all_end.append(desired_val_end)
#           y2_all_begin.append(desired_val_begin)
#    #print(np.array(y2_all).shape)
#        y2_all_end = np.array(y2_all_end).T
#        y2_all_begin = np.array(y2_all_begin).T   
#        final_mask = np.zeros((1028, 1024))
#        final_mask[0:900-rows, :] = np.flipud(measured_val[0:900-rows, :])
#        final_mask[100:rows+100, :] = measured_val
#        final_mask[rows+100:1028, :] = np.flipud(measured_val[672:800, :])
#        print(final_mask.shape)
#        
#        #plt.show()
#        prnu_all.append(final_mask)
#    print(np.array(prnu_all).shape)
#    processed_image = create_final_image(np.array(prnu_all))
#    plt.figure()
#    ax = plt.gca()
#    image = ax.imshow(processed_image, cmap='nipy_spectral', origin='lower')
#    plt.title('TEMPO PRNU Mask, '+ quads_name[quads])
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("right", size="5%", pad=0.05)
#    plt.colorbar(image, cax=cax)
#    plt.grid(False)
#    plt.show()
#    
#    np.savetxt(r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\spectro_prnu_interpol.csv', np.array(processed_image), delimiter=",")
#    
        
        
    
    
    
if __name__ == "__main__":
    main()
    
    
    
    
    