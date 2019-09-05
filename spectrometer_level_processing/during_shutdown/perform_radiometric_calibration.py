# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:34:45 2017

@author: nmishra
    This function creates TEMPO Radiance calibration coefficients. Steps includes
    1. Read the sphere radiance for different source positions
    2. Read the TEMPO pixel to wavelen map.
    3. Perform Spectral Interpolation of Sphere Radiance and TEMPO Pixel to wavelen map.
        I.e Create a pixel to sphere radiance Map
    4. Read the processed TEMPO DSS Data for UV and Visible Channels.
    5. Scale the DSS Data by the Sphere Radiance Map derive in 3.
    6. Compare the BATC derived Radiance Cal. to LaRC derived Radiance Cal.
    """

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d

def create_radiance_plots(rad_cal_dir, radiance):
    """
    This function is used to create plots for sphere radiances
    """
    dframe = pd.read_csv(os.path.join(rad_cal_dir, radiance), delimiter=",")
    sphere_wavelen = dframe['Wavelength'].values
    sphere_radiances = np.array([dframe['Radiance_UV_UG11+Y'].values,
                                 dframe['Radiance_UV_UG11-Y'].values,
                                 dframe['Radiance_UV_UG11_center'].values,
                                 dframe['Radiance_VISIR_NOfilt+Y'].values,
                                 dframe['Radiance_VISIR_NOfilt-Y'].values,
                                 dframe['Radiance_VISIR_NOfilt_Center'].values])

    plt.figure(1)
    plt.plot(sphere_wavelen, sphere_radiances[0, :], 'r', label='Radiance_UV_UG11+Y')
    plt.plot(sphere_wavelen, sphere_radiances[1, :], 'b', label='Radiance_UV_UG11-Y')
    plt.plot(sphere_wavelen, sphere_radiances[2, :], 'g', label='Radiance_UV_UG11_center')
    plt.plot(sphere_wavelen, sphere_radiances[3, :], 'm', label='Radiance_VISIR_NOfilt+Y')
    plt.plot(sphere_wavelen, sphere_radiances[4, :], 'k', label='Radiance_VISIR_NOfilt-Y')
    plt.plot(sphere_wavelen, sphere_radiances[5, :], 'y', label='Radiance_VISIR_NOfilt_center')
    plt.xlim(290, 740)
    plt.grid(linestyle=':')
    plt.ylim(0, 40)
    plt.legend(loc='best')
    plt.title('Measured Radiances for Different Source Positions')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Radiance(uW/cm^2/nm/SR)')
    plt.show()

    plt.figure(2)
    plt.plot(dframe['Wavelength'], dframe['Radiance+Y'], 'r')
    plt.plot(dframe['Wavelength'], dframe['Radiance-Y'], 'b')
    plt.plot(dframe['Wavelength'], dframe['Radiance_center'], 'g')
    plt.xlim(290, 740)
    plt.grid(linestyle=':')
    plt.ylim(0, 40)
    plt.legend(loc='best')
    plt.title('Measured Radiances for Different Source Positions (After UG11 Fitering)')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Radiance(uW/cm^2/nm/SR)')
    plt.show()

    sphere_radiances = np.array([dframe['Radiance+Y'], dframe['Radiance-Y'],
                                 dframe['Radiance_center']])
    unct_est = 100*np.std(sphere_radiances, axis=0, ddof=0)/np.mean(sphere_radiances, axis=0)
    plt.plot(dframe['Wavelength'], unct_est, 'k')
    plt.xlim(290, 740)
    plt.grid(linestyle=':')
    plt.ylim(-0.01, 2)
    plt.title('Variation in Radiance for Different Sphere Positions')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Coffient of Variation [100*Std. Dev/Mean] ')
    plt.show()


def create_sph_rad_to_wavelen_map(rad_cal_dir, radiance, pixel_to_wavelen_map,
                                  radiance_save_dir):

    """ This function performs the spectral interpolation of sphere radiance
    to TEMPO Pixel to wavelen map. Sphere position are defined at +Y,-Y and
    center for borth UV and VIS channelss.
    """

    pixel_to_wavelen_map = np.genfromtxt(os.path.join(rad_cal_dir, pixel_to_wavelen_map),
                                         delimiter=',')
    pixel_to_wavelen_map = np.flipud(pixel_to_wavelen_map)
    # flip upside down because BATC tend to have lower wavleengths on TOP
    dframe = pd.read_csv(os.path.join(rad_cal_dir, radiance), delimiter=",")
    sphere_wavelen = dframe['Wavelength'].values
    sphere_radiances = np.array([dframe['Radiance_VISIR_NOfilt+Y'].values,
                                 dframe['Radiance_VISIR_NOfilt-Y'].values,
                                 dframe['Radiance_VISIR_NOfilt_Center'].values])

    radiance_all = []
    for i in range(0, sphere_radiances.shape[0]):
        radiance_wavelen_map = []
        for wav_l in range(0, pixel_to_wavelen_map.shape[1]):
            tempo_wvl = pixel_to_wavelen_map[:, wav_l]
            sphere_radiance = sphere_radiances[i, :]
            # Now let's perform spectral interpolation
            fit_params_rad = interp1d(sphere_wavelen, sphere_radiance,
                                      kind='linear')
            fitted_val_rad = fit_params_rad(tempo_wvl)
            radiance_wavelen_map.append(fitted_val_rad)
        radiance_all.append(radiance_wavelen_map)
    radiance_all = np.array(radiance_all)
    csv_file_names = [r'Radiance+Y_all_VIS.csv',
                      r'Radiance-Y_all_VIS.csv',
                      r'Radiance_center_VIS.csv']
    np.savetxt(os.path.join(radiance_save_dir, csv_file_names[0]),
               radiance_all[0, :, :].T, delimiter=',', fmt='%1.3f')
    np.savetxt(os.path.join(radiance_save_dir, csv_file_names[1]),
               radiance_all[1, :, :].T, delimiter=',', fmt='%1.3f')
    np.savetxt(os.path.join(radiance_save_dir, csv_file_names[2]),
               radiance_all[2, :, :].T, delimiter=',', fmt='%1.3f')

    return



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
    plt.title(title)
    divider = make_axes_locatable(fig_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(image, cax=cax)
    plt.grid(False)
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close('all')


def main():
    """
    This is the main function that does all the processing and calling out various routines
    that ultimately derive the radiometric calibration coefficients.
   """
    home_dir = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer'
    rad_cal_dir = os.path.join(home_dir, r'Radiometric_Calibration')
    pixel_to_wavelen_map = r'pixel_to_wavelen_map.csv'
    sphere_radiance = r'sphere_radiance.csv'
#    dss_files_uv = [r'mean_radiance_DSScenter_UV', r'mean_radiance_DSSminus_UV',\
#                    r'mean_radiance_DSSplus_UV', r'stitched_radiance_image_max_UV']
#    dss_files_vis = [r'mean_radiance_DSScenter_VIS', r'mean_radiance_DSSminus_VIS',\
#                    r'mean_radiance_DSSplus_VIS', r'stitched_radiance_image_max_VIS']

    # define a directory to save the sphere radiance and the radiance to pixel wavelen map
    #create_radiance_plots(rad_cal_dir, sphere_radiance)
    radiance_save_dir = os.path.join(rad_cal_dir, r'Sphere_radiance_map')
    if not os.path.exists(radiance_save_dir):
        os.makedirs(radiance_save_dir)

    create_sph_rad_to_wavelen_map(rad_cal_dir, sphere_radiance,
                                  pixel_to_wavelen_map, radiance_save_dir)

    # Notice this function doesnt return anything. Instead it saves
    #the sphere to wavelen map as a csv file for each sphere postion.
    #The code below will read the CSV files as deisred.


if __name__ == "__main__":
    main()
