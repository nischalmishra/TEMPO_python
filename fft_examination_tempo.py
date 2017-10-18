# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 12:23:32 2017

@author: nmishra
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
def perform_bias_subtraction_ave(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors for both the active quads and trailing overclocks
    # The trailing overclocks are averaged and the average offset is subtracted
    # from the active quad. This is done for both, ping and pong
    """ Remove offset from active quads. Take care of ping-pong by breaking
        Quads and overclocks into even and odd
        """
    # sepearate out even and odd detectors
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
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
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (nx_quad, ny_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad

def calculate_fft(active_image):
    """Compares FFT before and after PRNU correction
    """
    
    single_line = active_image[600, :]
    print(len(single_line))
#    plt.plot(single_line)
#    plt.show()
#    cc
    fft_val = np.fft.fft(single_line)    
    fft_freq = np.fft.fftfreq(len(single_line))
    i = fft_freq>0
    dc_val = fft_val[i]
    dc_comp = dc_val[0]
    print(np.abs(fft_val[i])[0])  
    plt.plot(1/fft_freq[i], np.abs(fft_val[i]/np.abs(dc_comp)))
    #plt.legend(['Original Image'])
    plt.title('FFT along Spectral Direction, Int time =125047 microsecs, Quad A')
    plt.ylabel('Amplitude (DC Normalized)',fontsize=14,fontweight="bold")
    plt.xlabel('Spatial Frequency (Pixel Units)')
    #plt.xlim(0, 25)
    plt.ylim(0, 0.5) 
    #plt.legend(loc= 'best')
    plt.grid(True, linestyle=':')
  
    plt.show()
    

def calculate_fft_full_frame(full_frame):
    """
    This function calculates the FFT of the given full frame image.
    It creates magnitude & phase plot of the full frame and saves it in a user
    defined directory. While this function doesnt have any analytic purpose, it
    provides a sanity check for the orientation of the quads and frequency
    distribution, if I may, of the quads.

    Input : Full_frame image, quads lengths in each direction
    Output : Image plot, magnitude and phase plot of the full frame image.
    TO DO : Remove the hard coded values such as plot directory. In future the
    paths can come from the main function.
    """
    fft_image = np.fft.fft2(full_frame)
    fshift = np.fft.fftshift(fft_image)
    plt.figure()
    plt.subplot(131), plt.imshow(np.real(full_frame), cmap='seismic')
    plt.title('Input Image'), plt.xticks([]), plt.yticks([])
    plt.subplot(132), plt.imshow(20*np.log10(np.abs(fshift)), cmap='seismic')
    plt.title('Magnitude Spectrum of FFT'), plt.xticks([]), plt.yticks([])
    plt.subplot(133), plt.imshow(np.angle(fshift), cmap='seismic')
    plt.title('Phase Spectrum of FFT'), plt.xticks([]), plt.yticks([])
    plt.show()
    fig_handle = plt.gcf()
    fig_handle.set_size_inches(10.0, 6.0)
    plt.savefig(r"C:\Users\nmishra\Desktop\full_frame_fft.png",
                bbox_inches="tight")  
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    if 'Integration_Sweep' in file_path1:
        saturated_collects = ['FT6_SHORT_INT_0.dat.sav','FT6_LONG_INT_130018.dat.sav', 'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']
        
    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
    nominal_int_files = [items for items in all_int_files if items not in saturated_collects]
    quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
    for i in range(0, 4): # for the 4 quads
        i=3   
        for data_files in nominal_int_files:
            print(data_files)
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)

            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                #string1 = 'VA_'
                #string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                #string1 = 'Integ_time_'
                #string2 = 'Int.time = '+ str(int_time)
            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data
            title = 'FFT along Spatial Direction, ' + quads[i]
            label = '*** Each line is one integration time'
            all_full_frame = IDL_variable.q
            #print(all_full_frame.shape)            
            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)
            active_quad_A = quad_A[4:1028, 10:1034]
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)  
            perform_bias_subtraction_ave(active_quad_A, tsoc_A)
            calculate_fft(active_quad_A)
            cc
             # spatial FFT
            single_line = active_quad_A[600, :]
            fft_val = np.fft.fft(single_line)    
            fft_freq = np.fft.fftfreq(len(single_line))
            j = fft_freq>0
            dc_val = fft_val[j]
            dc_comp = dc_val[0]
            print(np.abs(fft_val[j])[0])  
            plt.plot(1/fft_freq[j], np.abs(fft_val[j]/np.abs(dc_comp)), 
                     label=label)
            #plt.legend(['Original Image'])
            plt.title(title)
            plt.ylabel('Amplitude (DC Normalized)', fontsize=14, 
                       fontweight="bold")
            plt.xlabel('Spatial Frequency (Pixel Units)')
            plt.xlim(0, 20)
            plt.ylim(0, 0.2)           
            plt.grid(True, linestyle=':')
  
        #plt.legend(loc= 'best')
        plt.show()
    cc
if __name__ == "__main__":
    main()
     