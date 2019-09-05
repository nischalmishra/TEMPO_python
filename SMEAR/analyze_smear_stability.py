# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:38:50 2017

@author: nmishra
"""
# This script uses the PRNU map derived from integration time sweep data to some
# of the images acquired during intesnity sweep.

import os
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def get_oultier_mask():
    file_path = r'C:\Users\nmishra\Workspace\TEMPO\outlier_mask'
    file_name = 'final_outlier_mask.csv'
    outlier_mask = genfromtxt(file_path+'/'+file_name, delimiter=',')
    rows, cols = outlier_mask.shape
    mask_A = outlier_mask[0: rows/2, 0:cols/2]
    mask_B = outlier_mask[rows/2:rows, 0:cols/2]
    mask_C = outlier_mask[rows/2:rows, cols/2:cols]
    mask_D = outlier_mask[0:rows/2, cols/2:cols]
    outlier_mask = [mask_A, mask_B, mask_C, mask_D]
    return outlier_mask
    

def filter_outlier_median(quads):
 
    if np.array(quads).ndim ==3:
        ndims, nx_quad, ny_quad = quads.shape
    elif np.array(quads).ndim ==2:      
        ndims=1
        nx_quad, ny_quad = quads.shape
    else:
        nx_quad= 1
        ndims=1
        ny_quad = len(quads)
        
    hist_data = np.reshape(quads,(ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    #print(outlier_filtered_data)
    return outlier_filtered_data

def perform_bias_subtraction_ave (active_quad, trailing_overclocks):
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
    
def perform_linearity(active_quad, quad_name):
    
    # perform linearity based on look up table. 
        
    quad = quad_name.split(' ')[1]
    active_quad[active_quad<0] = 0
    
    path = r'C:\Users\nmishra\Workspace\TEMPO\Linearity_testing\Ping_pong_included\plots_integration_sweep\Look_up_table'
    linearity_file= 'Linearity_look_up_table_final_DN.csv'
    linearity_file = pd.read_csv(path+'/'+ linearity_file)
    # for even quad
    lut_even =  linearity_file[['Input Signal (DN)', quad+str(1)]].dropna()   
    lut_even['DN'] = lut_even.pop('Input Signal (DN)').astype(int)
    lut_even['Value'] = lut_even.pop(quad+str(1)).astype(float)
    lut_even = lut_even.set_index('DN') 
    even_detector_active_quad = active_quad[:, ::2]
    nx_half_quad, ny_half_quad = even_detector_active_quad.shape
    even_detector_active_quad = np.reshape(np.array(even_detector_active_quad), (nx_half_quad*ny_half_quad, 1))
    df = pd.DataFrame(data=even_detector_active_quad)
    df.columns = ['DN']
    datin_even = df['DN'].astype(int)
    linearized_quad_even = datin_even.apply(lambda x:lut_even.ix[x])
    even_pixels_quad = np.reshape(linearized_quad_even.values,(nx_half_quad, ny_half_quad))    
    df = None
    
    lut_odd =  linearity_file[['Input Signal (DN)', quad+str(2)]].dropna()   
    lut_odd['DN'] = lut_odd.pop('Input Signal (DN)').astype(int)
    lut_odd['Value'] = lut_odd.pop(quad+str(2)).astype(float)
    lut_odd = lut_odd.set_index('DN')    
    odd_detector_active_quad = active_quad[:, 1::2]
    
    nx_half_quad, ny_half_quad = odd_detector_active_quad.shape
    odd_detector_active_quad = np.reshape(np.array(odd_detector_active_quad), (nx_half_quad*ny_half_quad, 1))
    df = pd.DataFrame(data=odd_detector_active_quad)
    df.columns = ['DN']
    datin_odd = df['DN'].astype(int) 
    
    linearized_quad_odd = datin_odd.apply(lambda x: lut_odd.ix[x])
    odd_pixels_quad = np.reshape(linearized_quad_odd.values,(nx_half_quad, ny_half_quad))
    # now let's put quad in place
    
    nx_quad,ny_quad = active_quad.shape
    linearized_quad = np.array([[0]*ny_quad]*nx_quad)
    linearized_quad = np.reshape(linearized_quad, (nx_quad, ny_quad))
    linearized_quad[:, ::2] = even_pixels_quad 
    linearized_quad[:, 1::2] = odd_pixels_quad 
    return linearized_quad

def perform_smear_subtraction(active_quad, int_time):
    # the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis. 
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO 
    # tft = 8ms 
    tFT = 8*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)   
      
    smear_subtracted_quad = active_quad - smear_factor[None, :]
#    print(np.shape(smear_factor[None, :]))
#    print(np.shape(active_quad))
#    cc
    #column_average = np.mean(active_quad, axis=0)
    return smear_subtracted_quad, smear_factor

def smear_calculation_overclocks(overclocks, trailing_overclocks, active_quad):
    # Ok perform smear subtraction using parallel overclocks
    # seperate out the even and odd
    overclocks= overclocks[2:, :]
    nx_quad, ny_quad = overclocks.shape
    even_detector_smear_quad = overclocks[:, ::2]       
    odd_detector_smear_quad = overclocks[:, 1::2]
    even_detector_bias = trailing_overclocks[ :, ::2]
    odd_detector_bias = trailing_overclocks[:, 1::2]    
    smear_quad = np.array([[0]*ny_quad]*nx_quad)
    even_smear = even_detector_smear_quad- np.mean(filter_outlier_median(even_detector_bias))    
    odd_smear = odd_detector_smear_quad- np.mean(filter_outlier_median(odd_detector_bias))
    smear_quad = np.reshape(smear_quad,(nx_quad, ny_quad))
    smear_quad[:, ::2] = even_smear
    smear_quad[:, 1::2] = odd_smear     
    avg_smear = np.mean(smear_quad, axis=0)
    smear_subtracted_quad = active_quad - avg_smear[None, :] 
    
    return smear_subtracted_quad, avg_smear
       

def plot_smear_signal(smear_signal,column_average, title, figure_name, colors):
    plt.figure()
    plt.plot(column_average[3:], smear_signal[3:], '.',markersize=3,color=colors, linestyle='dotted')
    #print((smear_signal.shape))
    plt.title(title ,fontsize=12,fontweight="bold")
    plt.xlabel('Column Average Signal (DN)')
    plt.ylabel('Smear Signal Estimate (DN)')
    plt.grid(True, linestyle=':')
    plt.ylim(0,700)
    plt.show()
    
    #plt.xlim(800, 16000)    
    #plt.savefig(figure_name,dpi=100,bbox_inches="tight")    
    plt.close('all')


def plot_smear_signal_comp(smear_signal_active_quad, smear_signal_overclock, title, figure_name):
    plt.figure()
    plt.plot(smear_signal_active_quad, '.', markersize=1, color='red', linestyle='dotted', label='using active quad ')
    #plt.plot(smear_signal_overclock, '.', markersize=1, color='blue', linestyle='dotted', label = 'using parallel overclocks')
    plt.title(title ,fontsize=12,fontweight="bold")
    plt.xlabel('Spatial Index (#)')
    plt.ylabel('Smear Signal Estimate (DN)')
    plt.grid(True, linestyle=':')
    plt.ylim(0.9, 1)
    #plt.text( 50, 350, '** Note: Row 200',)
    plt.legend(loc='upper right')  
    plt.show()
    cc      
    #plt.savefig(figure_name,dpi=100,bbox_inches="tight")     #plt.show()
    #plt.show()
    #cc
    
    plt.close('all') 


def plot_smear_signal_correction_comp(active_quad,active_quad_corr, title, figure_name):
    plt.figure()
    row = [50, 600]
    ylim = [1000,16000]
    for i in range(0, 2):
        plt.plot(active_quad[row[i], :], '.', markersize=1, color='red', linestyle='dotted', label='Original Image ')
        plt.plot(active_quad_corr[row[i], :], '.', markersize=1, color='blue', linestyle='dotted', label='after SMEAR removal')
        plt.title(title, fontsize=12, fontweight="bold")
        plt.xlabel('Spatial Index (#)')
        plt.ylabel('Smear Corrected Image Signal (DN)')
        plt.grid(True, linestyle=':')
        plt.ylim(0.9, 1)
        plt.text( 50, ylim[i]-0.16*ylim[i], '** Note: Slice of row ' + str(row[i]), color='blue')
        plt.xlim(0,1200)
        plt.legend(loc='best')  
        #plt.savefig(figure_name+str(row[i])+'_slice.png',dpi=100,bbox_inches="tight") 
        #plt.show()
        #cc        
        plt.close('all') 


    
    
def create_image(image_data, title, figure_name):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)
    plt.show()
    cc    
    
    #plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    #plt.show()
    #cc
    plt.close('all')   
    
    
def create_hist(image, title, figure_name) : 
    
    if np.array(image).ndim ==2:      
        
        nx_quad, ny_quad = image.shape
    else:
        nx_quad= 1        
        ny_quad = len(image)
        #print(ny_quad)
        #cc
    
    label = 'Mean = '+ str(round(np.mean(image), 3)) +'%'            
    plt.figure(figsize=(8, 5))
    plt.hist(np.reshape(image, (nx_quad* ny_quad, 1)), 100, facecolor='red', label=label)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.xlim(-1, 10)
    #plt.ylim(0, 40000)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' % difference [100*(Overclock-Active Quad)/Active Quad]  ', fontsize=12,
              fontweight="bold")
    plt.title(title)
    #plt.show()
    #plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    #plt.show()
    #cc
    plt.close('all')
    
    
    
def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2
    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark'
    if 'Integration_Sweep' in file_path1:
    #if 'Bar_target' in file_path1:
        saturated_collects = [ 'FT6_LONG_INT_130018.dat.sav',#'FT6_SHORT_INT_0.dat.sav',
                              'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']

    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
   
  
   
    nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files] # for cross talk image
    
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark\saved_quads\smear_dark_image'
    PRNU_directory = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads'
    
    for i in range(0, 4): # for the 4 quads
        i=0
        column_average_all = [ ]
        column_average_all_overclock = []
        smear_signal_all = [ ]
        smear_signal_all_overclock =  []
        systematic_diff = []
        active_quad_all = []
        diff_image_mean = [ ]
        diff_image_std = [ ]
        int_time_all = [ ]
        smear_signal_all = [ ]        
        for data_files in nominal_int_files[1:]:
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            IDL_variable = readsav(data_file)
            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                string1 = 'VA_'
                string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '           
            
            int_time_all.append(int_time)
            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data

            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            color = ['red','blue', 'green','magenta']
            all_full_frame = IDL_variable.q
#            print(np.max(all_full_frame))
#            cc
            quad_dir = quads[i]
            PRNU_maps = quads[i] +'/'+'Final_PRNU\Filtered_PRNU'
            PRNU_file = [each for each in os.listdir(os.path.join(PRNU_directory, PRNU_maps)) \
                        if each.endswith('_Final_PRNU.csv')]
            PRNU_mask = os.path.join(PRNU_directory, PRNU_maps, PRNU_file[0])
            PRNU = genfromtxt(PRNU_mask, delimiter=',')
            all_PRNU= np.ones((1024, 1024))
            all_PRNU[9:1000, :] = PRNU

            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)       
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)            
            smear_oclck_A = np.mean(quad[:, 1028:1046, 10:1034], axis=0)
            tsoc_smear = np.mean(quad[:, 1028:1046, 1034:1056], axis=0)
        
            #------perform bias subtraction using trailing overclocks----------
            active_quad_A = perform_bias_subtraction_ave(quad_A[4:1028, 10:1034], tsoc_A)
            
            
            # perform non-lineairty correction using look up table
            linearized_quad_A = perform_linearity(active_quad_A, quads[i])
            
            # Skip lineairty if you are using a bright homegeneous photon transfer data
            #linearized_quad_A = active_quad_A
            
            # Perform dark current subtraction using dark data
            
            
            
            quad_save = 'saved_quads'
            save_dir_image = os.path.join(save_dir, quad_dir,quad_save)
            if not os.path.exists(save_dir_image):
               os.makedirs(save_dir_image)
            title = 'Bar Target  Image, '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
           # print(np.min(linearized_quad_A))
           # cc
            #create_image(linearized_quad_A , title, figure_name)
            #linearized_quad_A = linearized_quad_A - 0.0015*np.fliplr(linearized_quad_A)
           
            # now estimate smear using two methods and then subtract out the smear
            # First use the overclock method. Remember to subtract offset from the smear overclock too  
            smear_subtracted_quad, smear_signal_overclock = smear_calculation_overclocks(smear_oclck_A, tsoc_smear, linearized_quad_A)
            
            # Now the method of active quad
            
            active_quad_A_corr, smear_signal_active_quad = perform_smear_subtraction(linearized_quad_A, int_time)
            smear_signal_all.append(smear_signal_active_quad)
            
            active_quad_A_corr = active_quad_A_corr
            #smear_subtracted_quad = smear_subtracted_quad/all_PRNU
            #---------save active quad smear corrected image-----------------------
            quad_save_active = 'saved_quads/active_method'
            save_dir_image = os.path.join(save_dir, quad_dir,quad_save_active)
            if not os.path.exists(save_dir_image):
               os.makedirs(save_dir_image)
            title = 'Smear Corrected Image (Using Active quad Method)\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image + '/'+ string1 + str(int_time) + '_image.png'                                               
            #create_image(active_quad_A_corr, title, figure_name) 
            
            
            #-------------- saved overclock smear corrected image---------------
            quad_save_overclock = 'saved_quads/overclock_method'
            save_dir_image_oclck = os.path.join(save_dir, quad_dir,quad_save_overclock)
            if not os.path.exists(save_dir_image_oclck):
               os.makedirs(save_dir_image_oclck)
            title = 'Smear Corrected Image (Using Overclock Method)\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image_oclck + '/'+ string1 + str(int_time) + '_image.png'                                               
            #create_image(smear_subtracted_quad, title, figure_name) 
            
            #smear_subtracted_quad= None # clear memory
            
            #------------- plot difference between two calculated smears-----------
            plot_diff_dir_name = 'smear_diff_betn_methods/diff_plot'
            plot_diff_dir = os.path.join(save_dir, quad_dir,plot_diff_dir_name)
            if not os.path.exists(plot_diff_dir):
               os.makedirs(plot_diff_dir)  
            title = 'Estimate of SMEAR Signal using two methods \n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= plot_diff_dir + '/'+ string1 + str(int_time) + '_diff_plot.png'
            plot_smear_signal_comp(smear_signal_active_quad, smear_signal_overclock, title, figure_name)
            cc
           
            
            quad_save_correction_comp = 'saved_quads/active_method/Correction_comp'
            save_dir_image_comp = os.path.join(save_dir, quad_dir,quad_save_correction_comp)
            if not os.path.exists(save_dir_image_comp):
               os.makedirs(save_dir_image_comp)
            title = 'Profile Along Spatial Direction \n'+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= save_dir_image_comp + '/'+ string1 + str(int_time) + '_plot'
             
            #plot_smear_signal_correction_comp(linearized_quad_A, active_quad_A_corr, title, figure_name)
             
            # Calculate differences for systematic uncertainty
            
            
            # Apply PRNU
            #active_quad_A =  np.true_divide(active_quad_A, PRNU)
            #smear_subtracted_quad = np.true_divide(smear_subtracted_quad[9:1000, :], PRNU)
            
            # Now let's take a difference between the two smear subtracted quads
            #print(smear_subtracted_quad)
            #difference_image = smear_subtracted_quad - active_quad_A_corr
           # difference_smear = (smear_signal_overclock - smear_signal_active_quad)
  
            diff_dir_name = 'smear_diff_betn_methods/diff_image'
            diff_dir = os.path.join(save_dir, quad_dir,diff_dir_name)
            if not os.path.exists(diff_dir):
               os.makedirs(diff_dir)
            
            title = 'Difference between two smear corrected methods\n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'  
            figure_name= diff_dir + '/'+ string1 + str(int_time) + '_diff_image.png'
            #create_image(difference_image, title, figure_name)
            
            # let us plot the histogram of percent difference now
            #pct_diff_image = 100 * np.true_divide( difference_image, smear_subtracted_quad)
           # pct_diff_smear = 100 * np.true_divide( difference_smear, smear_signal_active_quad)
            
            hist_dir = 'smear_diff_betn_methods/diff_image_hist'
            diff_hist_dir = os.path.join(save_dir, quad_dir, hist_dir)
            if not os.path.exists(diff_hist_dir):
               os.makedirs(diff_hist_dir)
            figure_name = diff_hist_dir + '/'+ string1 + str(int_time) + '_diff_hist.png'
            title = 'Histogram of  diff. between two smear corrected methods\n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            #create_hist(filter_outlier_median(pct_diff_image), title, figure_name)
            
            hist_dir_plot = 'smear_diff_betn_methods/diff_plot_hist'
            diff_hist_dir_plot = os.path.join(save_dir, quad_dir, hist_dir_plot)
            if not os.path.exists(diff_hist_dir_plot):
               os.makedirs(diff_hist_dir_plot)
            figure_name = diff_hist_dir_plot + '/'+ string1 + str(int_time) + '_diff_hist_plot.png'
            title = 'Histogram of  difference between two SMEAR Estimates\n' + quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            #create_hist((pct_diff_smear), title, figure_name)       
           
            #diff_image_mean.append(np.mean(filter_outlier_median(pct_diff_image)))
            #print(diff_image_mean)
            #diff_image_std.append(np.std(filter_outlier_median(difference_image))/np.mean(filter_outlier_median(difference_image)))
            
           # Ok let's estimate systematic uncertainty
           #________________________________________________________________
            row = [50, 60, 100, 120, 150, 180, 200]
            for rows in row:
                # for Quad A
                profile1 = np.mean(active_quad_A_corr[row, 430:450])
                profile2 = np.mean(active_quad_A_corr[row, 550:570])
                profile3 = np.mean(active_quad_A_corr[row, 480:520])
                residual_diff = (profile3-(profile1+profile2)/2)/(profile1+profile2)/2
                
                # for quad B
#                profile1 = np.mean(active_quad_A_corr[row, 200:220])
#                profile2 = np.mean(active_quad_A_corr[row, 280:290])
#                profile3 = np.mean(active_quad_A_corr[row, 240:265])
#                residual_diff = (profile3 -(profile1+profile2)/2)/(profile1+profile2)/2
#                int_time_all.append(int_time)
                
                 # for quad c
#                profile1 = np.mean(active_quad_A_corr[row, 200:220])
#                profile2 = np.mean(active_quad_A_corr[row, 300:320])
#                profile3 = np.mean(active_quad_A_corr[row, 240:275])
#                residual_diff = (profile3 -(profile1+profile2)/2)/(profile1+profile2)/2
#                int_time_all.append(int_time)
                
#                 profile1 = np.mean(active_quad_A_corr[row, 300:450])
#                 profile2 = np.mean(active_quad_A_corr[row, 520:550])
#                 profile3 = np.mean(active_quad_A_corr[row, 480:530])
#                 residual_diff = (profile3 -(profile1+profile2)/2)/(profile1+profile2)/2
#                 #int_time_all.append(int_time)
#                 print(profile1)
#                 print(profile2)
#                 print(profile3)
#                 print(residual_diff)
#                 cc
                systematic_diff.append(residual_diff)
         
#_____________________________________________________________________________        
        # save smear file
        
        
        
        smear_diff_save_dir = 'smear_diff_betn_methods/diff_data'
        save_diff_data =  os.path.join(save_dir, quad_dir, smear_diff_save_dir)
        #print(save_diff_data)
        
        if not os.path.exists(save_diff_data):
               os.makedirs(save_diff_data)
        
        csv_name_systematic = save_diff_data + '/'+ 'systematic_diff.csv'   
        int_name_all  = save_diff_data + '/'+ 'all_int_time.csv'
        csv_name_all_smear = save_diff_data + '/'+ 'all_smear_data.csv'
        
        np.savetxt(csv_name_systematic, np.array(systematic_diff), delimiter = ",")
        np.savetxt(int_name_all, np.array( int_time_all), delimiter = ",")
        np.savetxt(csv_name_all_smear, np.array(smear_signal_all), delimiter =",")
        cc
       
            
if __name__ == "__main__":
    main()
    