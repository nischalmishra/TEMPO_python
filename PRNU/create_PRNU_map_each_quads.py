# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

The purpose of this code is to generate PRNU map for each TEMPO quads. 
Over the course of time this script will be modified as we learn more about the
TEMP CCD behavior

"""
from random import randint
import os
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd


#*****************************************************************************
def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def perform_bias_subtraction_ave (active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    spec_pix, spat_pix = active_quad.shape
    bias_subtracted_quad = np.array([[0]*spec_pix]*spat_pix)
    even_detector_bias = trailing_overclocks[ :, ::2]
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
    #print(active_quad.shape)
    #print(smear_factor.shape)
    smear_subtracted_quad = active_quad - smear_factor[None, :]
    return smear_subtracted_quad, smear_factor

def calculate_dark_current(image, i, int_time):
    """ Calculate the dark current based off the dark data
    It takes the filename of Light data and searches for the mathing integration
    time   in the dark data directory, find the required dark data, subtracts off 
    the offset and smear and computes the dark current
    """
    dark_data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Dark'
    data_path_name_split = image.split('_')
    #print(data_path_name_split)
    all_int_files = [each for each in os.listdir(dark_data_dir) \
                         if each.endswith('_'+data_path_name_split[-1])]   
    print(all_int_files)
    
    dark_data_file = os.path.join(dark_data_dir, all_int_files[0])
    IDL_variable = readsav(dark_data_file)            
    all_full_frame = IDL_variable.q 
    quad = all_full_frame[:, i, :, :]
    active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)               
    tsoc = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
    bias_subtracted_quad = perform_bias_subtraction_ave(active_quad, tsoc)
    smear_subtracted_quad, smear_signal = perform_smear_subtraction(bias_subtracted_quad[10:1000, :], int_time)
    return smear_subtracted_quad



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


def plot_smear_signal(smear_signal, title, figure_name):
    plt.figure()
    plt.plot(smear_signal, color='magenta', linestyle='dotted')
    #print((smear_signal.shape))
    plt.title(title ,fontsize=12,fontweight="bold")
    plt.xlabel('Spatial Pixel Indices')
    plt.ylabel('Smear Factor')
    plt.grid(True, linestyle=':')
    plt.ylim(700,900)
#    plt.show()
#    cc
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    plt.close('all')


def fill_vignetted_pixels(subset_quad):
    rows, cols = subset_quad.shape
    final_image = [[0 for x in range(1024)] for y in range(1028)]
    final_image = np.array(final_image)
    x1_rows = list(range(100, rows+100))    
    new_points_end = list(range(rows+100, 1028)) 
    new_points_begin = list(range(0, 100))
    print(rows, cols)
    print(np.array(final_image).shape)   
    final_image[0:100, :] = subset_quad[0:100, :] 
    final_image[900:, :] = subset_quad[672:, :]
    plt.plot(x1_rows, subset_quad[:, 100], 'r')
    plt.plot(new_points_begin, final_image[new_points_begin, 100], 'g')
    plt.plot(new_points_end, final_image[new_points_end, 100], 'm')
    plt.grid(True, linestyle=':')
    plt.show()
    cc
    
    
    
                        


def interpolate_vignetted_pixels(subset_quad):
    # Let's interpolate forward and backward to fill the vginetted rows
    
    #active_quad_A[15:990, :] slice taken for PRNU
    rows, cols = subset_quad.shape
    print(rows, cols)
    x1_rows = list(range(15, rows+15))    
    new_points_end = list(range(rows+15-1, 1028))   
    new_points_begin = list(range(0, 16))    
    y2_all_end = []
    y2_all_begin = []
    for i in range(0, cols):
        fit = np.polyfit(x1_rows, subset_quad[:, i], 10)
   
        line = np.poly1d(fit)
        val = np.polyval(fit, x1_rows)
     
        desired_val_end = line(new_points_end)
        desired_val_begin = line(new_points_begin)
#        plt.figure()
#        plt.plot(x1_rows, val,'r.', label='5th order polyfit')
#        plt.plot(x1_rows, subset_quad[:, i],'g', label='actual data')
#        plt.plot(new_points_begin, desired_val_begin,'b', label='extrapolated (beginning)')
#        plt.plot(new_points_end, desired_val_end,'m', label='extrapolated (end)')
#        plt.grid(True, linestyle=':')
#        plt.legend(loc='best', ncol=1)
#        plt.title('Example of extrapolation to fill vignetted pixels')
#        plt.ylabel('Raw Image - Offset - SMEAR (DN)')
#        plt.xlabel('Spectral Indices (#)')
#        plt.show()
#        
#        plt.figure()
#        plt.plot(100*(val-subset_quad[:, i])/subset_quad[:, i],'k.')
#        plt.ylabel('%Diff')
#        plt.xlabel('Spectral Indices (#)')
#        plt.title('Residue of model and the data')
#        plt.ylim([-1, 1])
#        plt.grid(True, linestyle=':')
#        plt.show()
#        cc
        y2_all_end.append(desired_val_end)
        y2_all_begin.append(desired_val_begin)
    y2_all_end = np.array(y2_all_end).T
    y2_all_begin = np.array(y2_all_begin).T  
    #print(y2_all_begin.shape)
    final_mask = np.zeros((1028, 1024))
    final_mask[0:16, :] = y2_all_begin
    final_mask[16:rows+16, :] = subset_quad
    final_mask[rows+14:, :] = y2_all_end  
    #print(final_mask.shape)
#    plt.plot(final_mask[:, 100],'r')
#    plt.plot(final_mask[:, 500], 'g')
#    plt.show()
    #cc
    return final_mask


def create_image(image_data, title, figure_name):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)
#    plt.show()
#    cc
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    plt.close('all')
  
def create_hist(PRNU_map, title, figure_name) : 
    nx_quad, ny_quad = PRNU_map.shape
    label = 'Mean = '+ str(round(np.mean(PRNU_map), 3)) + \
            '\n Median = '+ str(round(np.median(PRNU_map), 3)) 
            #'\n Std. = '+ str(round(np.std(PRNU_map), 3))+ \
            #'\n Max = '+ str(round(np.max(PRNU_map), 3)) + \
            #'\n Min = '+ str(round(np.min(PRNU_map), 3))
    plt.hist(np.reshape(PRNU_map, (nx_quad* ny_quad, 1)), 400, facecolor='red', label=label)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.xlim(0.8, 1.1)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel('PRNU', fontsize=12,
              fontweight="bold")
    plt.title(title)
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    plt.close('all')
 
def create_striping_metric(sample_rows, sample_cols, active_quads_filtered, title, figure_name):
    nx_quad, ny_quad = active_quads_filtered.shape
    striping_metric_all_rows = []
    striping_metric_all_cols = []
    for i in range(0, nx_quad):
        x = list(active_quads_filtered[i,:])       
        xdiff = [(x[n]-0.5*(x[n+1]+x[n-1]))/x[n] for n in range(1, len(x)-1)]
        striping_metric_all_rows.append(xdiff)
        
    
    nx, ny = np.array(striping_metric_all_rows).shape
    
    for j in range(0, ny_quad):
        y = list(active_quads_filtered[:, j])
        ydiff = [(y[n]-0.5*(y[n+1]+y[n-1]))/y[n] for n in range(1, len(y)-1)]
        striping_metric_all_cols.append(ydiff)
    nx1, ny1 = np.array(striping_metric_all_cols).shape 
          
    label = ['Spatial Direction', 'Spectral  Direction' ]
    #sample_rows = randint(0,nx_quad)
    #sample_cols = randint(0, ny_quad)
    nrows=2
    ncols=1
    fig, axrr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    colors = ['red', 'blue']
    text1 = label [0] + ' (Row = ' + str(sample_rows) +')'
                        
    text2 = label[1] + ' (Column = ' + str(sample_cols)+')'        
    
      
    axrr[0].hist(striping_metric_all_rows[sample_rows]*100 ,
                     50, histtype='stepfilled',color=colors[0], alpha=0.7, label=text1)
    axrr[0].set_ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    axrr[0].set_xlabel('Striping Metric in Spatial Direction (%)', fontsize=12,
              fontweight="bold")
    legend = axrr[0].legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    axrr[0].set_xlim(-0.02, 0.02)
    axrr[0].grid(True, linestyle=':')
    axrr[0].set_title(title, fontsize=12, fontweight="bold")
    
    
    
    axrr[1].hist(striping_metric_all_cols[sample_cols]*100 ,
                     50, histtype='stepfilled',color=colors[1], alpha=0.7, label=text2)
    axrr[1].set_ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    axrr[1].set_xlabel('Striping Metric in Spatial Direction (%)', fontsize=12,
              fontweight="bold")
    legend = axrr[1].legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    axrr[1].set_xlim(-0.02, 0.02)
    axrr[1].grid(True, linestyle=':') 
    #plt.show()
    #cc
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")    
    plt.close('all')
    return striping_metric_all_rows, striping_metric_all_cols
    
def create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name):
      label = ['Spatial Direction', 'Spectral Direction'  ]
      nrows=2
      ncols=1
      fig, axrr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
      axrr[0].plot(100*np.abs(np.array(striping_metric_all_rows)), '.', label=label[0])
      axrr[0].set_xlabel('Pixels Indices (#)',fontsize=14, fontweight="bold" )
      axrr[0].set_ylabel('Spectral Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[0].set_title(title, fontsize=14,fontweight="bold")
      axrr[0].set_xlim(0, 1050)
      axrr[0].set_ylim(0, 14)     
      axrr[0].grid(True, linestyle=':')
      
      axrr[1].plot(100*np.abs(np.array(striping_metric_all_cols)), '.', label=label[1])
      axrr[1].set_xlabel('Pixels Indices (#)',fontsize=14, fontweight="bold" )
      axrr[1].set_ylabel('Spatial Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[1].grid(True, linestyle=':')
      axrr[1].set_ylim(0, 14)
      axrr[1].set_xlim(0, 1050)     
      axrr[1].grid(True, linestyle=':')
#      plt.show()
#      cc
      fig.savefig(figure_name, dpi=100, bbox_inches="tight")
  

def plot_few_PRNU_vals(PRNU_map, figure_name, title):
    plt.figure()
    plt.plot(PRNU_map[:, 10], 'r', label = 'Spatial Index 1')
    plt.plot(PRNU_map[:, 55],'k', label = 'Spatial Index 55')
    plt.plot(PRNU_map[:, 100], 'g', label = 'Spatial Index 100')
    plt.plot(PRNU_map[:, 381],color='orange', label = 'Spatial Index 381')
    plt.plot(PRNU_map[:, 500], 'b', label = 'Spatial Index 500')
    plt.plot(PRNU_map[:, 1000], 'm', label = 'Spatial Index 1000')
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.title(title)
    plt.ylabel('PRNU Value (Ratio)')
    plt.ylim([0.90, 1.05])
    plt.xlabel('Pixel Indices (#)')
    #plt.show()    
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    #plt.show() 
    
def plot_few_PRNU_vals_spatial(PRNU_map, figure_name, title):
    plt.figure()
    plt.plot(PRNU_map[:, 10], 'r', label = 'Spectral Index 1')
    plt.plot(PRNU_map[:, 55],'k', label = 'Spectral Index 55')
    plt.plot(PRNU_map[:, 100], 'g', label = 'Spectral Index 100')
    plt.plot(PRNU_map[:, 381],color='orange', label = 'Spectral Index 381')
    plt.plot(PRNU_map[:, 500], 'b', label = 'Spectral Index 500')
    plt.plot(PRNU_map[:, 1000], 'm', label = 'Spectral Index 1000')
    plt.grid(True, linestyle=':')
    plt.legend(loc='best')
    plt.title(title)
    plt.ylabel('PRNU Value (Ratio)')
    plt.ylim([0.90, 1.05])
    plt.xlabel('Pixel Indices (#)')       
    plt.savefig(figure_name, dpi=100, bbox_inches="tight") 
    #plt.show() 
      
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
        not_used_collects = ['FT6_SHORT_INT_25048.dat.sav','FT6_SHORT_INT_0.dat.sav',
                             'FT6_SHORT_INT_10038.dat.sav', 'FT6_SHORT_INT_19981.dat.sav',
                             'FT6_SHORT_INT_50000.dat.sav', 'FT6_SHORT_INT_54971.dat.sav',
                             'FT6_SHORT_INT_34990.dat.sav','FT6_SHORT_INT_30019.dat.sav',
                             'FT6_SHORT_INT_45028.dat.sav','FT6_SHORT_INT_4971.dat.sav',
                             'FT6_SHORT_INT_39961.dat.sav','FT6_SHORT_INT_54971.dat.sav',
                             'FT6_SHORT_INT_15009.dat.sav', 'FT6_SHORT_INT_10038.dat.sav'
                             'FT6_SHORT_INT_19981.dat.sav', 'FT6_SHORT_INT_15009dat.sav',
                             'FT6_LONG_INT_130018.dat.sav','FT6_LONG_INT_134990.dat.sav', 
                             'FT6_LONG_INT_139961.dat.sav','FT6_LONG_INT_145028.dat.sav', 
                             'FT6_LONG_INT_149999.dat.sav','FT6_LONG_INT_154970.dat.sav',
                             'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                             'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                             'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                             'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                             'FT6_LONG_INT_199999.dat.sav', 'FT6_LONG_INT_125047.dat.sav']
    
    all_int_files = [each for each in os.listdir(file_path1) \
                     if each.endswith('.dat.sav')]
                                                        
    nominal_int_files = [items for items in all_int_files if items not in not_used_collects]
    #print(nominal_int_files)

    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Analysis_5th_order'
    
    for i in range(0, 4): # for the 4 quads
        #i=1         
        all_PRNU_map = [ ]
        for data_files in nominal_int_files:      
            
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
            
            # read the dark data for dark current subtraction
            # perform bias removal using serial overclocks for both dark data and the photon transfer data                
            
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q
            quad = all_full_frame[:, i,:, :]
            quad_A = np.mean(quad[:,:,:], axis=0)                
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)
            
            create_image(quad_A[4:1028, 10:1034], title='check', figure_name='check')                             
            active_quad_A = perform_bias_subtraction_ave(quad_A[4:1028, 10:1034], tsoc_A)
             # removed the vignetted zones and perform linearity
            linearized_quad_A = perform_linearity(active_quad_A[10:1000, :], quads[i]) 
            print(linearized_quad_A.shape)
            
            active_quad_A, smear_signal = perform_smear_subtraction(linearized_quad_A, int_time)
            dark_current = calculate_dark_current(data_files, i, int_time)
            active_quad_A = active_quad_A - dark_current                 
            # Ok, let's interpolate all the pixels to fullframe
            active_quad_A = interpolate_vignetted_pixels(active_quad_A)
           #*******************************************************************************************
            quad_dir = quads[i]
            orig_image_dir = 'saved_quads'
            plot_dir = os.path.join(save_dir, quad_dir, orig_image_dir)
            if not os.path.exists(plot_dir):
               os.makedirs(plot_dir) 
           
           #Let's plot the image of the active quad 
            figure_name= plot_dir+'/'+ string1+ str(int_time)+'_new_small.png'
            title = 'Active Pixels in '+quads[i]+', ' + string2 + str(int_time)
            create_image(active_quad_A, title, figure_name)

             # Ok, let's plot the smear signal
            smear_dir_name = 'smear_plot'
            smear_dir = os.path.join(save_dir, quad_dir,smear_dir_name)
            if not os.path.exists(smear_dir):
               os.makedirs(smear_dir)
            figure_name = smear_dir+ '/'+ string1+ str(int_time)+'_new_small.png'
            title = 'Estimate of Mean Smear, '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs' 
            plot_smear_signal(smear_signal, title, figure_name)
            
            #*********************************************************************
            # Toss out the outliers 
            #active_quads_filtered = filter_outlier_median(active_quads)                
            #******************************************************************
            
            #******************************************************************
            # The now divide each pixel by scene mean. This will be one of the
            #ways to generate a PRNU mask
            mean_val = np.mean(filter_outlier_median(active_quad_A))
            
            if mean_val == 0:
                mean_val = np.median(active_quad_A)
            PRNU_map = np.true_divide(active_quad_A, mean_val)
            
            PRNU_plot_name = 'PRNU_plots_spectral_dir'
            PRNU_dir = os.path.join(save_dir, quad_dir, PRNU_plot_name)
            if not os.path.exists(PRNU_dir):
               os.makedirs(PRNU_dir)
            figure_name = PRNU_dir + '/'+ string1+ str(int_time)+'_PRNU_plot' 
            title = 'PRNU values along spectral direction\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'                               
            plot_few_PRNU_vals(PRNU_map, figure_name, title)
            
            PRNU_plot_name = 'PRNU_plots_spatial_dir'
            PRNU_dir = os.path.join(save_dir, quad_dir, PRNU_plot_name)
            if not os.path.exists(PRNU_dir):
               os.makedirs(PRNU_dir)
            figure_name = PRNU_dir + '/'+ string1+ str(int_time)+'_PRNU_plot' 
            title = 'PRNU values along satial direction\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'                               
            plot_few_PRNU_vals_spatial(PRNU_map.T, figure_name, title)
            
                     
            #******************************************************************

            #****************STRIPING METRIC SECTION**************************
            
            # let us plot the striping metric for each integration time
            # begin with histograms
            striping_directory = 'striping_metric_hist'
            striping_met_directory = os.path.join(save_dir, quad_dir, striping_directory)
            if not os.path.exists(striping_met_directory):
                os.makedirs(striping_met_directory)
            figure_name = striping_met_directory +'/'+ string1+ str(int_time)+'_new_small.png'
            title = 'Typical Striping Metric In TEMPO FPS Level Data, '+quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            nx_quad, ny_quad = active_quad_A.shape
            sample_rows = randint(0, nx_quad)
            sample_cols = randint(0, ny_quad)
            striping_metric_all_rows, striping_metric_all_cols = create_striping_metric(sample_rows, sample_cols, active_quad_A, title, figure_name)
            
            # now the actual plot
            striping_plot_directory = 'striping_metric_plots'
            striping_met_plot_directory = os.path.join(save_dir, quad_dir, striping_plot_directory)
            if not os.path.exists(striping_met_plot_directory):
                os.makedirs(striping_met_plot_directory)
            figure_name = striping_met_plot_directory +'/'+ string1+ str(int_time)+'_new_small.png'
            title =  'Striping metric plot, '+quads[i] +', '+ string2+ str(int_time) + r" $\mu$" +'secs'
            create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name)
            
            # now create a striping metric image. THe purspose is to understand
            # the most homogeneous part of the CCD
            striping_image_directory_name = 'striping_image_plot'
            spatial_striping_image_directory = os.path.join(save_dir, quad_dir,
                                                            striping_image_directory_name,
                                                            'Spatial')
            
            spectral_striping_image_directory = os.path.join(save_dir,quad_dir,
                                                            striping_image_directory_name,
                                                            'Spectral')
            if not os.path.exists(spatial_striping_image_directory):
                os.makedirs(spatial_striping_image_directory)
                
            if not os.path.exists(spectral_striping_image_directory):
                os.makedirs(spectral_striping_image_directory)
                
            figure_name = spatial_striping_image_directory +'/'+ string1+ str(int_time)+'_new_small.png'
            title = 'Spatial Striping Metric Image, '+quads[i]+', '+ string2 + str(int_time)+ r" $\mu$" +'secs' 
            create_image(np.abs(striping_metric_all_rows), title, figure_name)
            
            figure_name = spectral_striping_image_directory +'/'+ string1+ str(int_time)+'_new_small.png'
            title = 'Spectral Striping Metric Image, '+quads[i]+ ', '+ string2 + str(int_time) + r" $\mu$" +'secs'
            create_image(np.abs(np.array(striping_metric_all_cols).T), title, figure_name)

            # lets' plot the PRNU map
            # Let's create a separate directory
            PRNU_image_dir = 'PRNU_map'
            PRNU_hist_dir = 'PRNU_hist'
            image_dir = os.path.join(save_dir,quad_dir, PRNU_image_dir)
            if not os.path.exists(image_dir):
                os.makedirs(image_dir)
            hist_dir = os.path.join( save_dir, quad_dir, PRNU_hist_dir)
            if not os.path.exists(hist_dir):
                os.makedirs(hist_dir)
   
            #plot the PRNU map
            figure_name  = image_dir+'/'+ string1+ str(int_time)+'_new_small.png'
            title = 'PRNU Map, '+quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            create_image(PRNU_map, title, figure_name)
            
            csv_file_name = image_dir+'/'+ quads[i]+string1 + str(int_time)+'_new_small.csv'
            np.savetxt(csv_file_name, np.array(PRNU_map), delimiter=",")
          
            # plot the histograms of PRNU map
            figure_name = hist_dir+'/'+ string1+ str(int_time)+'._new_small.png'
            title = 'Histogram of PRNU Map, '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            create_hist(PRNU_map, title, figure_name)
            all_PRNU_map.append(PRNU_map)
           
            # Ok, take average of all the PRNU maps and compute a final mask
        
        final_PRNU = np.mean(np.array(all_PRNU_map), axis=0)
        final_PRNU_std = np.std(np.array(all_PRNU_map), axis=0)
    
        PRNU_ave_dir = 'Final_PRNU/Full_quad_PRNU'  
        quad_dir = quads[i]
        final_image_dir = os.path.join(save_dir, quad_dir, PRNU_ave_dir)
        if not os.path.exists(final_image_dir):
            os.makedirs(final_image_dir) 
        
        title = 'Average PRNU Map, ' + quads[i]
        figure_name = final_image_dir+'/'+ 'final_mask_image'+'.png'
        create_image(final_PRNU, title, figure_name) 
        
        title = ' Uncertainty associated with average PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'unct_final_mask_image'+'.png'
        create_image(100* final_PRNU_std/final_PRNU, title, figure_name)
        
        title = ' Histogram of Average PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'final_mask_hist'+'.png' 
        create_hist(final_PRNU, title, figure_name)
        
        title = ' Histogram of Uncertainty associated with final PRNU Map, '+ quads[i]
        figure_name = final_image_dir+'/'+ 'unct_final_mask_hist'+'.png' 
        #create_hist(100* final_PRNU_std/final_PRNU, title, figure_name)
        
        csv_file_name = final_image_dir+'/'+ quads[i]+'_Final_PRNU.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU), delimiter=",")
        
        csv_file_name = final_image_dir+'/'+ quads[i]+'_Final_PRNU_Std.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU_std), delimiter=",")
        
if __name__ == "__main__":
    main()
