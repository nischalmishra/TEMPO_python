# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:43:46 2017

@author: nmishra

The purpose of this code is to identify outliers on TEMPO CCDs. To begin with,
it takes the raw focal plane plane data in binay format.
This bascailly takes Dave F. IDL variables and does the analysis in python
The hardcoded variables are standard to tempo instruments.
The tests are based on both the photon transfer data as well as well as dark
current data

"""
from random import randint
import os
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
import scipy


#*****************************************************************************
def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def perform_bias_subtraction(active_quad, trailing_overclocks):
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
    avg_bias_even = np.mean(even_detector_bias)  
    odd_detector_bias = trailing_overclocks[:, 1::2]
    odd_samples = odd_detector_bias[:, 4:]
    rows, cols = odd_samples.shape   
    odd_detector_bias = odd_samples[:, 0:cols-1] 
    avg_bias_odd = np.mean(odd_detector_bias)
    even_detector_active_quad = active_quad[:, ::2]     
    odd_detector_active_quad = active_quad[:, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (nx_quad, ny_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad



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


def create_image(image_data, title, figure_name):
    plt.figure()
    ax = plt.gca()
    image = ax.imshow(image_data, cmap='nipy_spectral', origin='lower')
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)
    #plt.show()
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")
    plt.close('all')
  
def create_hist(PRNU_map, title, figure_name) : 
    nx_quad, ny_quad = PRNU_map.shape
    label = 'Mean = '+ str(round(np.mean(PRNU_map), 3)) + \
            '\n Median = '+ str(round(np.median(PRNU_map), 3)) + \
            '\n Std. = '+ str(round(np.std(PRNU_map), 3))+ \
            '\n Max = '+ str(round(np.max(PRNU_map), 3)) + \
            '\n Min = '+ str(round(np.min(PRNU_map), 3))
    plt.hist(np.reshape(PRNU_map, (nx_quad* ny_quad, 1)), 100, facecolor='red', label=label)
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
      axrr[0].set_xlabel('Spatial Pixels',fontsize=14, fontweight="bold" )
      axrr[0].set_ylabel('Spatial Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[0].set_title(title, fontsize=14,fontweight="bold")
      axrr[0].set_xlim(0, 1050)
      axrr[0].set_ylim(0, 14)     
      axrr[0].grid(True, linestyle=':')      
      axrr[1].plot(100*np.abs(np.array(striping_metric_all_cols)), '.', label=label[1])
      axrr[1].set_xlabel('Spectral Pixels',fontsize=14, fontweight="bold" )
      axrr[1].set_ylabel('Spectral Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[1].grid(True, linestyle=':')
      axrr[1].set_ylim(0, 14)
      axrr[1].set_xlim(0, 1050)     
      axrr[1].grid(True, linestyle=':')
      #plt.show()
      #cc
      fig.savefig(figure_name, dpi=100, bbox_inches="tight")
        
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
                             'FT6_LONG_INT_199999.dat.sav']
  
    all_int_files = [each for each in os.listdir(file_path1) \
                         if each.endswith('.dat.sav')]                                                   
    nominal_int_files = [items for items in all_int_files if items not in not_used_collects]
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\PRNU_Quads'
    
    for i in range(0, 4): # for the 4 quads
        #i=1         
        all_PRNU_map = [ ]
        for data_files in nominal_int_files:      
            print(data_files)            
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
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q       
            quad = all_full_frame[:, i, :, :]
            quad_A = np.mean(quad[:, :, :], axis=0)                
            tsoc_A = np.mean(quad[:, 4:1028, 1034:1056], axis=0)                             
            active_quad_A = perform_bias_subtraction_ave(quad_A[4:1028, 10:1034], tsoc_A)
            active_quad_A, smear_signal = perform_smear_subtraction(active_quad_A[9:1000, :], int_time)

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
            smear_dir = os.path.join(save_dir, quad_dir, smear_dir_name)
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
            mean_val = np.mean(active_quad_A)
            if mean_val == 0:
                mean_val = np.median(active_quad_A)
            PRNU_map = np.true_divide(active_quad_A, mean_val)
            #active_quads = np.true_divide(active_quad_A, PRNU_map)
                     
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
            hist_dir = os.path.join(save_dir, quad_dir, PRNU_hist_dir)
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
        PRNU_ave_dir = 'Final_PRNU/Filtered_PRNU'  
        quad_dir = quads[i]
        final_image_dir = os.path.join(save_dir, quad_dir, PRNU_ave_dir)
        if not os.path.exists(final_image_dir):
            os.makedirs(final_image_dir) 
        
        title = quads[i]+' Average PRNU Map (Integration sweep)'
        figure_name = final_image_dir+'/'+ 'final_mask_image'+'._new_small.png'
        create_image(final_PRNU, title, figure_name) 
        
        title = quads[i]+' Histogram of Average PRNU Map (Integration sweep)'
        figure_name = final_image_dir+'/'+ 'final_mask_hist'+'_new_small.png' 
        create_hist(final_PRNU, title, figure_name)
        
        csv_file_name = final_image_dir+'/'+ quads[i]+'_Check_Small_Final_PRNU.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU), delimiter=",")
        #cc
if __name__ == "__main__":
    main()
