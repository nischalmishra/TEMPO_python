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

import os
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io.idl import readsav
import matplotlib.pyplot as plt

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
    nx_quad,ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[0]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[ :, ::2]
    avg_bias_even = np.mean(even_detector_bias, axis=1)  
    odd_detector_bias = trailing_overclocks[:, 1::2]
    avg_bias_odd = np.mean(odd_detector_bias, axis=1)
    even_detector_active_quad = active_quad[:, ::2]     
    odd_detector_active_quad = active_quad[:, 1::2]    
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, None] 
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, None] 
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ny_quad, nx_quad))    
    bias_subtracted_quad[:, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad


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
 
def create_striping_metric(active_quads_filtered, title, figure_name):
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
          
    label = ['Spectral Direction', 'Spatial Direction' ]
    colors = ['red', 'blue']
    text1 = 'Max = ' + str(round(100*np.max(np.abs(striping_metric_all_rows)),2))+'%, ' + \
             'Mean = ' + str(round(100*np.mean(np.abs(striping_metric_all_rows)),2))+'%'             
    text2 = 'Max = ' + str(round(100*np.max(np.abs(striping_metric_all_cols)),2))+'%, ' + \
              'Mean = ' + str(round(100*np.mean(np.abs(striping_metric_all_cols)),2))+'%'            
    
    plt.figure(figsize=(8,6))
    plt.hist(np.reshape(100*np.abs(striping_metric_all_rows),(nx*ny,1)) ,
                     100, histtype='stepfilled',color=colors[0], alpha=0.7, label=label[0])
    
    plt.hist(np.reshape(100*np.abs(striping_metric_all_cols),(nx1*ny1,1)), 
                       100, color=colors[1], histtype='stepfilled', alpha=0.7, label=label[1])
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel('Striping Metric (%)', fontsize=12,
              fontweight="bold")
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    plt.xlim(-0.05, 2)
    plt.grid(True, linestyle=':')
    plt.title(title)
    plt.text(1.25, 750000, text1, color=colors[0])
    plt.text(1.25,600000, text2, color=colors[1] )
    
    plt.savefig(figure_name,dpi=100,bbox_inches="tight")    
    plt.close('all')
    return striping_metric_all_rows, striping_metric_all_cols
    
def create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name):
      label = ['Spectral Direction', 'Spatial Direction' ]
      nrows=2
      ncols=1
      fig, axrr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
      axrr[0].plot(100*np.abs(np.array(striping_metric_all_rows)), '.', label=label[0])
      axrr[0].set_xlabel('Spectral Pixels',fontsize=14, fontweight="bold" )
      axrr[0].set_ylabel('Spectral Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[0].set_title(title, fontsize=14,fontweight="bold")
      axrr[0].set_xlim(0, 2100)
      axrr[0].set_ylim(0, 14)     
      axrr[0].grid(True, linestyle=':')
      
      axrr[1].plot(100*np.abs(np.array(striping_metric_all_cols)), '.', label=label[1])
      axrr[1].set_xlabel('Spatial Pixels',fontsize=14, fontweight="bold" )
      axrr[1].set_ylabel('Spatial Striping Metric (%)',fontsize=14, fontweight="bold" )
      axrr[1].grid(True, linestyle=':')
      axrr[1].set_ylim(0, 14)
      axrr[1].set_xlim(0, 2100)     
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
    file_path2 = r'F:\TEMPO\Data\GroundTest\FPS\Intensity_Sweep\Light\Saved_quads'
    filepaths = [file_path1, file_path2]
    for file_path in filepaths:
        
        all_PRNU_map = [] 
        if 'Integration_Sweep' in file_path:
            saturated_collects = ['FT6_LONG_INT_130018.dat.sav','FT6_LONG_INT_134990.dat.sav', 
                        'FT6_LONG_INT_139961.dat.sav','FT6_LONG_INT_145028.dat.sav', 
                        'FT6_LONG_INT_149999.dat.sav','FT6_LONG_INT_154970.dat.sav',
                        'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                        'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                        'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                        'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                        'FT6_LONG_INT_199999.dat.sav']
        
        elif 'Intensity_Sweep' in file_path:
            saturated_collects = ['162_OP_INT_118000.dat.sav','164_OP_INT_118000.dat.sav',
                                  '166_OP_INT_118000.dat.sav','168_OP_INT_118000.dat.sav',
                                  '170_OP_INT_118000.dat.sav','172_OP_INT_118000.dat.sav',
                                  '174_OP_INT_118000.dat.sav','176_OP_INT_118000.dat.sav',
                                  '178_OP_INT_118000.dat.sav','180_OP_INT_118000.dat.sav',
                                  '182_OP_INT_118000.dat.sav','184_OP_INT_118000.dat.sav',
                                  '186_OP_INT_118000.dat.sav','188_OP_INT_118000.dat.sav',
                                  '190_OP_INT_118000.dat.sav','192_OP_INT_118000.dat.sav',
                                   '194_OP_INT_118000.dat.sav','196_OP_INT_118000.dat.sav',
                                   '198_OP_INT_118000.dat.sav','200_OP_INT_118000.dat.sav',
                                   '202_OP_INT_118000.dat.sav','254_OP_INT_118000.dat.sav', 
                                   '252_OP_INT_118000.dat.sav', '250_OP_INT_118000.dat.sav',
                                   '248_OP_INT_118000.dat.sav', '246_OP_INT_118000.dat.sav',
                                   '244_OP_INT_118000.dat.sav', '242_OP_INT_118000.dat.sav',]
        
        
        all_int_files = [each for each in os.listdir(file_path) \
                         if each.endswith('.dat.sav')]
                                                           
        nominal_int_files = [items for items in all_int_files if items not in saturated_collects]
    #    print(nominal_int_files)
    #    cc
        save_dir = 'PRNU'
        orig_image_dir = 'saved_quads'
        plot_dir = os.path.join(file_path, save_dir, orig_image_dir)
        if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        for data_files in nominal_int_files:           
                
                data_path_name_split = data_files.split('_')  
                print(data_files)  
                
               
                if 'Intensity_Sweep' in file_path:
                     int_time = data_path_name_split[0]
                     string1 = 'VA_'
                     string2 = 'VA Setting = '
                else:
                    int_time = round(int(data_path_name_split[-1].split('.')[0]))
                    string1 = 'Integ_time_'
                    string2 = 'Int.time = '
                    
                data_file = os.path.join(file_path, data_files)  
                IDL_variable = readsav(data_file)            
                all_full_frame = IDL_variable.quads 
                # identify the regions to create PRNU map and create a sample image                            
                quad_A = all_full_frame[:, 0, :, :]
                tsoc_A = np.mean(all_full_frame[:,0, 4:1028, 1034:1056], axis=0)
                active_quad_A = perform_bias_subtraction_ave(quad_A, tsoc_A)
                
                quad_B = np.mean(all_full_frame[:, 1,4:1028, 10:1034], axis=0)
                tsoc_B = np.mean(all_full_frame[:,1, 4:1028, 1034:1056], axis=0)
                active_quad_B = perform_bias_subtraction_ave (quad_B, tsoc_B)
                
                quad_C = np.mean(all_full_frame[:, 2, 4:1028, 10:1034], axis=0)
                tsoc_C = np.mean(all_full_frame[:,2, 4:1028, 1034:1056], axis=0)
                active_quad_C = perform_bias_subtraction_ave (quad_C, tsoc_C)
                
                quad_D = np.mean(all_full_frame[:, 3, 4:1028, 10:1034], axis=0)
                tsoc_D = np.mean(all_full_frame[:,3, 4:1028, 1034:1056], axis=0)
                active_quad_D = perform_bias_subtraction_ave (quad_D, tsoc_D)
                
                #lower_quads = np.concatenate((active_quad_A[20:900, 10:1034], active_quad_B[20:900, 10:1034]), axis=1)
                #upper_quads = np.concatenate((active_quad_D[20:900, 10:1034], active_quad_C[20:900, 10:1034]), axis=1)
                lower_quads_full_frame =   np.concatenate((quad_A, np.fliplr(quad_B)), axis=1) 
                upper_quads_full_frame = np.concatenate((np.flipud(quad_D), np.rot90(quad_C, 2)), axis=1) 
                active_quads_full_frame = np.concatenate((lower_quads_full_frame, upper_quads_full_frame), axis=0)
                
                # let's crop the image. There are masking issues at the boundary between two CCDs.
                # So I cropped off 24 pixels in spectral direction for each quads.
                lower_quads = np.concatenate((active_quad_A[0:1000, :], np.fliplr(active_quad_B[0:1000, :])), axis=1) 
                upper_quads = np.concatenate((np.flipud(active_quad_D[0:1000, :]), np.rot90(active_quad_C[0:1000, :], 2)), axis=1) 
                active_quads = np.concatenate((lower_quads, upper_quads), axis=0)
      
               #*******************************************************************************************
               # Let's plot the image of the active quad 
                figure_name  = plot_dir+'/'+ string1 + str(int_time)+'.png'
                title = 'Active Quad Image (Before PRNU), ' + string2 + str(int_time)
                # plot the image
                create_image(active_quads_full_frame, title, figure_name)            
               
                #**********************************************************************
                # Toss out the outliers 
                active_quads_filtered = filter_outlier_median(active_quads)                
                #******************************************************************
                
                #******************************************************************
                # The now divide each pixel by scene mean. This will be one of the
                #ways to generate a PRNU mask
                PRNU_map = active_quads/np.median(active_quads)
                #******************************************************************
    
    
                #****************STRIPING METRIC SECTION**************************
                
                # let us plot the striping metric for each integration time
                # begin with histograms
                striping_directory = 'striping_metric_hist'
                striping_met_directory = os.path.join(file_path, save_dir, striping_directory)
                if not os.path.exists(striping_met_directory):
                    os.makedirs(striping_met_directory)
                figure_name = striping_met_directory +'/'+ string1+ str(int_time)+'.png'
                title = 'Striping metric, ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
                striping_metric_all_rows, striping_metric_all_cols = create_striping_metric(active_quads_filtered, title, figure_name)
                
                # now the actual plot
                striping_plot_directory = 'striping_metric_plots'
                striping_met_plot_directory = os.path.join(file_path, save_dir, striping_plot_directory)
                if not os.path.exists(striping_met_plot_directory):
                    os.makedirs(striping_met_plot_directory)
                figure_name = striping_met_plot_directory +'/'+ string1+ str(int_time)+'.png'
                title =  'Striping metric plot, ' + string2+ str(int_time) + ', image size = '+ str(active_quads.shape)
                create_striping_metric_plot(striping_metric_all_rows, striping_metric_all_cols, title, figure_name)
                
                # now create a striping metric image. THe purspose is to understand
                # the most homogeneous part of the CCD
                striping_image_directory_name = 'striping_image_plot'
                striping_image_directory = os.path.join (file_path, save_dir,striping_image_directory_name )
                if not os.path.exists(striping_image_directory):
                    os.makedirs(striping_image_directory)
                figure_name = striping_image_directory +'/'+ string1+ str(int_time)+'.png'
                title = 'Striping Metric Image, ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
                #create_image(np.abs(striping_metric_all_cols), title, figure_name)
                
                
                # lets' plot the PRNU map
                # Let's create a separate directory
                PRNU_image_dir = 'PRNU_map'
                PRNU_hist_dir = 'PRNU_hist'
                image_dir = os.path.join(file_path, save_dir, PRNU_image_dir)
                if not os.path.exists(image_dir):
                    os.makedirs(image_dir)
                hist_dir = os.path.join(file_path, save_dir, PRNU_hist_dir)
                if not os.path.exists(hist_dir):
                    os.makedirs(hist_dir)
                
                
                #plot the PRNU map
                figure_name  = image_dir+'/'+ string1+ str(int_time)+'.png'
                title = 'PRNU Map, ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
                create_image(PRNU_map, title, figure_name)
                
                csv_file_name = image_dir+'/'+ string1 + str(int_time)+'.csv'
                np.savetxt(csv_file_name, np.array(PRNU_map), delimiter=",")
               
                
                # plot the histograms of PRNU map
                figure_name = hist_dir+'/'+ string1+ str(int_time)+'.png'
                title = 'Histogram of PRNU Map, ' + string2 + str(int_time) + ', image size = '+ str(active_quads.shape)
                create_hist(PRNU_map, title, figure_name)
                           
    
                all_PRNU_map.append(PRNU_map)
                PRNU_map = None
                all_full_frame = None
                
        # Ok, take average of all the PRNU maps and compute a final mask
        final_PRNU = np.mean(all_PRNU_map, axis=0)
        
        PRNU_ave_dir = 'Final_PRNU'   
        final_image_dir = os.path.join(file_path, save_dir, PRNU_ave_dir)
        if not os.path.exists(final_image_dir):
            os.makedirs(final_image_dir)  
        title = 'Average PRNU Map (Intensity sweep)'
        figure_name = final_image_dir+'/'+ 'final_mask_image'+'.png'
        create_image(final_PRNU, title, figure_name) 
        
        title = 'Histogram of Average PRNU Map (Intensity sweep)'
        figure_name = final_image_dir+'/'+ 'final_mask_hist'+'.png' 
        create_hist(final_PRNU, title, figure_name)
        
        csv_file_name = final_image_dir+'/'+ 'Final_PRNU.csv'
        np.savetxt(csv_file_name, np.array(final_PRNU), delimiter=",")
        
        #title = 'Histogram of Average Striping Metric (Integration Time Sweep)'
        #figure_name = final_image_dir + '/'+ 'final_striping_metric'+'.png'
        #create_striping_metric (final_PRNU, title, figure_name)

    
if __name__ == "__main__":
    main()
