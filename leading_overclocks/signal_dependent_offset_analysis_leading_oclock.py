# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:22:55 2017

@author: nmishra
"""

import os
import numpy as np

#from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
#from outlier_detection import identify_saturation,\
#                              reject_outlier_median,\
#                              reject_outlier_mean,\
#                              create_outlier_mask,\
#                              create_final_mask,\
#                              create_ORed_mask

def get_size(filename):
    """
    This function reads the filename and passes to the main function
    TO DO : N/A
    """
    fileinfo = os.stat(filename)
    return fileinfo

def filter_outlier_median(quads):
    if np.array(quads).ndim == 3:
        ndims, nx_quad, ny_quad = quads.shape
    else:
        ndims = 1
        nx_quad, ny_quad = quads.shape
    hist_data = np.reshape(quads, (ndims*nx_quad*ny_quad, 1))
    diff = abs(hist_data - np.median(hist_data)) # find the distance to the median
    median_diff = np.median(diff) # find the median of this distance
    measured_threshold = diff/median_diff if median_diff else 0.
    outlier_filtered_data = hist_data[measured_threshold < 5.]
    #print(outlier_filtered_data)
    return outlier_filtered_data

def perform_smear_subtraction(active_quad, int_time):
    # the underlying assumption in smear subtraction is that the dark current
    #in the storage region is really small and hence neglected from the analysis.
    #typically, Csmear = tFT / (ti+ tFT) * (AVG[C(w)] - DCStor * tRO
    # tft = 8ms
    tFT = 8*10**(3)
    ti = int_time
    smear_factor = (tFT / (ti+ tFT))* np.mean(active_quad, axis=0)
    #print(smear_factor.shape)
    #cc
    smear_subtracted_quad = active_quad - smear_factor[None, :]
#    print(np.shape(smear_factor[None, :]))
#    print(np.shape(active_quad))
#    cc
    #column_average = np.mean(active_quad, axis=0)
    return smear_subtracted_quad, smear_factor

def plot_hist(even_samples, odd_samples, title, figure_name):

    if np.array(even_samples).ndim == 3:
        num_dims, nx_quad, ny_quad = even_samples.shape

    elif np.array(even_samples).ndim == 2:

        nx_quad, ny_quad = even_samples.shape
        num_dims = 1
    else:
        nx_quad = 1
        ny_quad = len(even_samples)
        num_dims = 1


    if np.array(odd_samples).ndim == 3:
        num_dims1, nx_quad1, ny_quad1 = odd_samples.shape

    elif np.array(odd_samples).ndim == 2:

        nx_quad1, ny_quad1 = odd_samples.shape
        num_dims1 = 1
    else:
        nx_quad1 = 1
        ny_quad1 = len(odd_samples)
        num_dims1 = 1

    #mean_diff = np.mean(odd_samples)- np.mean(even_samples)
    #text1 = 'Mean Diff (Odd-Even) = ' +  str(round(mean_diff, 2)) +'DN'
    #text2 = 'Uncertainty(Even) = '+ round(100*np.std(even_samples)/np.mean(even_samples))+'%'
    #text3 = 'Uncertainty(Odd) = '+ round(100*np.std(even_samples)/np.mean(even_samples))+'%'
    #print(text2)
    #print(text3)
    plt.figure(figsize=(8, 5))
    plt.hist(np.reshape(even_samples, (num_dims*nx_quad* ny_quad, 1)),
             facecolor='red', label='Even Lines Samples')
    plt.hist(np.reshape(odd_samples, (num_dims1*nx_quad1* ny_quad1, 1)),
             facecolor='blue', label='Odd Lines Samples')
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=1, shadow=True,
                        prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    #plt.xlim(750, 850)
    #plt.ylim(0, 1000)
   # plt.text(815, 300, text1)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
               fontweight="bold")
    plt.xlabel('Signal Counts (DN)', fontsize=12,
               fontweight="bold")
    plt.title(title)

    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.close('all')


def perform_bias_subtraction(active_quad, trailing_overclocks):
    # sepearate out even and odd detectors
    ndims, nx_quad, ny_quad = active_quad.shape
    bias_subtracted_quad = np.array([[[0]*ndims]*ny_quad]*nx_quad)
    even_detector_bias = trailing_overclocks[:, :, ::2]
    avg_bias_even = np.mean(even_detector_bias, axis=2)
    odd_detector_bias = trailing_overclocks[:, :, 1::2]
    avg_bias_odd = np.mean(odd_detector_bias, axis=2)
    even_detector_active_quad = active_quad[:, :, ::2]
    odd_detector_active_quad = active_quad[:, :, 1::2]
    bias_subtracted_quad_even = even_detector_active_quad - avg_bias_even[:, :, None]
    bias_subtracted_quad_odd = odd_detector_active_quad - avg_bias_odd[:, :, None]
    bias_subtracted_quad = np.reshape(bias_subtracted_quad, (ndims, ny_quad, nx_quad))
    bias_subtracted_quad[:, :, ::2] = bias_subtracted_quad_even
    bias_subtracted_quad[:, :, 1::2] = bias_subtracted_quad_odd
    return bias_subtracted_quad

def filter_outlier_trailing(even_samples_avg, odd_samples_avg):
    """ Filter outlier for trailing overclocks"""
    even_samples_avg = even_samples_avg[:, 5:]
    odd_samples = odd_samples_avg[:, 5:]
    rows, cols = odd_samples.shape
    odd_samples_avg = odd_samples[:, 0:cols-1]
    return np.mean(even_samples_avg, axis=1), np.mean(odd_samples_avg, axis=1)




def plot_few_tsocs(even_samples_avg, odd_samples_avg, figure_name, title):
    # let's take the mean tsoc for 100 frames

    nrows = 2
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    #even_samples_avg = even_samples_avg[:, 3:]
    ax[0].plot(even_samples_avg[:, 0], 'b.', label='First')
    ax[0].plot(even_samples_avg[:, 1], 'g.', label='Second')
    ax[0].plot(even_samples_avg[:, 2], 'r.', label='Third')
    ax[0].plot(even_samples_avg[:, 3], '.', color='#9467bd',label='Fourth')
    ax[0].plot(even_samples_avg[:, 4], '.', color='#bcbd22',label='Fifth')
    legend = ax[0].legend(loc='best', ncol=1, shadow=True,
                          prop={'size':16}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    
    ax[0].set_title(title+' (Even Lines)', fontsize=12, fontweight='bold')
    #ax[0].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    ax[0].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')
    ax[0].grid(True, linestyle=':')
    #ax[0].set_ylim(800, 900)
    # Now for the odd
    #print(np.max(even_samples_avg))
#
    odd_samples = odd_samples_avg
    rows, cols = odd_samples.shape
    #odd_samples_avg = odd_samples[:, 2:]
    #print(np.max(odd_samples_avg))
    ax[1].plot(odd_samples_avg[:, 0], 'b.', label='First')
    ax[1].plot(odd_samples_avg[:, 1], 'g.', label='Second')
    ax[1].plot(odd_samples_avg[:, 2], 'r.', label='Third')
    ax[1].plot(odd_samples_avg[:, 3], '.', color='#9467bd',label='Fourth')
    ax[1].plot(odd_samples_avg[:, 4], '.', color='#bcbd22',label='Fifth')
    #print(np.std(odd_samples_avg[:,6])/np.mean(odd_samples_avg[:,6]))
    ax[1].set_title(title+' (Odd Lines)', fontsize=12, fontweight='bold')
    ax[1].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    ax[1].set_ylabel('Serial Overclock Signal (DN)', fontsize=12, fontweight='bold')
#    legend = ax[1].legend(loc='best', ncol=1, shadow=True,
#                          prop={'size':16}, numpoints=1)
#    legend.get_frame().set_edgecolor('wheat')
    plt.savefig(figure_name, dpi=100, bbox_inches="tight")
    plt.grid(True, linestyle=':')
    plt.show()
    cc
    plt.close('all')
    return np.mean(even_samples_avg, axis=1), np.mean(odd_samples_avg, axis=1)


def plot_comparison_overclocks(even_samples_mean_tsoc, odd_samples_mean_tsoc,
                               even_samples_mean, odd_samples_mean,
                               title, figure_name):
    """ plot the line mean of leading and trailing overclocks
    """
    nrows = 2
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)

    diff = 100* (np.mean(even_samples_mean_tsoc) -
                 np.mean(even_samples_mean)) /np.mean(even_samples_mean_tsoc)
    #print(diff)
    text1 = '% Mean Diff = ' + str(round(diff, 2)) +'%'

    diff1 = 100* (np.mean(odd_samples_mean_tsoc) - 
                  np.mean(odd_samples_mean))/np.mean(odd_samples_mean_tsoc)
    #print(diff)
    text2 = '% Mean Diff = ' + str(round(diff1, 2)) + '%'

    ax[0].plot(even_samples_mean_tsoc, '.', color='red', label='Trailing Overclocks')
    ax[0].plot(even_samples_mean, '.', color='blue', label='Leading Overclocks')
    ax[0].set_title(title+' (Even Lines)', fontsize=18, fontweight='bold')
    ax[0].grid(True, linestyle=':')
    ax[0].set_ylim(800, 830)
    #ax[0].set_xlabel('Pixel indices (#)', fontsize=12, fontweight='bold')
    ax[0].set_ylabel('Serial Overclock Signal (DN)', fontsize=17, fontweight='bold')
    legend = ax[0].legend(loc='lower right', ncol=1, shadow=True,
                          prop={'size':16}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    ax[0].text(200, 802, text1, fontsize=18, fontweight='bold')
    ax[0].tick_params(axis='both', which='major', labelsize=14)
    ax[1].plot(odd_samples_mean_tsoc, '.', color='red', label='Trailing Overclocks')
    ax[1].plot(odd_samples_mean, '.', color='blue', label='Leading Leading')
    ax[1].set_title(title+' (Odd Lines)', fontsize=18, fontweight='bold')
    ax[1].grid(True, linestyle=':')
    ax[1].set_xlabel('Pixel indices (#)', fontsize=17, fontweight='bold')
    ax[1].set_ylabel('Serial Overclock Signal (DN)', fontsize=17, fontweight='bold')
    legend = ax[1].legend(loc='lower right', ncol=1, shadow=True,
                          prop={'size':17}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)
    ax[1].set_ylim(820, 850)
    ax[1].text(200, 822, text2, fontsize=18, fontweight='bold')
    ax[1].tick_params(axis='both', which='major', labelsize=14)
    plt.savefig(figure_name, dpi=200, bbox_inches="tight")
    #plt.show()
    plt.close('all')


    #cc


def main():
    """
    Tme main function
    """
    #nx_quad = 1056 # For Tempo
    #ny_quad = 1046 # For Tempo
    #nlat = nx_quad*2
    #nspec = ny_quad*2




    file_path1 = r'F:\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\Saved_quads'
    save_dir = r'C:\Users\nmishra\Workspace\TEMPO\Leading_Serial_overclocks\Integration_time'
    all_int_files = [each for each in os.listdir(file_path1) \
                         if each.endswith('.dat.sav')]
    if 'Integration_Sweep' in file_path1:
        saturated_collects = ['FT6_LONG_INT_130018.dat.sav', 'FT6_SHORT_INT_0.dat.sav',
                              'FT6_LONG_INT_134990.dat.sav',
                              'FT6_LONG_INT_139961.dat.sav', 'FT6_LONG_INT_145028.dat.sav',
                              'FT6_LONG_INT_149999.dat.sav', 'FT6_LONG_INT_154970.dat.sav',
                              'FT6_LONG_INT_160037.dat.sav', 'FT6_LONG_INT_165008.dat.sav',
                              'FT6_LONG_INT_169980.dat.sav', 'FT6_LONG_INT_175047.dat.sav',
                              'FT6_LONG_INT_180018.dat.sav', 'FT6_LONG_INT_184989.dat.sav',
                              'FT6_LONG_INT_189960.dat.sav', 'FT6_LONG_INT_195027.dat.sav',
                              'FT6_LONG_INT_199999.dat.sav']
    elif 'Intensity_Sweep' in file_path1:
        saturated_collects = ['162_OP_INT_118000.dat.sav', '164_OP_INT_118000.dat.sav',
                              '166_OP_INT_118000.dat.sav', '168_OP_INT_118000.dat.sav',
                              '170_OP_INT_118000.dat.sav', '172_OP_INT_118000.dat.sav',
                              '174_OP_INT_118000.dat.sav', '176_OP_INT_118000.dat.sav',
                              '178_OP_INT_118000.dat.sav', '180_OP_INT_118000.dat.sav',
                              '182_OP_INT_118000.dat.sav', '184_OP_INT_118000.dat.sav',
                              '186_OP_INT_118000.dat.sav', '188_OP_INT_118000.dat.sav',
                              '190_OP_INT_118000.dat.sav', '192_OP_INT_118000.dat.sav',
                              '194_OP_INT_118000.dat.sav', '196_OP_INT_118000.dat.sav',
                              '198_OP_INT_118000.dat.sav', '200_OP_INT_118000.dat.sav',
                              '202_OP_INT_118000.dat.sav']


    nominal_int_files = [items for items in all_int_files
                         if not items.endswith(tuple(saturated_collects))
                         if items in all_int_files]

    for i in range(0, 4):

        #dframe1 = pd.DataFrame()
        #dframe2 = pd.DataFrame()
        #dframe3 = pd.DataFrame()
        all_int_time = []
        active_quad_even_all = []
        active_quad_odd_all = []
        #active_quad_even_all_outlier_filt = []
        #active_quad_odd_all_outlier_filt = []
        tsoc_even_all = []
        tsoc_odd_all = []
        #tsoc_even_all_outlier_filt = []
        #tsoc_odd_all_outlier_filt = []
        #unct_spectral_even = []
        #unct_spectral_odd = []
        #nominal_int_files = [nominal_int_files[0], nominal_int_files[1]]

        for data_files in nominal_int_files:
            data_path_name_split = data_files.split('_')
            data_file = os.path.join(file_path1, data_files)
            print(data_file)

            IDL_variable = readsav(data_file)
            if 'Intensity_Sweep' in file_path1:
                int_time = data_path_name_split[0]
                string1 = 'VA_'
                string2 = 'VA Setting = '
            else:
                int_time = round(int(data_path_name_split[-1].split('.')[0]))
                string1 = 'Integ_time_'
                string2 = 'Int.time = '
            #print(int_time)
            all_int_time.append(int_time)

            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q
            quads = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
            all_full_frame = IDL_variable.q
            quad = all_full_frame[:, i, :, :]
            tsoc_all = quad[:, 4:1028, 0:10]
            trailing_tsoc = quad[:, 4:1028, 1034:1056]

            active_quad = np.mean(quad[:, 4:1028, 10:1034], axis=0)
            #active_quad, smear = perform_smear_subtraction(active_quad, int_time)
            active_quad[active_quad == 16383] = 'nan'
            active_quad_even = active_quad[:, ::2]
            active_quad_even = np.nanmean(active_quad_even, axis=1)
            active_quad_odd = np.nanmean((active_quad[:, 1::2]), axis=1)
            #active_quad_even_outlier_filt = np.mean((active_quad[:, ::2]), axis=1)
            #active_quad_odd_outlier_filt = np.mean((active_quad[:, 1::2]), axis=1)
            active_quad_even_all.append(active_quad_even)
            active_quad_odd_all.append(active_quad_odd)
           # active_quad_even_all_outlier_filt.append(active_quad_even_outlier_filt)
            #active_quad_odd_all_outlier_filt.append(active_quad_odd_outlier_filt)

            quad_dir = quads[i]

            #----------------------------------------------------------------#
            # Ok, let's plot the histogram of saved quads

            all_frames_hist = 'all_frames_hist'
            save_dir_image = os.path.join(save_dir, quad_dir, all_frames_hist)
            if not os.path.exists(save_dir_image):
                os.makedirs(save_dir_image)

            # separate out even and odd lines
            even_samples_all = tsoc_all[:, :, ::2]
            odd_samples_all = tsoc_all[:, :, 1::2]

            even_samples_trailing = trailing_tsoc[:, :, ::2]
            odd_samples_trailing = trailing_tsoc[:, :, 1::2]


            even_samples_avg = np.mean(even_samples_all, axis=0)
            odd_samples_avg = np.mean(odd_samples_all, axis=0)

            even_samples_avg_trailing = np.mean(even_samples_trailing, axis=0)
            odd_samples_avg_trailing = np.mean(odd_samples_trailing, axis=0)


#            tsoc_even_all.append(np.mean((even_samples_avg)))
#            tsoc_even_all_outlier_filt.append(np.mean(filter_outlier_median(even_samples_avg)))
#            tsoc_odd_all.append(np.mean((odd_samples_avg)))
#            tsoc_odd_all_outlier_filt.append(np.mean(filter_outlier_median(odd_samples_avg)))

            title = 'Histogram of Serial Overclocks (All 100 Frames)\n '\
                                                    + quads[i]+', ' + string2 + str(int_time)#+ r" $\mu$" +'secs'
            figure_name = save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
            #plot_hist(even_samples_all, odd_samples_all, title, figure_name)

            avg_frames_hist = 'avg_frames_hist'
            save_dir_image = os.path.join(save_dir, quad_dir, avg_frames_hist)
            if not os.path.exists(save_dir_image):
                os.makedirs(save_dir_image)
            title = 'Histogram of Serial Overclocks (Avg. of 100 Frames)\n '+ quads[i]+', ' + string2 + str(int_time)#+ r" $\mu$" +'secs'
            figure_name = save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
            #plot_hist(even_samples_avg, odd_samples_avg,  title, figure_name)

            tsoc_profile = 'tsoc_plot_outliers_filtered'
            save_tsoc_profile = os.path.join(save_dir, quad_dir, tsoc_profile)
            if not os.path.exists(save_tsoc_profile):
                os.makedirs(save_tsoc_profile)
            figure_name = save_tsoc_profile + '/'+ string1 + str(int_time) + '_image.png'
            title = 'Profile of Serial Overclocks '+ quads[i]+', ' + string2 + str(int_time)#+ r" $\mu$" +'secs'
            even_samples_mean, odd_samples_mean = plot_few_tsocs(even_samples_avg, odd_samples_avg, figure_name, title)
            # save average in spatial direction

            comp_plot = 'leading_trailing_comparisons'
            save_dir_image = os.path.join(save_dir, quad_dir, comp_plot)
            if not os.path.exists(save_dir_image):
                os.makedirs(save_dir_image)
            even_samples_mean_tsoc, odd_samples_mean_tsoc = filter_outlier_trailing(even_samples_avg_trailing, odd_samples_avg_trailing)

            title = 'Comparison between Trailing and Leading Serial Overclocks\n '+ quads[i]+', ' + string2 + str(int_time)+ r" $\mu$" +'secs'
            figure_name = save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
            plot_comparison_overclocks(even_samples_mean_tsoc, odd_samples_mean_tsoc,
                                        even_samples_mean, odd_samples_mean,
                                        title, figure_name )

            #title = 'Histogram of Serial Overclocks (Even and Odd lines)\n '+ quads[i]+', ' + string2 + str(int_time)#+ r" $\mu$" +'secs'
            #figure_name = save_dir_image + '/'+ string1 + str(int_time) + '_image.png'
            #plot_hist(even_samples_mean, odd_samples_mean,  title, figure_name)

            tsoc_even_all.append(even_samples_mean)
            tsoc_odd_all.append(odd_samples_mean)



            #odd_even_lines = [even_samples_avg, odd_samples_avg]
            #odd_even_lines_name = ['Even Lines', 'Odd Lines']
#            for k in np.arange(0, 2):
#                # Lets make directory of
#                #k=1
#                median_plot = r'Outlier_median_tsoc'+'/'+ odd_even_lines_name[k]
#                folder_name_hist = 'Hist_plot'
#                folder_name_mask = 'Mask_plot'
#                folder_name_sat = 'Saturation_plot'
#                save_median_hist = os.path.join(save_dir, quad_dir, median_plot, folder_name_hist)
#                if not os.path.exists(save_median_hist):
#                    os.makedirs(save_median_hist)
#                save_median_mask = os.path.join(save_dir, quad_dir, median_plot, folder_name_mask)
#                if not os.path.exists(save_median_mask):
#                    os.makedirs(save_median_mask)
#                save_sat_mask = os.path.join(save_dir, quad_dir, median_plot, folder_name_sat)
#                if not os.path.exists(save_sat_mask):
#                    os.makedirs(save_sat_mask)
#
#
#                figure_name = save_sat_mask + '/'+ string1 + str(int_time) + '_image.png'
#                sat_pixels = identify_saturation(odd_even_lines[k], int_time, figure_name)
#                outlier_filt_med, outlier_med = reject_outlier_median(odd_even_lines[k])
#                if len(outlier_filt_med) == 1:
#                    outlier_med = 0
#                else:
#                    outlier_med = outlier_med[1]
#
#                title = 'Binary Outlier Mask of Trailing Overclocks (' + quads[i]+', '+ odd_even_lines_name[k]+', '+ string2+str(int_time)#+' micro secs)'
#
#                figure_name = save_median_mask + '/'+ string1 + str(int_time) + '_image.png'
#                median_mask = create_outlier_mask(odd_even_lines[k], outlier_med,
#                                                  title, figure_name)
#
#                title = 'Histogram of Trailing Overclocks (' + quads[i]+', '+ odd_even_lines_name[k]+', '+ string2+ str(int_time)+')'
#                title1 = 'outliers = '+ str(outlier_med.shape[0])
#                figure_name = save_median_hist + '/'+ string1 + str(int_time) + '_image.png'
#                xlim_low = [800, 825]
#                xlim_high = [815, 840]

                #plot_hist_image(odd_even_lines[k], title, title1, outlier_filt_med, outlier_med, figure_name, xlim_low[k], xlim_high[k])

#        print(all_int_time)
#        print(active_quad_even_all_outlier_filt)
#        print(active_quad_odd_all_outlier_filt)
#        print(tsoc_even_all_outlier_filt)
#        print(tsoc_odd_all_outlier_filt)
#
        #dframe1 = pd.DataFrame(
#                  {'Avg_Active_Quad_even' : active_quad_even_all,
#                   'Avg_Active_Quad_odd' : active_quad_odd_all,
#                   'Avg_tsoc_Quad_even' : tsoc_even_all,
#                   'Avg_tsoc_Quad_odd': tsoc_odd_all
#                  })
#        dframe2 = pd.DataFrame(
#                  {'Int_time.' : all_int_time,
#                   'Avg_Active_Quad_even' : active_quad_even_all_outlier_filt,
#                   'Avg_Active_Quad_odd' : active_quad_odd_all_outlier_filt,
#                   'Avg_tsoc_Quad_even' : tsoc_even_all_outlier_filt,
#                   'Avg_tsoc_Quad_odd': tsoc_odd_all_outlier_filt
#                    })
#        dframe3 = pd.DataFrame(
#                  {'Int_time.' : all_int_time,
#                   'Unct_tsoc_Quad_even' : unct_spectral_even,
#                   'Unct_tsoc_Quad_odd': unct_spectral_odd
#                    })
        #data_to_be_saved = np.concatenate((active_quad_even_all, active_quad_odd_all,
                                          # tsoc_even_all, tsoc_odd_all), axis=0)
        csv_save_dir = os.path.join(save_dir, quad_dir)
        if not os.path.exists( csv_save_dir):
            os.makedirs( csv_save_dir)
        #csv_file_name1 = csv_save_dir +'/'+ quads[i]+'_Signal_dependent_offset_analysis.csv'
        #np.savetxt(csv_file_name1, np.array(data_to_be_saved).T, delimiter=',', fmt='%1.2f')
        #csv_file_name2 = quads[i]+'_Signal_dependent_offset_outlier_filt.csv'
        #csv_file_name3 = quads[i]+'_Signal_dependent_offset_unct.csv'
        #dframe1.to_csv(csv_save_dir+'/'+csv_file_name1)
        #dframe2.to_csv(csv_save_dir+'/'+csv_file_name2)
        #dframe3.to_csv(csv_save_dir+'/'+csv_file_name3)

if __name__ == "__main__":
    main()
    
