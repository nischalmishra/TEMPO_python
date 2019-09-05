# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 10:58:23 2017

@author: nmishra
"""
import os
import random
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

#from scipy import stats as sc
import seaborn as sns
#import matplotlib.mlab as mlab

#******************************************************************************

def plot_full_frame_image(full_frame, image_title, collection_type, frames,
                          int_time, plot_dir):
    """
    Function to make a plot frame image including the overclock pixels.
    The purpose is to check if the quads are alligned properly or not.
    TO DO: Once the path to save plots is known beforehand, pass the location
    through the main function
    """
    int_time = str(int_time)+' micro secs' # for integ_sweep
    #int_time = str(int_time) # for instensity Sweep
    title = image_title + ' ('+ collection_type + ' Integration, ' + int_time+')'
    #title = image_title + ' ('+ collection_type + ' Intensity, ' + int_time+')'
    plt.figure()
    image = plt.imshow(full_frame, cmap='bwr', interpolation='none',
                       origin='lower')
    cbar = plt.colorbar(image)
    plt.title(title, fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    plt.savefig(plot_dir+'/'+ collection_type+'_'+\
                    frames+'.png', dpi=100, bbox_inches="tight")
    plt.close('all')
#******************************************************************************

def  plot_each_quad(quads, image_title, collection_type, frames, plot_dir):
    """
    Makes image of each quad
    """
    title = image_title + ' ('+ collection_type + ' Integration)'
    plt.figure()
    plt.grid(False)
    plt.subplot(221), plt.imshow(quads[3], cmap='gray')
    plt.colorbar()
    plt.title('Quad D'), plt.xticks([]), plt.yticks([])
    plt.subplot(222), plt.imshow(quads[2], cmap='gray')
    plt.colorbar()
    plt.title('Quad C'), plt.xticks([]), plt.yticks([])
    plt.subplot(223), plt.imshow(quads[0], cmap='gray')
    plt.colorbar()
    plt.title('Quad A'), plt.xticks([]), plt.yticks([])
    plt.subplot(224), plt.imshow(quads[1], cmap='gray')
    plt.colorbar()
    plt.title('Quad B'), plt.xticks([]), plt.yticks([])
    plt.suptitle(title)
    plt.savefig(plot_dir+'/'+ collection_type+'_'+\
                    frames+'.png', dpi=100, bbox_inches="tight")

    plt.close('all')
#******************************************************************************
def calculate_std_dev(full_frame):
    """
    This function calculated the standard deviation in spatial and
    spectral direction including the overclocks and plots them in a single
    figure usinf dual y-axis plot

    Input : Full_frame image
    Output : spectral and spatial standard deviation expressed as Coefficient
    of Variation, in a 2D- matrix

    To DO : Remove overclock pixels (both in x and y direction)

    """
    std_spectral = (np.std(full_frame, axis=0)/np.mean(full_frame, axis=0)*100)
    std_spatial = (np.std(full_frame, axis=1)/np.mean(full_frame, axis=1)*100)
    std_spectral = std_spectral.round(2)
    std_spatial = std_spatial.round(2)
    #plt.figure()
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(std_spectral, 'g-', label='spectral direction')
    ax2.plot(std_spatial, 'b-', label='spatial direction')
    ax1.set_xlabel('pixel#')
    ax1.set_xticks(np.arange(0, 2400, 300))
    ax1.set_yticks(np.arange(0, 21, 4))
    ax2.set_yticks(np.arange(0, 21, 4))
    ax1.set_ylabel('% Uncertainty [1-Sigma stdev/Mean)*100] ', color='g',
                   fontweight='bold')
    ax2.set_ylabel('% Uncertainty [1-Sigma stdev/Mean)*100]', color='b',
                   fontweight='bold')
    ax1.set_title('Coefficient of Variation along spectral & spatial direction',
                  fontsize=12, fontweight='bold')
    ax1.legend(loc=2)
    ax2.legend(loc=0)
    ax1.set_xlim(xmin=0)
    ax1.set_ylim(ymin=0)
    ax2.set_ylim(ymin=0)
    ax1.grid(True, linestyle=':')
    plt.show()
    #plt.close('all')
    return (std_spectral, std_spatial)

#*****************************************************************************
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
    plt.savefig(r"C:\Users\nmishra\Desktop\test\full_frame_fft.png",
                bbox_inches="tight")

#******************************************************************************
def calculate_fft_each_quad(quads, nx_quad, ny_quad):
    """
    This function calculates the FFT of each quads in spatial and spectral
    direction. The spatial direction FFT and spatial direction FFTS are plotted
    in two figures. Each figure contains 4 plots of each quads. The motivation
    is to check if similar frequency components exist in each quad.

    Note: NumPy's Fourier coefficients are N times larger than expected.
    Follow FFT equation to check what it means

    Input : Full_frame image, quads lengths in each direction
    Output : Image plot, magnitude and phase plot of the full frame image.
    TO DO : Remove the hard coded values such as plot directory. In future the
    paths will be provided from the main function
    """
    #quad_names = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
    #color = ['blue', 'green', 'red', 'purple']
    #nrows = 2 # to create 2*2 plot for each quads in each direction
    #ncols = 2
    k = 0
    # Number fo sample points
    #N = nx_quad * ny_quad
    #fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 12))
    #fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                       # wspace=0.3, hspace=.25)
    i = 0
    for quad in quads:
        #ax = axes[int(k / ncols)][int(k % ncols)]
        quadi = np.asarray(quads[k])
        fft_s = np.fft.fft2(quadi)
        signal_PSD = np.abs(fft_s)**2
        freqx = np.fft.fftfreq(signal_PSD.shape[0])
        freqy = np.fft.fftfreq(signal_PSD.shape[1])
        i = freqx > 0
        j = freqy > 0
        plt.figure()
        plt.plot(freqx[i], 20*np.log10(signal_PSD[i]))
        plt.show()
        plt.figure()
        plt.plot(freqy[j], 20*np.log10(signal_PSD[j]))
        plt.show()

def moving_average_filter(data, size):
    """
    Computes the moving average on a 1-D data.
    """
    # Note moving averaging can be done either using this cumulative sum
    #approach or convolution approach. Convolution approach uses FFT
    #and hence the speed decreases as the the length of data increases.
    #http://stackoverflow.com/questions/13728392/moving-average-or-running-mean/27681394#27681394
    cum_sum = np.cumsum(np.insert(data, 0, 0)) # python indexing starts at zero
    running_mean = (cum_sum[size:] - cum_sum[:-size]) / size
    return running_mean

def median_filtering(data, window):
    """
    Apply a length(defined by window) median filter to 1-D array  data.
    Also remember the window size or the kernel size must be odd length
    """
    assert window % 2 == 1 # Median filter length must be odd
    assert data.ndim == 1  # Input must be one dimensional
    med_smoothed_data = sp.signal.medfilt(data, window)
    return med_smoothed_data

def random_num_generator(length):
    """ Generates random number from 0 to length
    """
    return random.sample(range(0, length), 3) # 3 refers to # of numbers


#******************************************************************************
def plot_hist_each_quad(quads, title, collection_type, frames):
    """
    This function calculates the normalized histogram or pdf of each quad with
    bin size of 200. Since the quads contain the overclock bits too, the x-axis
    on the histogram has been stretched to see the true histogram.

    Input: full_frame image, nx_quad and ny_quad
    Output: Histogram plot of each qyads.
    TO DO : Remove the hard coded values such as plot directory. In future the
    paths will be provided from the main function
    """
    title = title + ' ('+ collection_type + ' Integration)'
    nx_quad = quads.shape[1]
    print(nx_quad)
    ny_quad = quads.shape[2]
    print(ny_quad)
    quad_names = ['Quad A', 'Quad B', 'Quad C', 'Quad D']
    color = ['blue', 'green', 'red', 'purple']
    nrows = 2 # to create 2*2 plot for each quads in each direction
    ncols = 2
    k = 0
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)
    for quad in quads:
        ax = axes[int(k / ncols)][int(k % ncols)]
        #(mean, sigma) = sc.fit(quads[k])
        mean = np.mean(quads[k])
        sigma = np.std(quads[k])
        med = np.median(quads[k])
        max_val = np.max(quads[k])
        min_val = np.min(quads[k])
        label = 'Mean = '+ str(round(mean, 1)) + \
                '\n Median = '+ str(round(med, 2)) + \
                '\n Std. = '+ str(round(sigma, 1))+ \
                '\n Max = '+ str(round(max_val, 2)) + \
                '\n Min = '+ str(round(min_val, 1))
        sns.set_context("talk")
        with sns.axes_style("darkgrid"):
            ax.hist(np.reshape(quads[k], (nx_quad*ny_quad, 1)),
                    200, normed=0, facecolor=color[k], alpha=0.75,
                    label=label)

          # Add the best fit line
        #y = mlab.normpdf( bins, mu, sigma)
        #l = ax.plot(bins, y, '--o',color='k',linewidth=0,markersize=0)
            ax.tick_params(axis='x', pad=10)
            ax.grid(True, linestyle=':')
            legend = ax.legend(loc='best', ncol=3, shadow=True,
                               prop={'size':12}, numpoints=1)
            legend.get_frame().set_edgecolor('r')
            legend.get_frame().set_linewidth(2.0)
            ax.set_ylabel('Frequency (# of pixels)', fontsize=15,
                          fontweight="bold")
            #ax.set_xlim(10000, 14000)
            ax.set_xlabel('Counts (DNs)', fontsize=14, fontweight="bold")
            ax.set_title(str(quad_names[k]), fontsize=14, fontweight="bold")
            k += 1
    for del_l in range(k, ncols* nrows):
        fig.delaxes(axes[del_l / ncols][del_l % ncols])
    plt.suptitle(title, fontsize=18, fontweight='bold')
    figure_name = r'C:\Users\nmishra\Desktop\test\Quad_based'+'/'+\
                    collection_type+'_'+\
                    frames+'.png'
    fig.savefig(figure_name, dpi=100)
    plt.close('all')



def plot_hist_image(quads, tit, collection_type, frames,
                    outlier_filt_data, outlier_dets, int_time, plot_dir):
    """
    This function calculates the normalized histogram or pdf of each quad with
    bin size of 200. Since the quads contain the overclock bits too, the x-axis
    on the histogram has been stretched to see the true histogram.

    Input: full_frame image, nx_quad and ny_quad
    Output: Histogram plot of each qyads.
    TO DO : Remove the hard coded values such as plot directory. In future the
    paths will be provided from the main function
    """
    int_time = str(int_time)+' micro secs'
    int_time = str(int_time)
    title = tit + ' ('+ collection_type + ', Int. time = ' + int_time+')'
    #title = tit + ' ('+ collection_type + ' Intensity, ' + int_time+')' # for intensity sweep

    #title = tit + ' ('+ collection_type + ' Integration)'
    num_outlier = 'outliers = '+ str((np.array(outlier_dets)).shape[1])

    title1 = tit + ' ('+ collection_type + ', Int. time = ' + int_time+ ', ' \
             +num_outlier+')'
    #title1 = tit + ' ('+ collection_type + ' Intensity' +', ' +num_outlier+')' #for instesity sweep

    nx_quad = quads.shape[0]
    #print(nx_quad)
    ny_quad = quads.shape[1]
   # print(ny_quad)

    nrows = 2 # to create 2*2 plot for each quads in each direction
    ncols = 1

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
    fig.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.3, hspace=.25)

    #(mean, sigma) = sc.fit(quads[k])
    mean = np.mean(quads)
    sigma = np.std(quads)
    med = np.median(quads)
    max_val = np.max(quads)
    min_val = np.min(quads)
    label = 'Mean = '+ str(round(mean, 1)) + \
            '\n Median = '+ str(round(med, 2)) + \
            '\n Std. = '+ str(round(sigma, 2))+ \
            '\n Max = '+ str(round(max_val, 2)) + \
            '\n Min = '+ str(round(min_val, 1))

    mean1 = np.mean(outlier_filt_data)
    sigma1 = np.std(outlier_filt_data)
    med1 = np.median(outlier_filt_data)
    max_val1 = np.max(outlier_filt_data)
    min_val1 = np.min(outlier_filt_data)
    label1 = 'Mean = '+ str(round(mean1, 1)) + \
            '\n Median = '+ str(round(med1, 2)) + \
            '\n Std. = '+ str(round(sigma1, 2))+ \
            '\n Max = '+ str(round(max_val1, 2)) + \
            '\n Min = '+ str(round(min_val1, 1))



    sns.set_context("talk")
    with sns.axes_style("darkgrid"):

        ax[0].hist(np.reshape(quads, (nx_quad*ny_quad, 1)),
                   50, normed=0, facecolor='red', alpha=0.75,
                   label=label)

        ax[0].tick_params(axis='x', pad=10)
        ax[0].grid(True, linestyle=':')
        legend = ax[0].legend(loc='best', ncol=3, shadow=True,
                              prop={'size':12}, numpoints=1)
        legend.get_frame().set_edgecolor('r')
        legend.get_frame().set_linewidth(2.0)
        ax[0].set_ylabel('Frequency (# of pixels)', fontsize=15,
                         fontweight="bold")
        #ax.set_xlim(10000, 14000)
        ax[0].set_xlabel('Counts (DNs)', fontsize=14, fontweight="bold")
        ax[0].set_title(title, fontsize=14, fontweight="bold")
        ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        # Now for filtered data

        ax[1].hist(outlier_filt_data, 50, normed=0, facecolor='blue',
                   label=label1)

        ax[1].tick_params(axis='x', pad=10)
        ax[1].grid(True, linestyle=':')
        legend = ax[1].legend(loc='best', ncol=3, shadow=True,
                              prop={'size':12}, numpoints=1)
        legend.get_frame().set_edgecolor('r')
        legend.get_frame().set_linewidth(2.0)
        ax[1].set_ylabel('Frequency (After outliers rejection)',
                         fontsize=15, fontweight="bold")
        #ax.set_xlim(10000, 14000)
        ax[1].set_xlabel('Counts (DNs)', fontsize=14, fontweight="bold")
        ax[1].set_title(title1, fontsize=14, fontweight="bold")
        ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        figure_name = plot_dir +'/'+collection_type+'_'+\
                      frames+'.png'
        fig.savefig(figure_name, dpi=100)
        plt.close('all')
