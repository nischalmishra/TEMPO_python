# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:49:55 2018

@author: nmishra
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    Main Function
    """
    data_dir = r'F:\TEMPO\Data\GroundTest\FPS\Spectrometer\Photon_Transfer_TVAC\saved_quads\saved_plots\SDO\updated'
    file_name_active_even = os.path.join(data_dir, ' quad_D_active_odd.csv')
    file_name_tsoc_even = os.path.join(data_dir, ' quad_D_tsoc_odd.csv')
    file_name_active_odd = os.path.join(data_dir, ' quad_D_active_even.csv')
    file_name_tsoc_odd = os.path.join(data_dir, ' quad_D_tsoc_even.csv')
    int_file = os.path.join(data_dir, 'integration_time.csv')

    quad_name1 = 'Quad A (Visible Channel, Sp. Index 1:1024)'
    quad_name2 = 'Quad B (Visible Channel, Sp. Index 1025:2048)'
    quad_name3 = 'Quad C (UV Channel, Sp. Index 1025:2048)'
    quad_name4 = 'Quad D (UV Channel, Sp. Index 1:1024)'


    active_data_odd = pd.read_csv(file_name_active_odd, delimiter=",")
    integration_time = pd.read_csv(int_file, delimiter=",")
    tsoc_data_odd = pd.read_csv(file_name_tsoc_odd, delimiter=",")
    integration_time.columns = ['Int. time']
    length = active_data_odd.shape[1]   
    active_data_odd = active_data_odd.values    
    tsoc_data_odd = tsoc_data_odd.values
    #print(tsoc_data_odd)



    active_data_even = pd.read_csv(file_name_active_even, delimiter=",")
    #integration_even = pd.read_csv(int_file, delimiter=",")
    tsoc_data_even = pd.read_csv(file_name_tsoc_even, delimiter=",")
    active_data_even = active_data_even.values
    tsoc_data_even = tsoc_data_even.values
    #print(tsoc_data_even)
    #
    active_dns_odd = []
    tsoc_dns_odd = []
    active_dns_even = []
    tsoc_dns_even = []

    for i in np.arange(0, length-1):
        #print(i)
        dframe1 = pd.DataFrame()
        dframe2 = pd.DataFrame()
        dframe1['Int_time'] = integration_time['Int. time'].values
        dframe1['Signal_Odd'] = np.array(active_data_odd)[:, i]
        dframe1['Signal_Even'] = np.array(active_data_even)[:, i]
        dframe1 = dframe1.groupby(['Int_time']).mean()
        active_dns_odd.append(dframe1['Signal_Odd'])
        active_dns_even.append(dframe1['Signal_Even'])

        dframe2['Int_time'] = integration_time['Int. time'].values
        dframe2['Signal_Odd'] = np.array(tsoc_data_odd)[:, i]

        dframe2['Signal_Even'] = np.array(tsoc_data_even)[:, i]
        dframe2 = dframe2.groupby(['Int_time']).mean()
        tsoc_dns_odd.append(dframe2['Signal_Odd'])
        tsoc_dns_even.append(dframe2['Signal_Even'])
    saturation_threshold = 0.90*16383
    active_dns_odd = np.array(active_dns_odd)
    tsoc_dns_odd = np.array(tsoc_dns_odd)
    rows, cols = active_dns_odd.shape
    active_dns_odd = np.reshape(active_dns_odd, [rows*cols, 1])
    valid_pixels = active_dns_odd < saturation_threshold
    tsoc_dns_odd = np.reshape(tsoc_dns_odd, [rows*cols, 1])
    valid_counts_odd = active_dns_odd[valid_pixels]
    valid_tsoc_odd = tsoc_dns_odd[valid_pixels]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(valid_counts_odd/(0.9*16383), valid_tsoc_odd/(0.9*16383), 'b.', markersize=1, label='Odd Lines')
    ax1.grid(linestyle=':')
    ax1.set_title('Active Area Signal Vs TSOC '+ quad_name4)
    ax1.set_xlabel(' Active Area Signal (Fraction of Full Well)')
    ax1.set_ylabel('TSOC Counts (Fraction of Full Well)')
    ax1.legend(loc='best')



    active_dns_even = np.reshape(active_dns_even, [rows*cols, 1])
    tsoc_dns_even = np.reshape(tsoc_dns_even, [rows*cols, 1])
    valid_counts_even = active_dns_even[valid_pixels]
    valid_tsoc_even = tsoc_dns_even[valid_pixels]
    #fig = plt.figure()
    ax2 = fig.add_subplot(212)
    ax2.plot(valid_counts_even, valid_tsoc_even, 'r.', markersize=1, label='Even Lines')
    ax2.grid(linestyle=':')
    #ax2.set_title('Active Area Signal Vs TSOC '+ quad_name1)
    ax2.set_xlabel(' Active Area Signal (Fraction of Full Well)')
    ax2.set_ylabel('TSOC Counts (Fraction of Full Well)')
    ax2.legend(loc='best')

    plt.show()
    plt.close('all')

#    print('QuadA')
#    file_name_active = os.path.join(data_dir, ' quad_A_active_odd.csv')
#    file_name_tsoc = os.path.join(data_dir, ' quad_A_tsoc_odd.csv')
#    int_file = os.path.join(data_dir, 'integration_time.csv')
#    active_data = pd.read_csv(file_name_active, delimiter=",")
#    integration_time = pd.read_csv(int_file, delimiter=",")
#    tsoc_data = pd.read_csv(file_name_tsoc, delimiter=",")
#    integration_time.columns = ['Int. time']
#    length = active_data.shape[1]
#    active_data = active_data.values
#    tsoc_data = tsoc_data.values
#    active_dns = []
#    tsoc_dns = []
#
#
#
#    for i in np.arange(0, length-1):
#        dframe1 = pd.DataFrame()
#        dframe2 = pd.DataFrame()
#        dframe1['Int_time'] = integration_time['Int. time'].values
#        dframe1['Signal'] = np.array(active_data)[:, i]
#        dframe1 = dframe1.groupby(['Int_time']).mean()
#        #dframe1.reset_index
#        active_dns.append(dframe1['Signal'])
#        dframe2['Int_time'] = integration_time['Int. time'].values
#        dframe2['Signal'] = np.array(tsoc_data)[:, i]
#        dframe2 = dframe2.groupby(['Int_time']).mean()
#        tsoc_dns.append(dframe2['Signal'])
#    saturation_threshold = 0.90*16383
#    active_dns = np.array(active_dns)
#    tsoc_dns = np.array(tsoc_dns)
#    #print(active_dns.shape)
#    #print(tsoc_dns.shape)
#    rows, cols = active_dns.shape
#    #print(rows*cols)
#    active_dns = np.reshape(active_dns, [rows*cols, 1])
#    valid_pixels = active_dns < saturation_threshold
#    tsoc_dns = np.reshape(tsoc_dns, [rows*cols, 1])
#    valid_counts = active_dns[valid_pixels]
#    valid_tsoc = tsoc_dns[valid_pixels]
#    fig = plt.figure()
#    ax2 = fig.add_subplot(211)
#    ax2.plot(valid_counts/(0.9*2**14), valid_tsoc/(0.9*16383), 'r.', markersize=1, label='Even Lines')
#    ax2.grid(linestyle=':')
#    ax2.set_title('Active Area Signal Vs. TSOC '+ quad_name2)
#    ax2.set_xlabel(' Active Area Signal (Fraction of Full Well)')
#    ax2.set_ylabel('Trailing Overclock Counts(Fraction of Full Well)')
#    ax2.legend(loc='best')
#    ax2.set_ylim(956, 976, 2)
#    plt.show()


if __name__ == "__main__":
    main()
