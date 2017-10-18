# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:30:34 2017

@author: nmishra
"""

import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns
def main():
    lower_quads = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads\Outlier_median\Lower_quads\Mask_plot'
    upper_quads = r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads\Outlier_median\Upper_quads\Mask_plot'
    mat_file_dir = 'mat_files_saved'
    
    lower_outliers = os.path.join(lower_quads,mat_file_dir)
    upper_outliers = os.path.join(upper_quads,mat_file_dir)
    
    data_path_lower = [each for each in os.listdir(lower_outliers) if each.endswith('.mat')]
    data_path_outer = [each for each in os.listdir(upper_outliers) if each.endswith('.mat')]
    print(len(data_path_lower))
    
    outliers_lower_quads =  [[0]*2048]*1024
    k=1
    for path in data_path_lower:
        print (os.path.join(lower_quads,mat_file_dir, path))
        
        read_outlier_mask_lower = sio.loadmat(os.path.join(lower_quads,mat_file_dir, path))
         
        read_outlier_mask_lower = read_outlier_mask_lower['mask']
        #print(np.where(read_outlier_mask_lower==1))
        #lower_outliers = outliers_lower_quads| read_outlier_mask_lower
        plt.figure(frameon=False)
        plt.grid(False)
        plt.imshow(np.invert(read_outlier_mask_lower), cmap='bwr', interpolation='none', origin='lower')
            
        final_outliers_lower = np.array(np.where(np.reshape(read_outlier_mask_lower,
                                          (1024*2048, 1))==1)).shape[1]
        plt.title('Quad A & B final outlier mask  (outliers = '+ str(final_outliers_lower)+')',
              fontsize=14)
        plt.savefig(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads\Outlier_median'+'/'+str(path)+'.png', dpi=100,
                bbox_inches="tight")
        k = k+1
    
    cc
    plt.figure(frameon=False)
    plt.imshow(lower_outliers, cmap='binary', interpolation='none', origin='lower')
        
    final_outliers_lower = np.array(np.where(np.reshape(lower_outliers,
                                          (1024*2048, 1))==1)).shape[1]
        
    plt.title('Quad A & B final outlier mask  (outliers = '+ str(final_outliers_lower)+')',
              fontsize=14)
    plt.xlabel('# of spatial pixels', fontsize=12)
    plt.ylabel('# of spectral pixels', fontsize=12)
    plt.grid(False)
    sns.despine(left=False, bottom=True, right=True)
    plt.savefig(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Dark\FPS_Dark\saved_quads\Outlier_median\save_masksfinal_mask_lower.png', dpi=100,
                bbox_inches="tight")
if __name__ == "__main__":
    main()