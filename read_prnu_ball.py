# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:18:31 2017

@author: nmishra
"""

def create_image(image_data):
    plt.figure()
    ax = plt.gca()
    cmap = cm.get_cmap('jet')
    image = ax.imshow(image_data, cmap=cmap, origin='lower')    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)            
    plt.colorbar(image, cax= cax)            
    plt.grid(False)    
    #plt.savefig(figure_name,dpi=95,bbox_inches="tight")
    plt.show()
    cc
    plt.close('all')   

import numpy as np
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm

hdf_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\batch_2017Jun20_TEMPO_PRNU_-20Tccd__46Tfpe_3pixSpectral_3pixSpatial.h5'
file = h5py.File(hdf_name, 'r')
prnu = file.get('prnu')
prnu = np.array(prnu).transpose()
print(prnu.shape)
quad_D = prnu[2:1030, 10:1034]
quad_C = prnu[2:1030, 1078:2102]
quad_A = prnu[1062:2090, 10:1034]
quad_B = prnu[1062:2090, 1078:2102]
prnu_map_lower = np.concatenate((quad_D, quad_C), axis=1)
prnu_map_upper = np.concatenate((quad_A, quad_B), axis=1)
prnu_map = np.concatenate((prnu_map_lower, prnu_map_upper), axis=0)
create_image(prnu_map)