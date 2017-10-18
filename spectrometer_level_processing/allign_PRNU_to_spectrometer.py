# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:22:32 2017

@author: nmishra
"""
import numpy as np
PRNU_A = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads\Quad A\Final_PRNU\Filtered_PRNU\Quad A_Check_Small_Final_PRNU.csv', delimiter=',')
PRNU_B = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads\Quad B\Final_PRNU\Filtered_PRNU\Quad B_Check_Small_Final_PRNU.csv', delimiter=',')
PRNU_C = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads\Quad C\Final_PRNU\Filtered_PRNU\Quad C_Check_Small_Final_PRNU.csv', delimiter=',')
PRNU_D = np.genfromtxt(r'C:\Users\nmishra\Workspace\TEMPO\Data\GroundTest\FPS\Integration_Sweep\Light\saved_quads\PRNU_quads\Quad D\Final_PRNU\Filtered_PRNU\Quad D_Check_Small_Final_PRNU.csv', delimiter=',')

# now make quad to full_frame spectrometer

PRNU_A_all = PRNU_B_all = PRNU_C_all = PRNU_D_all = np.ones((1024, 1024))
PRNU_A_all[9:1000, :] = PRNU_A
PRNU_B_all[9:1000, :] = PRNU_B
PRNU_C_all[9:1000, :] = PRNU_C
PRNU_D_all[9:1000, :] = PRNU_D
          
# At spectrometer level size of active quad is 1028, 1023

#'s backfill again

PRNU_A_all_s = PRNU_B_all_s = PRNU_C_all_s = PRNU_D_all_s = np.ones((1028, 1023))

PRNU_A_all_s[0:1024, :] = PRNU_A_all[:, 0:1023]
PRNU_B_all_s[0:1024, :] = PRNU_B_all[:, 0:1023]
PRNU_C_all_s[0:1024, :] = PRNU_C_all[:, 0:1023]
PRNU_D_all_s[0:1024, :] = PRNU_D_all[:, 0:1023]
            
lower_quad = np.concatenate((PRNU_D_all_s, np.fliplr(PRNU_C_all_s)), axis=1)
upper_quad = np.concatenate((np.flipud(PRNU_D_all_s), np.rot90(PRNU_C_all_s,2)), axis=1)
PRNU_full_frame_spectrometer = np.concatenate((lower_quad, upper_quad), axis=0)
print(PRNU_full_frame_spectrometer.shape)
final_mask_name = r'C:\Users\nmishra\Workspace\TEMPO\PRNU_map\Final_Map\PRNU_mask_active_regions.csv'
np.savetxt(final_mask_name, np.array(PRNU_full_frame_spectrometer), delimiter=",")
            