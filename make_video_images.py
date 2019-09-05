# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:13:22 2018

@author: nmishra
"""

import cv2
import os

image_folder = r'C:\Users\nmishra\Workspace\TEMPO_Spectrometer\MTF_analysis\390_nm_px\NS\pixel_response'
video_name = 'video.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
#print(images)
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*"MJPG"), 1, (width,height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))
    #video.release()
#xx

cv2.destroyAllWindows()
video.release()